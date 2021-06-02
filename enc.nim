import hts/bam
import hts/vcf
import json
import hts/fai
import hts/bgzf/bgzi
import hts/files
import pedfile
import zippy
import base64
import os
import random
import strutils
import strformat
import ./track

randomize()

const prefix = "data:application/gzip;base64,"

template stripChr*[T:string|cstring](s:T): string =
  if s.len > 3 and ($s).startswith("chr"): ($s)[3..<s.len] else: $s

proc sameChrom(a: string, b: string): bool =
  return stripChr(a) == stripChr(b)

proc check_chrom(targets: seq[Target|Contig], chrom:string): string =
  for t in targets:
    if sameChrom(t.name, chrom):
      return t.name
  raise newException(KeyError, &"[jigv] error chromosome: {chrom} not found in bam")

proc encode*(ibam:Bam, region:string): string =

  var obam:Bam
  var path = os.getTempDir() & "/" &  $(rand(int.high)) & ".bam"

  defer:
    discard os.tryRemoveFile(path)

  if not obam.open(path, mode="wb"):
    quit "could not open open bam"

  obam.write_header(ibam.hdr)

  # handle chr prefix stuff
  var chrom = region.split(':')[0]
  chrom = check_chrom(ibam.hdr.targets, chrom)
  var region = &"{chrom}:{region.split(':')[1]}"

  for aln in ibam.query(region):
    obam.write(aln)

  obam.close()
  result = prefix & base64.encode(path.readFile)

proc encode*(ivcf:VCF, region:string): string =

  var ovcf:VCF
  var path = os.getTempDir() & "/" &  $(rand(int.high)) & ".vcf.gz"

  defer:
    discard os.tryRemoveFile(path)

  if not ovcf.open(path, mode="wz"):
    quit "could not open open bam"

  var contigs = ivcf.contigs
  let chrom = check_chrom(contigs, region.split(":")[0])

  var new_header: seq[string]
  for l in ($(ivcf.header)).split("\n"):
    if l.startswith("##contig=") and &"ID={chrom}," notin l: continue
    new_header.add(l)

  var h:vcf.Header
  h.from_string(new_header.join("\n"))

  # TODO: can drop contig lines from header.
  ovcf.copy_header(h)
  doAssert ovcf.write_header

  for v in ivcf.query(region):
    doAssert ovcf.write_variant(v)

  ovcf.close()
  result = prefix & base64.encode(path.readFile)

proc encode*(s:string): string =
  return prefix & base64.encode(compress(s))

proc expand(region:string, dist:int): string =
  if dist == 0: return region
  var chromse = region.split(":")
  var se = chromse[1].split("-")
  return &"{chromse[0]}:{max(0, parseInt(se[0]) - dist)}-{parseInt(se[1]) + dist}"

proc encode*(fai:Fai, region:string, expand:int=10): string =

  var region = region.expand(expand)

  var s:string
  try:
    s = fai.get(region)
  except:
    s = fai.get(if not region.startswith("chr"): "chr" & region else: stripChr(region))
    stderr.write_line "[jigv] found fasta sequence with chr prefix change"
  var tmp = &">{region}\n{s}"
  result = prefix & base64.encode(compress(tmp))

type TrackFileType* {.pure.} = enum
  cytoband
  bed
  bed12

proc encode*(path:string, region:string, typ:TrackFileType): string =

  let chrom = stripChr(region.split(":")[0])
  let start = parseInt(region.split(":")[1].split("-")[0])
  let stop =  parseInt(region.split(":")[1].split("-")[1])

  case typ
  of TrackFileType.cytoband, TrackFileType.bed:
    var clines: seq[string]
    for line in path.hts_lines:
      var toks = line.split("\t")
      if not sameChrom(toks[0], chrom): continue
      # for bed format, we check the actual positions.
      if typ == TrackFileType.bed:
        let s = parseInt(toks[1])
        if s > stop: continue
        let e = parseInt(toks[2])
        if e < start: continue
      clines.add(line)
    clines.add("") # so we get a trailing newline
    var tmp = clines.join("\n")
    return prefix & base64.encode(compress(tmp, dataFormat=dfGzip, level=1))

  of TrackFileType.bed12:

    var bgz:BGZI
    if not bgz.open(path):
      stderr.write_line &"[jigv] warning: {path} should be compressed and with csi index. trying (slow) full search over lines instead."
      return path.encode(region, TrackFileType.bed)

    var chromstuff = region.split(":")
    var start = parseInt(chromstuff[1].split('-')[0])
    var stop = parseInt(chromstuff[1].split('-')[1])
    var clines:seq[string]
    for l in bgz.query(chromstuff[0], start - 1, stop):
      clines.add(l)
    var tmp = clines.join("\n") & "\n"
    return prefix & base64.encode(compress(tmp))

iterator encode*(variant:Variant, ivcf:VCF, bams:seq[Bam], fasta:Fai, samples:seq[pedfile.Sample], sample_i:int,
                 anno_files:seq[string], note:string="", n_backgrounds:int=0, flank:int=150): tuple[left:string, right:string] =
  let locus = &"{variant.CHROM}:{max(1, variant.start - flank)}-{variant.stop + flank}"
  var json:JsonNode = %* {
      "locus": locus,
      "reference": {"fastaURL": fasta.encode(locus) },
      "queryParametersSupported": true,
      "showChromosomeWidget": false,
    }

  var sample = samples[sample_i]

  var tracks: seq[Track]
  echo "fname:", ivcf.fname
  # TODO: set n_tracks
  var tr = Track(name:extractFileName(ivcf.fname), path:ivcf.fname, n_tracks:2, file_type:FileType.VCF, region:locus)
  var js = % tr
  js["url"] = % ivcf.encode(locus)
  echo js
  yield ($js, "")




when isMainModule:
  import pedfile

  var ibam:Bam
  if not ibam.open("/data/human/hg002.cram", fai="/data/human/g1k_v37_decoy.fa", index=true):
    quit "could not open cram"

  #echo encode(ibam, "1:22000000-22000900")

  var ivcf:VCF
  if not ivcf.open("/data/human/HG002_SVs_Tier1_v0.6.vcf.gz"):
    quit "could not open vcf"

  #echo encode(ivcf, "1:250000-3000000")

  var fa:Fai
  if not fa.open("/data/human/g1k_v37_decoy.fa"):
    quit "could not open fai"


  #echo fa.encode("chr5:474488-475489")
  #

  #echo encode("/home/brentp/src/igv-reports/examples/variants/cytoBandIdeo.txt", "chr5:474969-475009", TrackFileType.cytoband)
  #

  #echo encode("/home/brentp/src/igv-reports/examples/variants/refGene.sort.bed.gz", "chr5:474969-475009", TrackFileType.bed12)

  var ifiles:seq[string]
  for v in ivcf:
    if v.FILTER != "PASS": continue
    for tr in v.encode(ivcf, @[ibam], fa, @[Sample(id:"HG002", i:0)], 0, ifiles):
      discard
      #echo tr

