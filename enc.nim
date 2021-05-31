import hts/bam
import hts/vcf
import hts/fai
import hts/files
import zippy
import base64
import os
import random
import strutils
import strformat

randomize()

proc encode*(ibam:Bam, region:string): string =

  var obam:Bam
  var path = os.getTempDir() & "/" &  $(rand(int.high)) & ".bam"

  defer:
    discard os.tryRemoveFile(path)

  if not obam.open(path, mode="wb"):
    quit "could not open open bam"

  obam.write_header(ibam.hdr)

  for aln in ibam.query(region):
    obam.write(aln)

  obam.close()

  result = base64.encode(path.readFile)

proc encode*(ivcf:VCF, region:string): string =

  var ovcf:VCF
  var path = os.getTempDir() & "/" &  $(rand(int.high)) & ".vcf.gz"

  defer:
    discard os.tryRemoveFile(path)

  if not ovcf.open(path, mode="wz"):
    quit "could not open open bam"

  let chrom = region.split(":")[0]
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
  result = base64.encode(path.readFile)

proc encode*(fai:Fai, region:string): string =
  var s = fai.get(region)
  var tmp = &">{region}\n{s}"
  result = base64.encode(compress(tmp))

type FileType* {.pure.} = enum
  cytoband

template stripChr*[T:string|cstring](s:T): string =
  if s.len > 3 and ($s).startswith("chr"): ($s)[3..<s.len] else: $s

proc sameChrom(a: string, b: string): bool =
  return stripChr(a) == stripChr(b)

proc encode*(path:string, region:string, typ:FileType): string =

  let chrom = stripChr(region.split(":")[0])

  case typ
  of FileType.cytoband:
    var clines: seq[string]
    for line in path.hts_lines:
      var toks = line.split("\t")
      if not sameChrom(toks[0], chrom): continue
      clines.add(line)
    var tmp = clines.join("\n") & "\n"
    return base64.encode(compress(tmp))

when isMainModule:

  var ibam:Bam
  if not ibam.open("/data/human/hg002.cram", fai="/data/human/g1k_v37_decoy.fa", index=true):
    quit "could not open cram"

  #echo encode(ibam, "1:22000000-22000900")

  var ivcf:VCF
  if not ivcf.open("/data/human/HG002_SVs_Tier1_v0.6.vcf.gz"):
    quit "could not open vcf"

  #echo encode(ivcf, "1:250000-3000000")

  var fa:Fai
  if not fa.open("/data/human/Homo_sapiens_assembly38.fasta"):
    quit "could not open fai"


  #echo fa.encode("chr5:474488-475489")
  #

  echo encode("/home/brentp/src/igv-reports/examples/variants/cytoBandIdeo.txt", "chr5:474969-475009", FileType.cytoband)
