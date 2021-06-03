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
import tables
import sets
import random
import strutils
import pedfile
import strformat
import ./track
import argparse

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
    v.c.rid = 0
    doAssert ovcf.write_variant(v), $(region, $v)

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

proc get_samples(samples:seq[Sample], sample_i:int, max_samples:int): seq[Sample] =
  # given sample_i, get related samples to show in the plot
  let sample = samples[sample_i]
  result.add(sample)

  if sample.dad != nil and sample.dad.i >= 0:
    result.add(sample.dad)
  if sample.mom != nil and sample.mom.i >= 0:
    result.add(sample.mom)

  # add kids:
  for k in sample.kids:
    if result.len < max_samples:
      result.add(k)

  # add siblings
  if result.len > 1:
    for k in sample.siblings:
      if k.i >= 0 and k.i != sample.i and result.len < max_samples:
        result.add(k)

  # TODO: add unrelated samples.


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

proc encode*(variant:Variant, ivcf:VCF, bams:TableRef[string, Bam], fasta:Fai, samples:seq[pedfile.Sample], sample_i:int,
             anno_files:seq[string], note:string="", max_samples:int=5, flank:int=120, single_locus:string=""): JsonNode =
  # single_locus is used when we don't want to specify a vcf
  # TODO: if region is too large, try multi-locus:
  # https://igv.org/web/release/2.8.4/examples/multi-locus.html
  # small locus is for the initial view.
  var small_locus, locus: string
  if variant == nil:
    doAssert single_locus != "", "[jigv] expected single_locus to be specified since no variant was given"
    small_locus = single_locus
    locus = single_locus
  else:
    small_locus = &"{variant.CHROM}:{max(1, variant.start - 20)}-{variant.stop + 20}"
    # locus how much data we pull (and how far user can zoom out).
    locus = &"{variant.CHROM}:{max(1, variant.start - flank)}-{variant.stop + flank}"
  var json:JsonNode = %* {
      "locus": small_locus,
      "reference": {"fastaURL": fasta.encode(locus) },
      "queryParametersSupported": true,
      "showChromosomeWidget": false,
    }

  var tracks: seq[Track]
  # TODO: set n_tracks
  # TODO: set name to genotype + GQ + AD
  if ivcf != nil:
    var tr = Track(name:extractFileName(ivcf.fname), path:ivcf.encode(locus), n_tracks:2, file_type:FileType.VCF, region:locus)
    tracks.add(tr)

  for sample in samples:
    var ibam = bams[sample.id]
    var tr = Track(name:sample.id, path: ibam.encode(locus), n_tracks:2, file_type:FileType.BAM, region:locus)
    tracks.add(tr)

  json["tracks"] = %* tracks
  return json

proc samplename(ibam:Bam): string =
  var found = initHashSet[string]()
  for l in ($ibam.hdr).split('\n'):
    if not l.startswith("@RG"): continue
    for t in l.split('\t'):
      if t.startswith("SM:"):
        var p = t.split(':')
        if p[1] notin found:
          result = p[1]
          found.incl(p[1])
  doAssert found.len == 1, &"[jigv] found {found} sample names (SM read-group tags in bam header), expected exactly one."


proc get_html(): string =
  const templ = staticRead("jigv-template.html")
  if getEnv("JIGV_TEMPLATE") != "":
   if not existsFile(getEnv("JIGV_TEMPLATE")):
     stderr.write_line "[jigv] template not found in value given in JIGV_TEMPLATE: {getEnv(\"JIGV_TEMPLATE\")}"
   else:
     result = readFile(getEnv("JIGV_TEMPLATE"))
     return
  return templ


proc main*(args:seq[string]=commandLineParams()) =

  var p = newParser("samplename"):
    option("--sample", help="sample-id for proband or sample of interest (default is first vcf sample)")
    option("--sites", help="VCF containing variants of interest for --sample. if this contains ':', then it's used as a single region and the first bam/cram given is the sample of interest.")
    # TODO: option("--js", help="custom javascript to load")
    option("-g", "--genome-build", help="genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json). with this gene models will be drawn")
    option("--ped", help="pedigree file used to find relations for --sample")
    option("--fasta", help="path to indexed fasta file; required for cram files")
    arg("xams", nargs= -1, help="indexed bam or cram files for relevant samples. read-groups must match samples in vcf.")


  const max_samples = 5
  try:
    var opts = p.parse(args)
    if opts.help: quit 0
    if opts.xams.len == 0:
      stderr.write_line p.help
      quit "specify at least 1 bam or cram file"

    if opts.sites == "":
      stderr.write_line p.help
      quit "specify a vcf file to --sites"

    var bams: TableRef[string, Bam] = newTable[string, Bam]()
    var firstbam:Bam
    for xam in opts.xams:
      var ibam:Bam
      if not ibam.open(xam, fai=opts.fasta, index=true):
        quit &"[jigv] couldnt open {xam}"
      if bams.len == 0: firstbam = ibam
      bams[ibam.samplename] = ibam

    var ivcf:VCF
    var samples:seq[Sample]
    var sample_i:int

    if ':' notin opts.sites:
      if not ivcf.open(opts.sites):
        quit &"[jigv] could not open {opts.sites}"

      if opts.ped == "":
        if opts.sample != "":
          samples = @[Sample(id: opts.sample, i:0)]
        else:
          samples = @[Sample(id: ivcf.samples[0], i:0)]
        sample_i = 0

      else:
        samples = parse_ped(opts.ped).match(ivcf)
        sample_i = if opts.sample != "": ivcf.samples.find(opts.sample) else: 0
        samples = samples.get_samples(sample_i, max_samples).match(ivcf)

      var sample_ids: seq[string]
      for s in samples: sample_ids.add(s.id)
      ivcf.set_samples(sample_ids)

    else:
      sample_i = 0
      samples = @[Sample(id: firstbam.samplename, i:0)]

    var fa:Fai
    if opts.fasta != "":
      if not fa.open(opts.fasta):
        quit &"[jigv] could not open {opts.fasta}"

    var ifiles: seq[string]
    var sessions: seq[string]

    if ivcf != nil:
      for v in ivcf:
        var tracks = v.encode(ivcf, bams, fa, samples, sample_i, ifiles)
        var s = ($tracks).encode
        sessions.add(s)
        if sessions.len > 1000: break
    else:
      var v:Variant
      var tracks = v.encode(ivcf, bams, fa, samples, sample_i, ifiles, single_locus=opts.sites)
      var s = ($tracks).encode
      sessions.add(s)

    stderr.write_line opts.genome_build
    stderr.write_line &"[jigv] writing {sessions.len} regions to html"
    var meta_options = %* {
      "showChromosomeWidget": false,
      "search": false,
      #"sessionURL": % encode($(options)),
      "showCursorTrackingGuide": true,
      "showChromosomeWidget": false,
      "queryParametersSupported": true,
    }
    meta_options["sessions"] = %* sessions
    if opts.genome_build != "":
      meta_options["genome"] = % opts.genome_build

    var index_html = get_html().replace("<OPTIONS>", $(meta_options)).replace("<JIGV_CUSTOM_JS>", "")
    echo index_html

  except UsageError as e:
    stderr.write_line(p.help)
    stderr.write_line(getCurrentExceptionMsg())
    quit(1)


when isMainModule:
  main()
