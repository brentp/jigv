import hts/bam
import hts/vcf except Header
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

proc encode*(ibam:Bam, region:string, hitid:int=int.high, drop_tags:seq[string]= @["RG", "PG", "MD", "XS" ,"AS"]): string =

  var obam:Bam
  var path = os.getTempDir() & "/" &  $(rand(int.high)) & ".bam"

  defer:
    discard os.tryRemoveFile(path)

  if not obam.open(path, mode="wb"):
    quit "could not open open bam"

  # when writing small bams, most of the space is actually the bam header.
  # so we drop all PG lines and take only name and len from the SQ lines.
  # on a test-set, this drops the size of a generated html from 127MB to 31MB.
  # using hitid to drop unused SQs (e.g. HLA, random chroms) that are not seen
  # as mates for the given regions drops this further to ~5MB.
  var new_header: seq[string]
  var new_line: seq[string]
  var ntids:int = 0
  for line in ($(ibam.hdr)).split("\n"):
    if line.startswith("@PG"): continue
    if line.startswith("@SQ"):
      ntids += 1
      new_line.setLen(0)
      for t in line.split("\t"):
        if t[0] == '@' or t.startswith("SN") or t.startswith("LN"):
          new_line.add(t)
      new_header.add(new_line.join("\t"))
      if ntids > hitid: break
    else:
      new_header.add(line)

  var h:bam.Header = Header()
  h.from_string(new_header.join("\n"))

  obam.write_header(h)

  # handle chr prefix stuff
  let chromse = region.split(':')
  let chrom = check_chrom(ibam.hdr.targets, chromse[0])
  var region = &"{chrom}:{chromse[1]}"
  var highest_tid = 0

  if hitid == int.high:
    for aln in ibam.query(region):
      highest_tid = max(aln.mate_tid, max(aln.tid, highest_tid))
      obam.write(aln)
  else:
    highest_tid = hitid
    for aln in ibam:
      for tag in drop_tags:
        discard aln.delete_tag(tag)
      obam.write(aln)
  obam.close()

  # we can make this even smaller by dropping SQ values from
  # the header that we don't need.
  # so we recurse and drop extra sq lines.
  if highest_tid < h.targets.len - 20 and hitid == int.high:
    var ibam2:Bam
    if not ibam2.open(path):
        quit "could not open input tmp bam"
    return ibam2.encode(region, highest_tid, drop_tags)

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

proc encode*(fai:Fai, region:string, expand:int=200): string =

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
  indexed

proc get_samples(samples:seq[Sample], sample_i:int, max_samples:int): seq[Sample] =
  # given sample_i, get related samples to show in the plot
  let sample = samples[sample_i]
  sample.extra.add(Pair(key: "label", val: sample.id))
  result.add(sample)

  for (lbl, parent) in [("dad", sample.dad), ("mom", sample.mom)]:
    if parent == nil: continue
    var aff = if parent.affected: "(affected)" else: ""
    parent.extra.add(Pair(key: "label", val: &"{lbl}:{parent.id}{aff}"))
    result.add(parent)

  # add kids:
  for k in sample.kids:
    if result.len < max_samples:
      var aff = if k.affected: "(affected)" else: ""
      k.extra.add(Pair(key: "label", val: &"offspring:{k.id}{aff}"))
      result.add(k)

  # add siblings
  if result.len > 1:
    for k in sample.siblings:
      if k.i >= 0 and k.i != sample.i and result.len < max_samples:
        var aff = if k.affected: "(affected)" else: ""
        k.extra.add(Pair(key: "label", val: &"sibling:{k.id}{aff}"))
        result.add(k)
  # TODO: add unrelated samples.

proc encode*(path:string, region:string, typ:TrackFileType): string =

  let chrom = stripChr(region.split(":")[0])
  let start = parseInt(region.split(":")[1].split("-")[0])
  let stop =  parseInt(region.split(":")[1].split("-")[1])

  case typ
  of TrackFileType.cytoband, TrackFileType.bed, TrackFileType.bed12:
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

  of TrackFileType.indexed:

    var bgz:BGZI
    if not bgz.open(path):
      stderr.write_line &"[jigv] warning: {path} should be compressed and with csi index. trying (slow) full search over lines instead."
      return path.encode(region, TrackFileType.bed)

    var chromstuff = region.split(":")
    var start = parseInt(chromstuff[1].split('-')[0])
    var stop = parseInt(chromstuff[1].split('-')[1])
    var clines:seq[string]
    # linear search...
    if chromstuff[0] notin bgz.csi.chroms:
      if not chromstuff[0].startswith("chr"):
        chromstuff[0] = "chr" & chromstuff[0]
      elif len(chromstuff[0]) > 3:
        chromstuff[0] = chromstuff[0][3..<chromstuff[0].len]
    for l in bgz.query(chromstuff[0], start - 1, stop):
      clines.add(l)
    var tmp = clines.join("\n") & "\n"
    return prefix & base64.encode(compress(tmp))

proc get_display_name(v:Variant): string =
  var r = v.REF
  const max_len = 12
  if len(r) > max_len:
    r = r[0 ..< int(max_len/2)] & &".." & r[^int(max_len/2) ..< ^0]
  var alts:seq[string]
  for a in v.ALT:
    var a = a
    if len(a) > max_len:
      a = a[0 ..< int(max_len/2)] & &".." & a[^int(max_len/2) ..< ^0]
    alts.add(a)

  result = &"{v.CHROM}:{v.start + 1}({r}/{alts.join(\",\")})"

proc is_del(v:Variant): bool {.inline.} =
  if v.ALT[0][0] == '<':
    return v.ALT[0].startswith("<DEL")
  return v.REF.len > v.ALT[0].len

proc getAB(v:Variant): seq[float32] =
  if v.format.get("AB", result) != Status.OK:

    var ad: seq[int32]
    if v.format.get("AD", ad) != Status.OK:
      return
    result = newSeq[float32](v.n_samples)
    for i in 0..<v.n_samples:
      result[i] = ad[2*i+1].float32 / max(1, ad[2*i+1] + ad[2*i]).float32
  for ab in result.mitems:
    if ab < 0: ab = 0
    if ab > 1: ab = 1

proc ff(f:float32): string {.inline.} =
  if f > 0.01 or f == 0:
    result = &"{f:.2f}"
  else:
    result = &"{f:.3f}"
  if result.endswith(".00"):
    result = result[0 ..< ^3]

proc encode*(variant:Variant, ivcf:VCF, bams:TableRef[string, Bam], fasta:Fai, samples:seq[pedfile.Sample],
             cytoband:string="",
             anno_files:seq[string], note:string="", max_samples:int=5, flank:int=100, single_locus:string=""): JsonNode =
  # single_locus is used when we don't want to specify a vcf
  # TODO: if region is too large, try multi-locus:
  # TODO: handle slivar fields e.g. show that the variant is de novo or
  # compound-het
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
  stderr.write_line "locus:", locus
  var json:JsonNode = %* {
      "locus": small_locus,
      #"search": true,
      #"queryParametersSupported": true,
      #"showChromosomeWidget": false,
      "sampleNameViewportWidth": 512,
    }
  if fasta != nil:
    json["reference"] = %* {"fastaURL": fasta.encode(locus) }
    if cytoband != "":
      json["reference"]["cytobandURL"] = % cytoband.encode(locus, TrackFileType.cytoband)

  var tracks: seq[Track]
  var n_tracks = 0
  for s in samples:
    n_tracks += int(s.id in bams)

  var x: seq[int32]
  var GTs: Genotypes
  var GQs: seq[int32]
  var ABs: seq[float32]

  if ivcf != nil:
    let variant_name = variant.get_display_name
    var fname = ivcf.fname.splitFile.name
    if fname.endsWith(".vcf") or fname.endsWith(".bcf"):
      fname = fname[0 ..< ^4]

    var tr = Track(name:fname & "<br>" & variant_name , path:ivcf.encode(locus), n_tracks:n_tracks, file_type:FileType.VCF, region:locus)
    tracks.add(tr)
    GTs = variant.format.genotypes(x)
    discard variant.format.get("GQ", GQs)
    ABs = variant.getAB()

  for sample in samples:
    if sample.id notin bams: continue
    var ibam = bams[sample.id]
    var name: string
    try:
        name = sample["label"]
    except:
        name = sample.id & "<br>"
    if GTs.len > 0:
      name &= &" GT:<b>{GTs[sample.i]}</b>"
    if GQs.len > 0:
      name &= &" GQ:<b>{GQs[sample.i]}</b>"
    if ABs.len > 0:
      name &= &" AB:<b>{ff(ABs[sample.i])}</b>"

    var tr = Track(name:name, path: ibam.encode(locus), n_tracks:n_tracks, file_type:FileType.BAM, region:locus)
    tracks.add(tr)

  for a in anno_files:
    let ft = if a.endswith(".vcf") or a.endswith(".vcf.gz") or a.endswith(".bcf"): FileType.VCF else: FileType.BED
    var tr = Track(name: extractFileName(a), path:a.encode(locus, TrackFileType.indexed), file_type: ft, region:locus)
    tracks.add(tr)

  json["tracks"] = %* tracks
  for tr in json["tracks"]:
    if $(tr["type"]) == "\"alignment\"" and variant != nil:
      tr["sort"] = %* {
        "chr": $variant.CHROM,
        "position": variant.start + (if variant.is_del: 2 else: 1),
        "option": "BASE",
        "direction": "ASC", # with ASC, igv.js always puts the alt base first.
      }
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
   if not fileExists(getEnv("JIGV_TEMPLATE")):
     stderr.write_line "[jigv] template not found in value given in JIGV_TEMPLATE: {getEnv(\"JIGV_TEMPLATE\")}"
   else:
     result = readFile(getEnv("JIGV_TEMPLATE"))
     return
  return templ

iterator generate_sites(path_or_region:string): string =
  if ':' in path_or_region:
    yield path_or_region
  else:
    for l in path_or_region.hts_lines:
      if l[0] == '#': continue
      var toks = l.split("\t")
      if toks.len < 3:
        stderr.write_line &"[jigv] WARNING line: {l} not recogized as bed line. skipping"
      yield &"{toks[0]}:{parseInt(toks[1])+1}-{parseInt(toks[2])}"

proc main*(args:seq[string]=commandLineParams()) =

  var p = newParser("jigv"):
    option("--sample", help="sample-id for proband or sample of interest (default is first vcf sample)")
    option("--sites", help="VCF or BED containing variants of interest for --sample. If this contains ':', then it's used as a single region and the first bam/cram given is the sample of interest. If it ends with '.bed' or '.bed.gz' it's assumed to be BED format.")
    # TODO: option("--js", help="custom javascript to load")
    option("-g", "--genome-build", help="genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json).  If this is specified then the page will request fasta, ideogram and gene data from a server.")
    option("--cytoband", help="optional path to cytoband/ideogram file")
    option("--annotation", help="path to additional bed or vcf file to be added as a track; may be specified multiple times", multiple=true)
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

    if opts.genome_build != "" and opts.cytoband != "":
      stderr.write_line "[jigv] warning: when -g/--genome-build is specified the cytoband argument is not used"
      opts.cytoband = ""

    var bams: TableRef[string, Bam] = newTable[string, Bam]()
    var firstbam:Bam
    for xam in opts.xams:
      var ibam:Bam
      if not ibam.open(xam, fai=opts.fasta, index=true):
        quit &"[jigv] couldnt open {xam}"
      if bams.len == 0: firstbam = ibam
      bams[ibam.samplename] = ibam

    var
      ivcf:VCF
      ivcf2:VCF # we iterate over the other vcf at each step so need to handles
    var samples:seq[Sample]
    var sample_i:int

    if ':' notin opts.sites and not (opts.sites.endsWith(".bed") or opts.sites.endsWith(".bed.gz")):
      if not ivcf.open(opts.sites):
        quit &"[jigv] could not open {opts.sites}"
      if not ivcf2.open(opts.sites):
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
    if opts.fasta != "" and opts.genome_build == "":
      if not fa.open(opts.fasta):
        quit &"[jigv] could not open {opts.fasta}"

    var ifiles: seq[string] = opts.annotation
    var sessions: seq[string]

    if sample_i != 0:
      swap(samples[0], samples[sample_i])
      sample_i = 0

    # this is used in the browser so a user can link
    # to a specific location
    var loc2idx = newTable[string, int]()

    if ivcf != nil:
      for v in ivcf:
        var tracks = v.encode(ivcf2, bams, fa, samples, anno_files=ifiles, cytoband=opts.cytoband)
        if opts.genome_build != "":
          tracks["genome"] = % opts.genome_build
        var s = ($tracks).encode
        loc2idx[($(tracks["locus"].str)).replace(",", "")] = loc2idx.len
        sessions.add(s)
        if sessions.len >= 1000: break
    else:
      var v:Variant

      for locus in generate_sites(opts.sites):
        var tracks = v.encode(ivcf, bams, fa, samples, anno_files=ifiles, cytoband=opts.cytoband, single_locus=locus)
        if opts.genome_build != "":
          tracks["genome"] = % opts.genome_build
        var s = ($tracks).encode
        sessions.add(s)

    stderr.write_line opts.genome_build
    stderr.write_line &"[jigv] writing {sessions.len} regions to html"
    if ivcf != nil:
      ivcf.close()
      ivcf2.close()
    var meta_options = %* {
      "showChromosomeWidget": false,
      "search": true,
      #"sessionURL": % encode($(options)),
      "showCursorTrackingGuide": true,
      "showChromosomeWidget": false,
      "queryParametersSupported": true,
    }
    meta_options["sessions"] = %* sessions
    meta_options["loc2idx"] = %* loc2idx

    var index_html = get_html().replace("<OPTIONS>", pretty(meta_options)).replace("<JIGV_CUSTOM_JS>", "")
    echo index_html

  except UsageError as e:
    stderr.write_line(p.help)
    stderr.write_line(getCurrentExceptionMsg())
    quit(1)


when isMainModule:
  main()
