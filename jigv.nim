import os
import argparse
import strformat
import jester
import browsers

#import regex
import re
import hts
import httpcore
import hts/files
import tables
import json

type Track* = object
  file_type: FileType
  path*:string
  name*: string

var tracks*: seq[Track]
var index_html: string

proc isremote(path:string): bool =
  result = path.startswith("http") or path.startswith("ftp:")

proc file_type_ez(path:string): FileType =
  if path.isremote:
    if path.endswith(".bam"):
      return FileType.BAM
    elif path.endswith(".cram"):
      return FileType.CRAM
    elif path.endswith(".vcf") or path.endswith(".vcf.gz"):
      return FileType.VCF

    return FileType.UNKNOWN

  else:
    return path.file_type

proc get_type(T:Track): string =
  if T.path.endsWith(".cram") or T.path.endsWith(".bam"):
    return "alignment"
  case T.file_type
  of FileType.CRAM, FileType.BAM:
    return "alignment"
  of FileType.VCF, FileType.BCF:
    return "variant"
  else:
    if T.path.endsWith(".bed") or T.path.endsWith(".bed.gz") or T.path.endsWith(".bedgraph"):
      return "annotation"
    if T.path.endsWith(".bw") or T.path.endsWith(".wig") or T.path.toLowerAscii.endswith(".bigwig"):
      return "wig"
    raise newException(ValueError, "unknown file type for " & $T.file_type)

proc height(T:Track): int =
  if T.path.endsWith(".cram") or T.path.endsWith(".bam"):
    return 250
  case T.file_type
  of FileType.CRAM, FileType.BAM:
    return 250
  of FileType.VCF, FileType.BCF:
    return 60
  else:
    return 50

proc format(T:Track): string =
  if T.path.endsWith(".cram"):
    return "cram"
  if T.path.endsWith(".bam"):
    return "bam"
  case T.file_type
  of FileType.CRAM, FileType.BAM:
    return ($T.file_type).toLowerAscii
  of FileType.VCF, FileType.BCF:
    return "vcf"
  else:
    if T.path.endsWith(".bed") or T.path.endsWith(".bed.gz") or T.path.endsWith(".bedgraph"):
      return "bed"
    if T.path.endsWith(".bw") or T.path.endsWith(".wig") or T.path.toLowerAscii.endswith(".bigwig"):
      return "wig"
    raise newException(ValueError, "unknown format for " & $T.file_type)

proc url(t:Track): string =
  if t.path.isremote:
    return t.path
  result = &"/data/tracks/{extractFileName(t.path)}"

proc index_ext*(t:Track): string =
  if t.path.isremote:
    case splitFile(t.path).ext
    of ".bam": return ".bai"
    of ".cram": return ".crai"
    of ".vcf.gz", "bed.gz": return ".tbi"

  case t.file_type
  of FileType.CRAM:
    result = ".crai"
  of FileType.BAM:
    result = ".bai"
  of FileType.VCF:
    result = ".tbi"
  else:
    if t.path.endsWith(".bed.gz"):
      return ".tbi"
    raise newException(ValueError, "unsupported file type:" & $t.file_type)

proc indexUrl(t:Track): string =
  result = t.url & t.index_ext

proc `$`*(T:Track): string =
  result = &"""{{
  type: "{T.get_type}", format: "{T.format}",
  url: "{T.url}","""
  try:
    result &= """
  indexURL: "{T.indexUrl}","""
  except ValueError:
    discard

  if T.height != 0:
    result &= &"\n  \"height\": {T.height},"
  if T.file_type in  {FileType.CRAM, FileType.BAM}:
    result &= """
   alignmentRowHeight: 8,
   colorBy: "strand",
   coverageTrackHeight: 24,
   showSoftClips: true,
   viewAsPairs: true,"""
  elif T.file_type == FileType.VCF:
    result &= """
    visibilityWindow: 200000,"""
  result &= "\n}"

proc read_range(file_path:string, range_req:string, extra_bytes:int=0): (seq[tuple[key: string, value: string]], string) =
    ## the fasta request doesn't match the bam request so we use `extra_bytes` to
    ## return 1 extra byte

    let r = range_req
    let byte_range = r.split("=")[1].split("-")
    let offset = parseInt(byte_range[0])
    let length = parseInt(byte_range[1]) - offset + extra_bytes

    var fh:File
    if not open(fh, file_path):
       quit "couldnt open file:" & file_path
    let size = getFileSize(fh)

    fh.setFilePos(offset)
    # not sure why we need + 1 here...
    var data = newString(min(size, length+1))
    data[data.high] = 0.char
    let got = fh.readBuffer(data[0].addr.pointer, length)
    if got != data.high:
      #, $(got, length, "size:" & $size, data.high)
      data.setLen(got+1)

    let range_str = &"bytes {offset}-{offset + length}/{size}"
    let headers = @[(key:"Content-Type", value:"application/octet-stream"), (key:"Content-Range", value: range_str)]
    fh.close

    return (headers, data)

proc `%`*(T:Track): JsonNode =

  var fields = initOrderedTable[string, JsonNode](4)
  fields["type"] = % T.get_type
  if T.format != "wig":
    fields["format"] = % T.format
  fields["url"] = % T.url
  try:
    fields["indexURL"] = % T.indexUrl
  except ValueError:
    discard

  if fields["type"] == % "annotation":
    fields["displayMode"] = % "SQUISHED"
  fields["name"] = % T.name
  if T.height != 0:
    fields["height"] = % T.height
  if T.file_type in  {FileType.CRAM, FileType.BAM}:
    fields["alignmentRowHeight"] = % 8
    fields["colorBy"] = % "strand"
    fields["coverageTrackHeight"] = % 24
    fields["showSoftClips"] = % true
    fields["viewAsPairs"] = % true
  return JsonNode(kind: JObject, fields: fields)

type region = tuple[chrom:string, start:int, stop:int]


template stripChr*[T:string|cstring](s:T): string =
  if s.len > 3 and ($s).startswith("chr"): ($s)[3..<s.len] else: $s

proc get_first_variant(path:string): string =
  ## get the first variant in a VCF file
  var ivcf:VCF
  if not ivcf.open(path):
    return ""
  defer:
    ivcf.close

  for v in ivcf:
    var p = v.start
    return &"{v.CHROM}:{p - 79}-{p + 80}"

proc findNextFeature(ivcf:VCF, chrom:string, left:int, right:int, rep:int): region =
  var contig_i = 0
  for i, c in ivcf.contigs:

    if stripChr(c.name) == stripChr(chrom):
      result.chrom = c.name
      contig_i = i
      break

  if result.chrom == "":
    stderr.write_line &"[jigv] chrom not found {chrom}"
    return

  var center = int(0.5 + (left + right) / 2)

  for v in ivcf.query(&"{result.chrom}:{left}"):
    if v.start < center + 1: continue

    var desired_range = right - left
    if v.stop - v.start < 20:
      var vmid = int(0.5 + int(v.start + v.stop) / 2)

      var delta = int(0.5 + float(desired_range - int(v.stop - v.start)) * 0.5)

      return (chrom, vmid - delta, vmid + delta)


    else:
      var offset = 50
      if v.start.int - offset <= left: continue
      # then show left end centered in window
      return (chrom, int(v.start.float - desired_range / 2), int(v.start.float + desired_range / 2))


  var contigs = ivcf.contigs
  if contig_i == contigs.high: contig_i = 0 else: contig_i.inc

  # after incrementing to next chrom, we can continue trying.
  # `rep` avoids infinite loop
  if rep < contigs.len:
    return ivcf.findNextFeature(contigs[contig_i].name, 0, 0, rep+1)

proc findNextFeature(path:string, chrom:string, left:int, right:int): region =
  var ivcf:VCF
  if not ivcf.open(path):
    stderr.write_line &"[jigv] requested path {path} not found"

  defer: ivcf.close

  return ivcf.findNextFeature(chrom, left, right, 0)

router igvrouter:

  get "/nextfeature/":
    # user clicked right in browser. this finds next feature
    # in first track that's BED or VCF
    var locus = request.params["locus"].replace(",", "").split(":")
    var chrom = locus[0]
    var left = 0
    var right = 100
    if locus.len > 1:
      var se = locus[1].split("-")
      left = parseInt(se[0])
      if se.len > 1:
        right = parseInt(se[1])
    var loc: region
    for tr in tracks:
      if not tr.path.endsWith(".vcf.gz"): continue
      loc = tr.path.findNextFeature(chrom, left, right)
      break

    let headers = [(key:"Content-Type", value:"application/json")]
    if loc.chrom == "":
      var data = %* {
           "success": false
      }
      resp(Http200, headers, $data)
      return

    let data = %* {
      "success": true,
      "position": &"{loc.chrom}:{loc.start}-{loc.stop}",
    }
    resp(Http200, headers, $data)


  get "/data/tracks/@name":
    let name = extractFileName(@"name")
    var file_path: string
    var tr: Track
    for otr in tracks:
      if name == extractFileName(otr.path):
        tr = otr
        file_path = tr.path

    if file_path == "" and "Range" in request.headers.table:
      resp(Http404)
    elif file_path == "":
      # requesting index or other full file
      for tr in tracks:
        if name == extractFileName(tr.indexUrl):
          file_path = tr.path & tr.index_ext
          let data = file_path.readFile
          let headers = [(key:"Content-Type", value:"application/octet-stream")]
          resp(Http200, headers, data)
          return
    elif file_path != "" and "Range" notin request.headers.table and "range" notin request.headers.table:
      stderr.write_line "slurping file:", file_path
      let data = file_path.readFile
      let headers = [(key:"Content-Type", value:"application/octet-stream")]
      resp(Http200, headers, data)
      return


    let r = request.headers["Range"]

    var (headers, data) = file_path.read_range(r, extra_bytes=int(tr.format == "cram"))
    resp(Http206, headers, data)

  get re"/reference/(.+)":
    var path = request.matches[0]
    if path.endsWith(".fai"):
      let data = path.readFile
      let headers = [(key:"Content-Type", value:"application/octet-stream")]
      resp(Http200, headers, data)
      return

    if "range" notin request.headers.table and "Range" notin request.headers.table:
      quit "can't handle reference request without range"

    let r = request.headers["Range"]
    var (headers, data) = path.read_range(r, extra_bytes=1)
    resp(Http206, headers, data)


  get "/":
    resp index_html

const templ = staticRead("jigv-template.html")

proc main() =

  let p = newParser("jigv"):
    option("-r", "--region", help="optional region to start at")
    flag("-o", "--open-browser", help="open browser")
    option("-g", "--genome-build", default="hg38", help="genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json)")
    option("-f", "--fasta", default="", help="optional fasta reference file if not in hosted and needed to decode CRAM")
    option("-p", "--port", default="5001")
    # TODO: regions files
    arg("files", help="bam/cram/vcf/bed file(s) (with indexes)", nargs= -1)

  var args = p.parse()
  if args.help:
    quit(0)

  var first_vcf = -1
  for i, f in args.files:
    var fs = f.split("#")
    var tr = Track(name: extractFileName(f), path: f)
    if fs.len == 2:
      tr.path = fs[0]
      tr.name = fs[1]
    tr.file_type = tr.path.file_type_ez
    if first_vcf == -1 and tr.file_type == FileType.VCF:
      first_vcf = i


    tracks.add(tr)

  if args.region == "" and first_vcf != -1:
    args.region = get_first_variant(tracks[first_vcf].path)
    echo "set args.region to:", args.region

  var tmpl = templ
  if getEnv("JIGV_TEMPLATE") != "":
    if not existsFile(getEnv("JIGV_TEMPLATE")):
      stderr.write_line "[jigv] template not found in value given in JIGV_TEMPLATE: {getEnv(\"JIGV_TEMPLATE\")}"
    else:
      tmpl = readFile(getEnv("JIGV_TEMPLATE"))

  let options = %* {
      "showCursorTrackingGuide": true,
      "tracks": tracks,
      "queryParametersSupported": true,
    }
  if args.fasta != "":
    var rp = &"/reference/{args.fasta}"
    var rpi = &"/reference/{args.fasta}.fai"
    options["reference"] = %* {"fastaURL": rp, "indexURL": rpi, "id": extractFileName(args.fasta)}
  else:
    options["genome"] = % args.genome_build

  if args.region != "":
    options["locus"] = % args.region

  index_html = tmpl.replace("<OPTIONS>", pretty(options))
#   index_html = tmpl.replace("</head>", insert_js.replace("<OPTIONS>", pretty(options)) & "</head>")
#   index_html = index_html.replace("</body>", """</body><script type="text/javascript">jigv()</script>""")

  let settings = newSettings(port=parseInt(args.port).Port)
  var jester = initJester(igvrouter, settings)
  if args.open_browser:
    openDefaultBrowser(&"http://localhost:{args.port}/")
  jester.serve()

when isMainModule:
  main()
