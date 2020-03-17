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

proc file_type_ez(path:string): FileType =
  if path.startswith("http"):
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
    return 220
  case T.file_type
  of FileType.CRAM, FileType.BAM:
    return 220
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
  if t.path.startswith("http") or t.path.startswith("ftp"):
    return t.path
  result = &"/data/tracks/{extractFileName(t.path)}"

proc index_ext*(t:Track): string =
  if t.path.startswith("http") or t.path.startswith("ftp"):
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
   alignmentRowHeight: 10,
   colorBy: "strand",
   showSoftClips: true,
   viewAsPairs: true,"""
  elif T.file_type == FileType.VCF:
    result &= """
    visibilityWindow: 200000,"""
  result &= "\n}"

proc read_range(file_path:string, range_req:string): (seq[tuple[key: string, value: string]], string) =

    let r = range_req
    let byte_range = r.split("=")[1].split("-")
    let offset = parseInt(byte_range[0])
    let length = parseInt(byte_range[1]) - offset

    var fh:File
    if not open(fh, file_path):
       quit "couldnt open file:" & file_path
    let size = getFileSize(fh)

    fh.setFilePos(offset)
    # not sure why we need + 1 here...
    var data = newString(min(size, length)+1)
    data[data.high] = 0.char
    let got = fh.readBuffer(data[0].addr.pointer, length)
    doAssert got == data.high, $(got, length, "size:" & $size)

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
    fields["alignmentRowHeight"] = % 10
    fields["colorBy"] = % "strand"
    fields["showSoftClips"] = % true
    fields["viewAsPairs"] = % true
  return JsonNode(kind: JObject, fields: fields)

router igvrouter:


  get "/data/tracks/@name":
    let name = extractFileName(@"name")
    var file_path: string
    for tr in tracks:
      if name == extractFileName(tr.path):
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

    var (headers, data) = file_path.read_range(r)
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
    var (headers, data) = path.read_range(r)
    resp(Http206, headers, data)


  get "/":
    resp index_html

const insert_js = """
<script type="text/javascript">
var browser
function jigv() {
    let div = document.getElementById("jigv");
    let options = <OPTIONS>
    if(location.hash.substr(1)) {
      options.locus = location.hash.substr(1)
    } else {
      location.hash = options.locus
    }
    browser = igv.createBrowser(div, options).then(function(browser) {
       browser.on('locuschange', function (referenceFrame) {
            location.hash = referenceFrame.label
       });
    })

}

</script>
"""

const templ = staticRead("jigv-template.html")

proc main() =

  let p = newParser("jigv"):
    option("-r", "--region", help="optional region to start at", default="chr1")
    flag("-o", "--open-browser", help="open browser")
    option("-g", "--genome-build", default="hg38", help="genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json)")
    option("-f", "--fasta", default="", help="optional fasta reference file if not in hosted and needed to decode CRAM")
    option("-p", "--port", default="5001")
    # TODO: regions files
    arg("files", help="bam/cram/vcf/bed file(s) (with indexes)", nargs= -1)

  var args = p.parse()
  if args.help:
    quit(0)

  for f in args.files:
    var fs = f.split("#")
    var tr = Track(name: extractFileName(f), path: f)
    if fs.len == 2:
      tr.path = fs[0]
      tr.name = fs[1]
    tr.file_type = tr.path.file_type_ez
    tracks.add(tr)

  var tmpl = templ
  if getEnv("JIGV_TEMPLATE") != "":
    if not existsFile(getEnv("JIGV_TEMPLATE")):
      stderr.write_line "[jigv] template not found in value given in JIGV_TEMPLATE: {getEnv(\"JIGV_TEMPLATE\")}"
    else:
      tmpl = readFile(getEnv("JIGV_TEMPLATE"))

  let options = %* {
      "genome": args.genome_build,
      "showCursorTrackingGuide": true,
      "tracks": tracks,
      "queryParametersSupported": true,
    }
  if args.fasta != "":
    var rp = &"/reference/{args.fasta}"
    var rpi = &"/reference/{args.fasta}.fai"
    options["reference"] = %* {"fastaURL": rp, "indexURL": rpi}
  if args.region != "":
    options["locus"] = % args.region

  index_html = tmpl.replace("</head>", insert_js.replace("<OPTIONS>", pretty(options)) & "</head>")
  index_html = index_html.replace("</body>", """</body><script type="text/javascript">jigv()</script>""")

  let settings = newSettings(port=parseInt(args.port).Port)
  var jester = initJester(igvrouter, settings)
  if args.open_browser:
    openDefaultBrowser(&"http://localhost:{args.port}/")
  jester.serve()

when isMainModule:
  main()
