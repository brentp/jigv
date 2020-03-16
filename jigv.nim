import os
import argparse
import strformat
import jester
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

proc get_type(T:Track): string =
  case T.file_type
  of FileType.CRAM, FILE_TYPE.BAM:
    return "alignment"
  of FileType.VCF, FILE_TYPE.BCF:
    return "variant"
  else:
    raise newException(ValueError, "unknown file type for " & $T.file_type)

proc height(T:Track): int =
  case T.file_type
  of FileType.CRAM, FILE_TYPE.BAM:
    return 220
  of FileType.VCF, FILE_TYPE.BCF:
    return 80
  else:
    raise newException(ValueError, "unknown format for " & $T.file_type)

proc format(T:Track): string =
  case T.file_type
  of FileType.CRAM, FILE_TYPE.BAM:
    return ($T.file_type).toLowerAscii
  of FileType.VCF, FILE_TYPE.BCF:
    return "vcf"
  else:
    raise newException(ValueError, "unknown format for " & $T.file_type)

proc url(t:Track): string =
  result = &"static/tracks/{extractFileName(t.path)}"

proc indexUrl(t:Track): string =
  result = t.url
  case t.file_type
    of FileType.CRAM:
      result &= ".crai"
    of FileType.BAM:
      result &= ".bai"
    of FileType.VCF:
      result &= ".tbi"

    else:
      raise newException(ValueError, "unsupported file type:" & $t.file_type)


proc `$`*(T:Track): string =
  result = &"""{{
  type: "{T.get_type}", format: "{T.format}",
  url: "{T.url}",
  indexURL: "{T.indexUrl}","""
  if T.height != 0:
    result &= &"\n  \"height\": {T.height},"
  if T.file_type in  {FileType.CRAM, FILE_TYPE.BAM}:
    result &= """
   showSoftClips: true,
   viewAsPairs: true,"""
  result &= "\n}"


proc `%`*(T:Track): JsonNode =

  var fields = initOrderedTable[string, JsonNode](4)
  fields["type"] = % T.get_type
  fields["format"] = % T.format
  fields["url"] = % T.url
  fields["indexUrl"] = % T.indexUrl
  if T.height != 0:
    fields["height"] = % T.height
  if T.file_type in  {FileType.CRAM, FILE_TYPE.BAM}:
    fields["showSoftClips"] = % true
    fields["viewAsPairs"] = % true
  return JsonNode(kind: JObject, fields: fields)

router igvrouter:

  get "/data/@name":
    let name = extractFileName(@"name")
    #stderr.write_line "name:", name
    var file_path: string
    for tr in tracks:
      if name == extractFileName(tr.path):
        file_path = tr.path

    let r = request.headers["Range"]
    let byte_range = r.split("=")[1].split("-")
    let offset = parseInt(byte_range[0])
    let length = parseInt(byte_range[1])

    var fh:File
    if not open(fh, file_path):
       quit "couldnt open file:" & file_path
    let size = getFileSize(fh)

    fh.setFilePos(offset)
    # not sure why we need + 1 here...
    var data = newString(length+1)
    data[data.high] = 0.char
    doAssert length == fh.readBuffer(data[0].addr.pointer, length)

    let range_str = &"bytes {offset}-{offset + length}/{size}"
    let headers = [(key:"Content-Type", value:"application/octet-stream"), (key:"Content-Range", value: range_str)]
    resp(Http206, headers, data)

  get "/":
    #echo "resuting index"
    #resp.headers["Content-Type"] = "text/html"
    resp index_html

const insert_js = """
<script type="text/javascript">
function jigv() {
    let div = document.getElementById("jigv");
    let browser = igv.createBrowser(div, <OPTIONS>)
}

</script>
"""

const templ = staticRead("jigv-template.html")

proc main() =

  let p = newParser("igv-server"):
    option("-r", "--region", help="optional region to start at")
    option("-b", "--genome-build", default="hg38")
    option("-p", "--port", default="5001")
    # TODO: regions files
    arg("files", help="bam/cram/vcf file(s) (with indexes)", nargs= -1)

  var args = p.parse()
  if args.help:
    quit(0)

  for f in args.files:
      tracks.add(Track(path:f, file_type: f.file_type, name: extractFileName(f)))

  var tmpl = templ
  if getEnv("JIGV_TEMPLATE") != "":
    if not existsFile(getEnv("JIGV_TEMPLATE")):
      stderr.write_line "[jigv] template not found in value given in JIGV_TEMPLATE: {getEnv(\"JIGV_TEMPLATE\")}"
    else:
      tmpl = readFile(getEnv("JIGV_TEMPLATE"))

  let options = %* {
      "genome": args.genome_build,
      "showCursorTrackingGuide": true,
      "tracks": tracks
      }

  index_html = tmpl.replace("</head>", insert_js.replace("<OPTIONS>", $options) & "</head>")
  index_html = index_html.replace("</body>", """</body><script type="text/javascript">jigv()</script>""")

  let settings = newSettings(port=parseInt(args.port).Port)
  var jester = initJester(igvrouter, settings)
  jester.serve()

when isMainModule:
  main()
