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
import hts/bam
import hts/vcf
import tables
import json
import ./enc
import ./track

var igvtracks*: seq[Track]
var index_html:string

proc read_range(file_path:string, range_req:string, extra_bytes:int=0): (seq[tuple[key: string, val: string]], string) =
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
    let headers = @[(key:"Content-Type", val:"application/octet-stream"), (key:"Content-Range", val: range_str)]
    fh.close

    return (headers, data)


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
    for tr in igvtracks:
      if not tr.path.endsWith(".vcf.gz"): continue
      loc = tr.path.findNextFeature(chrom, left, right)
      break

    let headers = [(key:"Content-Type", val:"application/json")]
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
    for otr in igvtracks:
      if name == extractFileName(otr.path):
        tr = otr
        file_path = tr.path

    if file_path == "" and "Range" in request.headers.table:
      resp(Http404)
    elif file_path == "":
      # requesting index or other full file
      for tr in igvtracks:
        if name == extractFileName(tr.indexUrl):
          file_path = tr.path & tr.index_ext
          let data = file_path.readFile
          let headers = [(key:"Content-Type", val:"application/octet-stream")]
          resp(Http200, headers, data)
          return
    elif file_path != "" and "Range" notin request.headers.table and "range" notin request.headers.table:
      stderr.write_line "slurping file:", file_path
      let data = file_path.readFile
      let headers = [(key:"Content-Type", val:"application/octet-stream")]
      resp(Http200, headers, data)
      return


    let r = request.headers["Range"]

    var (headers, data) = file_path.read_range(r, extra_bytes=int(tr.format == "cram"))
    resp(Http206, headers, data)

  get re"/reference/(.+)":
    var path = request.matches[0]
    if path.endsWith(".fai"):
      let data = path.readFile
      let headers = [(key:"Content-Type", val:"application/octet-stream")]
      resp(Http200, headers, data)
      return

    if "range" notin request.headers.table and "Range" notin request.headers.table:
      quit "can't handle reference request without range"

    let r = request.headers["Range"]
    var (headers, data) = path.read_range(r, extra_bytes=1)
    resp(Http206, headers, data)


  get "/":
    resp index_html

proc fill(templ:string, args:auto, tracks:var seq[Track], options:JsonNode, first_vcf_track_index:int, flank=100) =
  # if we are here, then we are filling the template with each track with the
  # data urls.
  doAssert args.region != ""
  for tr in tracks.mitems:
    tr.region = args.region

  options["tracks"] = %* tracks
  options["locus"] = % args.region

  var fa:Fai
  if args.fasta != "" and fa.open(args.fasta):
    options["reference"] = %* {"fastaURL": fa.encode(args.region) }
    if args.cytoband != "":
      options["reference"]["cytobandURL"] = % args.cytoband.encode(args.region, TrackFileType.cytoband)
  else:
      options["genome"] = % args.genome_build

  var sessions = newSeq[string]()

  var meta_options = %* {
    "showChromosomeWidget": false,
    "search": false,
    #"sessionURL": % encode($(options)),
    "showCursorTrackingGuide": true,
    "showChromosomeWidget": false,
    "queryParametersSupported": true,
  }

  if first_vcf_track_index >= 0:
    var ivcf:VCF
    if not ivcf.open(tracks[first_vcf_track_index].path):
      raise newException(IOError, &"[jigv] couldn't open vcf path {tracks[first_vcf_track_index].path}")
    for v in ivcf:
      if v.stop - v.start > 10000: continue
      var locus = &"{v.CHROM}:{max(0, v.start-flank)}-{v.stop + flank}"
      stderr.write_line locus
      for tr in tracks.mitems:
        tr.region = locus
      options["reference"] = %* {"fastaURL": fa.encode(locus) }
      options["tracks"] = %*tracks
      options["locus"] = % locus
      sessions.add(encode($options))
      if sessions.len > 100:
        break

  meta_options["sessionURL"] = % sessions[0]
  meta_options["sessions"] = %* sessions
  var index_html = templ.replace("<OPTIONS>", pretty(meta_options))
  var js:string = args.js
  if js.endsWith(".js"): js = readFile(js)
  index_html = index_html.replace("<JIGV_CUSTOM_JS>", js)
  echo index_html



const templ = staticRead("jigv-template.html")

proc main() =

  let p = newParser("jigv"):
    option("-r", "--region", help="optional region to start viewing (will default to first variant in first vcf)")
    flag("-o", "--open-browser", help="open browser")
    option("-g", "--genome-build", default="hg38", help="genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json)")
    option("-f", "--fasta", default="", help="optional fasta reference file if not in hosted and needed to decode CRAM")
    option("-c", "--cytoband", default="", help="optional path to cytoband bed file")
    option("-p", "--port", default="", help="if this is not specified, jigv will try to set up the page so it doesn't need a server. (otherwise, a value like '5001' is a good choice.")
    option("--js", help="custom javascript to inject. will have access to `options` and `options.tracks`. if this ends in .js it is read as a file")
    arg("files", help="bam/cram/vcf/bed file(s) (with indexes)", nargs= -1)

  var argv = commandLineParams()
  if argv.len == 0: argv.add("-h")
  var args = p.parse(argv)
  if args.help:
    quit(0)

  var first_vcf = -1
  for i, f in args.files:
    var fs = f.split("#")
    var tr = Track(name: extractFileName(f), path: f, n_tracks:args.files.len)
    if fs.len == 2:
      tr.path = fs[0]
      tr.name = fs[1]
    tr.file_type = tr.path.file_type_ez
    if first_vcf == -1 and tr.file_type == FileType.VCF:
      first_vcf = i


    igvtracks.add(tr)

  args.region = args.region.replace(",", "")
  if args.region == "" and first_vcf != -1:
    args.region = get_first_variant(igvtracks[first_vcf].path)
    stderr.write_line "set args.region to:", args.region

  var tmpl = templ
  if getEnv("JIGV_TEMPLATE") != "":
    if not existsFile(getEnv("JIGV_TEMPLATE")):
      stderr.write_line "[jigv] template not found in value given in JIGV_TEMPLATE: {getEnv(\"JIGV_TEMPLATE\")}"
    else:
      tmpl = readFile(getEnv("JIGV_TEMPLATE"))

  let options = %* {
      "tracks": igvtracks,
      "queryParametersSupported": true,
    }

  if args.port == "":
    tmpl.fill(args, igvtracks, options, first_vcf)

  else:

    if args.fasta != "":
      var rp = &"/reference/{args.fasta}"
      var rpi = &"/reference/{args.fasta}.fai"
      options["reference"] = %* {"fastaURL": rp, "indexURL": rpi, "id": extractFileName(args.fasta)}
    else:
      options["genome"] = % args.genome_build

    if args.region != "":
      options["locus"] = % args.region

    index_html = tmpl.replace("<OPTIONS>", pretty(options))
    var js:string = args.js
    if js.endsWith(".js"): js = readFile(js)
    index_html = index_html.replace("<JIGV_CUSTOM_JS>", js)

    let settings = newSettings(port=parseInt(args.port).Port)
    var jester = initJester(igvrouter, settings)
    if args.open_browser:
      openDefaultBrowser(&"http://localhost:{args.port}/")
    jester.serve()

when isMainModule:
  main()
