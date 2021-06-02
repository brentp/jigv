import hts/files
import strutils
import tables
import json
import strformat
import hts/bam
import os

type Track* = object
  file_type*: FileType
  path*:string
  name*: string
  n_tracks*:int
  reference*:string
  # if region is specified, then we base64 encode the track, otherwise we
  # use the server.
  region*:string

proc isremote*(path:string): bool =
  result = path.startswith("http") or path.startswith("ftp:")

proc file_type_ez*(path:string): FileType =
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

proc get_type*(T:Track): string =
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

proc height*(T:Track): int =
  var mult = 1
  if T.n_tracks < 3:
    mult = 2
  if T.path.endsWith(".cram") or T.path.endsWith(".bam"):
    return 250 * mult
  case T.file_type
  of FileType.CRAM, FileType.BAM:
    return 250 * mult
  of FileType.VCF, FileType.BCF:
    return 60
  else:
    return 50

proc format*(T:Track): string =
  if T.path.split('#')[0].endsWith(".cram"):
    if T.region != "": return "bam"
    return "cram"
  if T.path.split('#')[0].endsWith(".bam"):
    return "bam"
  case T.file_type
  of FileType.CRAM, FileType.BAM:
    return ($T.file_type).toLowerAscii
  of FileType.VCF, FileType.BCF:
    return "vcf"
  else:
    if T.path.split('#')[0].endsWith(".bed") or T.path.split('#')[0].endsWith(".bed.gz") or T.path.split('#')[0].endsWith(".bedgraph"):
      return "bed"
    if T.path.split('#')[0].endsWith(".bw") or T.path.split('#')[0].endsWith(".wig") or T.path.split('#')[0].toLowerAscii.endswith(".bigwig"):
      return "wig"
    raise newException(ValueError, "unknown format for " & $T.file_type)

template isdata(path:string): bool =
  path.startswith("data:")

proc url(t:Track): string =
  if t.path.isremote or t.path.isdata:
    return t.path
  result = &"/data/tracks/{extractFileName(t.path)}"

#[
proc dataurl(t:Track): string =

  # region is the swithc. if we have it, we use it and encode
  # otherwise, we rely on the server
  if t.region == "": return t.url

  case t.format
  of "bam", "cram":
    var ibam:Bam
    if not ibam.open(t.path, fai=t.reference, index=true):
      quit &"couldnt open bam/cram file {t.path} with reference {t.reference}"
    result = ibam.encode(t.region)
    ibam.close()
  of "vcf", "bcf":
    var ivcf:VCF
    if not ivcf.open(t.path):
      quit &"couldnt open vcf/bcf file: {t.path}"
    result = ivcf.encode(t.region)
    ivcf.close()
  of "bed":
    # TODO: check for index and use bed12
    result = t.path.encode(t.region, TrackFileType.bed)

  else:
    return t.url
    ]#

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

proc indexUrl*(t:Track): string =
  result = t.url & t.index_ext

proc `%`*(T:Track): JsonNode =
  var fields = initOrderedTable[string, JsonNode](4)
  fields["type"] = % T.get_type
  if T.format != "wig":
    fields["format"] = % T.format
  fields["url"] = % T.url
  if T.region == "":
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

proc `$`*(T:Track): string =
  return $(%T)
