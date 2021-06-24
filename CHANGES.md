v0.1.7
======
+ breaking: change `--template` to `--prefix`

v0.1.5
======
+ add `--template` option to allow for each variant or site to be written to a separte js (or base64) file.
  this allows a user to scale to as many variants as they like.

v0.1.4
======
+ handle and report weird AD values

v0.1.3
======
+ expose --flank option so user can specify size of flanking region to extract
+ allow linking by chrom:pos:ref:alt (in addition to region chrom:start-stop)

v0.1.2
======
+ if no --sample use first affected sample from pedigree file


v0.1.1
======
+ support BED as input

v0.1.0
======
+ rewrite so no server is needed

v0.0.8
======
+ show help when run with no args (#7)
+ allow injection of custom javascript to affect `options` and `option.tracks` (#4)

v0.0.7
======
+ css from @brwnj
+ better positioning of variant within window when using ->, <- arrow keys for SVs
+ more concise popover for alignment track click
+ given a VCF and no `--region` specified, the viewer will open at the first variant
+ alignment track-height is larger when fewer tracks.

v0.0.6
======
+ more fixes for custom reference+cram

v0.0.5
======
+ styling and left/right buttons by Joe Brown
+ better aesthetics (sized to fit a trio of bams +variant track in a laptop screen)
+ update readme to include directions on how to generate screenshots
+ fix for custom reference + cram
+ bump igv.js version to 2.4.1

v0.0.4
======
+ navigation along VCF variants with left/right arrows.

v0.0.3
======
+ color alignments by strand and set height to 10 by default
+ guess type by file extension for remote files
