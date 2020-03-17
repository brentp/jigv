# jigv

[igv.js](https://github.com/igvteam/igv.js) is a great way to view aligments an other files. This provides a server
and some default configuration, javascript and HTML so you can do:

```
jigv --open-browser --region chr1:34566-34999 *.bam *.cram *.vcf.gz
```
With that, a server will start and a browser will open an igv.js viewer for your requested files.

# installation

grab a static linux binary from [releases](https://github.com/brentp/jigv/releases/latest)

# options

```
Usage:
  jigv [options] [files ...]

Arguments:
  [files ...]      bam/cram/vcf/bed{,.gz} file(s) (with indexes)

Options:
  -r, --region=REGION        optional region to start at (default: chr1)
  -o, --open-browser         automatically open default browser to view files
  -g, --genome-build=GENOME_BUILD
                             genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json) (default: hg38)
  -f, --fasta=FASTA          optional fasta reference file if not in hosted and needed to decode CRAM
  -p, --port=PORT            (default: 5001)
  -h, --help                 Show this help
```

# notes

+ this is likely insecure in many ways.
+ there will soon be a way to customize the options and javascript (but this probably covers 85% of use-cases as-is).
+ if you have some custom javascript used with igv.js, that is generally useful, please open an issue so I can add it.
