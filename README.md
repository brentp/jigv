# jigv

```
jigv *.bam *.cram *.vcf.gz
```
then go to localhost:5000 and see your files in [igv.js](https://github.com/igvteam/igv.js)

# options

```
Usage:
  jigv [options] [files ...]

Arguments:
  [files ...]      bam/cram/vcf file(s) (with indexes)

Options:
  -r, --region=REGION        optional region to start at (default: chr1)
  -o, --open-browser         automatically open default browser to view files
  -g, --genome-build=GENOME_BUILD
                             genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json) (default: hg38)
  -f, --fasta=FASTA          optional fasta reference file if not in hosted and need to decode CRAM
  -p, --port=PORT            (default: 5001)
  -h, --help                 Show this help
```
