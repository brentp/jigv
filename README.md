[igv.js](https://github.com/igvteam/igv.js) is a great way to view aligments and other genomic data. 
It requires that the files are hosted on a server with access to the original data.
javascript.

`jigv` encodes all variants, alignments, and annotations into a single HTML page that you can
to collaborators who don't have access to the cluster where your data is stored.

The resulting file is very fast to navigate; the left/right arrow keys advance to next/previous variants of interest.

Usage looks like:

```
jigv \
    --sample HG01053 \     # the sample of interest drawn in top panel
    --ped trio.ped \       # jigv will use this to also show parents and sibs of --sample
    --sites dn.vcf.gz \    # e.g. candidate de novo varants
    --fasta $reference_fasta \
    --annotation hg38.refGene.bed.gz \         # see: https://github.com/brentp/jigv/wiki/bed12
    --annotation /data/human/LCR-hs38.bed.gz \ # specify as many of these as needed.
   > denovos.html
```
With that, `denovos.html` will contain **all genomic data around variants of interest** embedded within it.

With this, we are able to encode 924 candidate *de novo* variants into a 31MB html file that includes
alignments for the proband, mom, and dad.

# installation

grab a static linux binary from [releases](https://github.com/brentp/jigv/releases/latest)

# notes and features

#### navigation

`jigv` supports navigation along variants in a VCF file using the left and right arrow keys. it uses the --sites VCF track
as the one to navigate through.

#### automated screenshots

given an html file at e.g. `denovos.html`, a PNG screenshot for `region=chr5:1410040-1412784` can be created with

```
google-chrome --window-size=1200,1200  --virtual-time-budget=10000 \
   --headless --run-all-compositor-stages-before-draw \
   --screenshot=${region/:/-}.png "denovos.html#$region"
```

This can be scripted with your favorite language for a set of regions.

# options

```
Usage:
  jigv [options] [xams ...]

Arguments:
  [xams ...]       indexed bam or cram files for relevant samples. read-groups must match samples in vcf.

Options:
  --sample=SAMPLE            sample-id for proband or sample of interest (default is first vcf sample)
  --sites=SITES              VCF containing variants of interest for --sample. if this contains ':', then it's used as a single region and the first bam/cram given is the sample of interest.
  -g, --genome-build=GENOME_BUILD
                             genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json).  If this is specified then the page will request fasta, ideogram and gene data from a server.
  --cytoband=CYTOBAND        optional path to cytoband/ideogram file
  --annotation=ANNOTATION    path to additional bed or vcf file to be added as a track; may be specified multiple times
  --ped=PED                  pedigree file used to find relations for --sample
  --fasta=FASTA              path to indexed fasta file; required for cram files
```

# limitations

+ this embeds **all** of the data in the HTML page. `jigv` tries to reduce the alignment data
  so that, for example 900 de novos variants with alignments for a trio generate an html file of
  only about 30 megabytes
+ if you have some custom javascript used with igv.js, that is generally useful, please open an issue so I can add it.
+ not all file types are supported


## See Also

+ [igv-reports](https://github.com/igvteam/igv-reports) by the IGV team does this as well. `jigv` adds extra features for
  working with pedigrees, reducing dependencies, and reducing the output size. That said, igv-reports may work well for your
  use-case.
