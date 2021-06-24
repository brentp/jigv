[igv.js](https://github.com/igvteam/igv.js) is a great way to view aligments and other genomic data. 
It requires that the files are hosted on a server with access to the original data.
javascript.

`jigv` encodes all variants, alignments, and annotations into a single HTML page (or as separate javascript file for
each variant) that you can to collaborators who don't have access to the cluster where your data is stored.

The resulting file is very fast to navigate; the left/right arrow keys advance to next/previous variants of interest.

Usage looks like:

```
jigv \
    --sample HG01053 \     # the sample of interest drawn in top panel
    --ped trio.ped \       # jigv will use this to also show parents and sibs of --sample
    --sites dn.vcf.gz \    # e.g. candidate de novo varants
    --fasta $reference_fasta \
    --annotation hg38.refGene.bed.gz \ # see: https://github.com/brentp/jigv/wiki/bed12
    --annotation LCR-hs38.bed.gz \     # specify as many of these as needed.
    /path/to/*.cram > denovos.html
```
With that, `denovos.html` will contain **all genomic data around variants of interest** embedded within it.

With this, we are able to encode 924 candidate *de novo* variants into a 31MB html file that includes
alignments for the proband, mom, and dad.

See [Examples](#Examples) for more

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
  jigv [options] [xams ...]

Arguments:
  [xams ...]       indexed bam or cram files for relevant samples. read-groups must match samples in vcf.

Options:
  --sample=SAMPLE            sample-id for proband or sample of interest (default is first vcf sample)
  --sites=SITES              VCF or BED containing variants of interest for --sample. If this contains ':', then it's used as a single region and the first bam/cram given is the sample of interest. If it ends with '.bed' or '.bed.gz' it's assumed to be BED format.
  -g, --genome-build=GENOME_BUILD
                             genome build (e.g. hg19, mm10, dm6, etc, from https://s3.amazonaws.com/igv.org.genomes/genomes.json).  If this is specified then the page will request fasta, ideogram and gene data from a server.
  --cytoband=CYTOBAND        optional path to cytoband/ideogram file
  --annotation=ANNOTATION    path to additional bed or vcf file to be added as a track; may be specified multiple times
  --ped=PED                  pedigree file used to find relations for --sample
  --prefix=PREFIX            if specified, encoded data for each region is written to it's own js file and no html is generated. this is a path prefix like: 'jigv_encoded/HG002/' where  where and `$site.js` will be added by jigv
  --prefix-raw               by default if --prefix is specified, the data is written to a javascript variable and includes the 'data:base64' prefix. if this option is also specified, then the raw base64 encoded data is written to the file.
  --fasta=FASTA              path to indexed fasta file; required for cram files
  --flank=FLANK              bases on either side of the variant or region to show (default: 100) (default: 100)
```
### prefix

By default, jigv will embed all data for all sites into a single html file. This can scale up to about 2000 variants
for a trio at which point the single file becomes quite large (~60MB).
Using `--prefix`, it's possible to send each encoded site to a separate file. These can be loaded on demand from a
server or sent to the user and loaded as a script on demand from the local file system.

Files written with `--prefix` but without `--prefix-raw' can be loaded from the local file-sytem of the user (or from a server depending on if they are viewing as `file:///...` on their local system or on a server. The code for that looks like this:

```Javascript

    var jigv_data; // define this here, but it is set when the script is loaded.
    var s = document.createElement("script")
    s.type = "text/javascript"
    s.src = "$sample/chr1-1746431-A-T.jigv.js" // path to whichever file made by jigv.
    s.onload = function(a, b) {
        jigv_data["showChromosomeWidget"] = false
        jigv_data["search"] = false
        jigv_data["queryParametersSupported"] = true
        jigv_data["showCursorTrackingGuide"] = true

        var div = document.getElementById("igv")
        meta_options = {"sessionURL": jigv_data }
        igv.createBrowser(div, meta_options)
    }

```

note that jigv writes the data as `jigv_data` so once the file is loaded, it can be accessed via that variable name.

# limitations

+ if you have some custom javascript used with igv.js, that is generally useful, please open an issue so I can add it.
+ not all file types are supported
+ the sites file must contain only the variants to plot; there is no filtering in `jigv`


## See Also

+ [igv-reports](https://github.com/igvteam/igv-reports) by the IGV team does this as well. `jigv` adds extra features for
  working with pedigrees, reducing dependencies, and reducing the output size. That said, igv-reports may work well for your
  use-case.

# examples

1. a set of bam files given a set of **regions in BED format**:

```
jigv \
    --fasta $fasta \
    --sites $bed \
    /path/to/*.bam > regions.html
```

This will use the first 8 bam files. Beyond that number the files are too large.


2. a set of **de novo** variants for a trio:

```
jigv \
    --sample $proband \     # the sample of interest drawn in top panel
    --ped trio.ped \        # jigv will use this to also show parents and sibs of --sample
    --sites dn.vcf.gz \     # file with candidate de novo varants
    --fasta $fasta \
    --annotation hg38.refGene.bed.gz \ # see: https://github.com/brentp/jigv/wiki/bed12
    --annotation LCR-hs38.bed.gz     \ # specify as many of these as needed.
    /path/to/*.cram > denovos.html
```

3. a **single region**

```
jigv \
    --sample $proband \     # the sample of interest drawn in top panel
    --ped trio.ped \        # jigv will use this to also show parents and sibs of --sample
    --sites "chr1:3453454-3453554" \  # single region.
    --fasta $fasta \
    --annotation hg38.refGene.bed.gz \ # see: https://github.com/brentp/jigv/wiki/bed12
    --annotation LCR-hs38.bed.gz     \ # specify as many of these as needed.
    /path/to/*.cram > denovos.html
```
