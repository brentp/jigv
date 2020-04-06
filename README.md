[igv.js](https://github.com/igvteam/igv.js) is a great way to view aligments and other genomic data. 
It requires that the files are hosted on a server, like apache or nginx and it requires writing html and
javascript.

In a single binary, `jigv` provides a server and some default configuration,
javascript, and HTML so you can do:

```
jigv --open-browser --region chr1:34566-34999 *.bam /path/to/some.cram my.vcf.gz
```
With that, a server will start and a browser will open an igv.js viewer for your requested files, 
similiar to how you'd do with the java version of IGV. Without the `--open-browser` option, you can
go to *http://$server:5001/* (where `$server` is often `localhost`).

# installation

grab a static linux binary from [releases](https://github.com/brentp/jigv/releases/latest)

# notes and features

#### navigation

`jigv` supports navigation along variants in a VCF file using the left and right arrow keys. it uses the first VCF track
as the one to navigate through.

#### labels 

by default a file like `/path/to/some.bam` will have a name of `some.bam` to change that send `"/path/to/some.bam#label"`
where anything after the `#` is used as the label.

#### remote servers 

This can be very useful on a `server` even if it's not open for browsing non-standard ports: in a screen start jigv with the desired files. 
Then from your remote machine connect with ssh tunneling:

```
# in a screen/tmux session
server: jigv "/path/to/kid.bam#proband" "/path/to/mom.bam#mom" "/path/to/dad.bam#dad" "/path/to/joint.vcf.gz"

# set up the tunnel on local machine (replace $server with the hostname you would normally ssh)
local: ssh -N -L localhost:5001:localhost:5001 $server
```

Then browse localhost:5001 on your local machine.


#### automated screenshots

given a server running at: localhost:5001, a PNG screenshot for `region=chr5:1410040-1412784` can be created with

```
google-chrome --window-size=1200,1200  --virtual-time-budget=10000 \
   --headless --run-all-compositor-stages-before-draw \
   --screenshot=${region/:/-}.png "http://localhost:5001/#$region"
```

This can be scripted with your favorite language for a set of regions.

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

# limitations

+ this is likely insecure in many ways.
+ there will soon be a way to customize the options and javascript (but this probably covers 85% of use-cases as-is).
+ if you have some custom javascript used with igv.js, that is generally useful, please open an issue so I can add it.
+ not all file types are supported

# examples

while this is most helpful when you don't already have a server, you can use it to quickly view files
that are public:
```
url=https://s3.amazonaws.com/1000genomes/phase3/data
jigv -g hg19 \
    -r "chr5:1,278,267-1,281,011" \
    -o "$url/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam#kid" \
       "$url/NA18637/exome_alignment/NA18637.mapped.ILLUMINA.bwa.CHB.exome.20121211.bam#mom"
```

