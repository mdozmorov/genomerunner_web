# GenomeRunner web server

[GenomeRunner web](http://www.intergarivegenomics.org/) is a tool for investigation of potential functional impact of sets of single nucleotide polymorphisms (SNPs) by considering their co-localization with fgenome annotation data (regulatory datasets). The philosophy behind GenomeRunner is that SNPs are not acting in isolation and may collectively alter regulatory/epigenomic features. Finding which regulatory features are affected may help to understand mechanisms of complex diseases from a holistic perspective.

GenomeRunner calculates enrichment p-values (Fisher's exact test) by evaluating whether a set of SNPs co-localizes with regulatory datasets more often that could happen by chance. GenomeRunner further performs three interpretation-oriented analyses:

1) Regulatory similarity analysis, aimed at identifying groups of SNP sets having similar functional impact;
2) Differential regulatory analysis, developed to identify regulatory signatures differentially enriched between the groups of SNPs;
3) Cell type regulatory enrichment analysis, designed to identify cell specificity of the regulatory enrichments. 

An example of GenomeRunner’s results can be found in the analysis of Sjogren’s syndrome GWAS ([Nature Genetics](http://www.nature.com/ng/journal/v45/n11/full/ng.2792.html) ), where it identified RFX5 transcription factor binding site as the most statistically significantly co-localized with the set of disease-associated SNPs.

Installation instructions

```bash
sudo ./setup.sh
```

will install all required packages. See the main documentation on
[https://mdozmorov.github.io/grdocs/index.html](https://mdozmorov.github.io/grdocs/index.html)

Use `R.GenomeRunner` to set up Shiny server, see [tips](https://github.com/mdozmorov/genome_runner/tree/shiny/web)

# Starting the server

```bash
gr-server -g hg19 -d /path/to/database/ -r /path/to/database/ -w 1 -p 8080
```
The "-g" argument specifies organism, the "-d" argument specifies path to the database of genome annotation data, the "-r" argument specifies path to the folder to output results, the optional "-w" argument specifies number of workers to run parallel jobs, the optional "-p" argument specifies at which port to start the server. Access the server locally at "http://localhost:8080"

# Command line usage

Prerequisites:
- Database of genome annotations, created with `gr-dbcreator`. Use full path. Example: /home/genomerunner/db_5.00_07-22-2015;
- Text file with full path(s) to BED file(s) containing genomic coordinates of SNPs of interest. Example: "foi.txt", where each line contains paths to BED files;
- Text file with full paths to genome annotation BED files. Example: "gf_encTfbs.txt", where each line contains paths to TFBS genome annotation BED files;
- The "background" file, a BED file containing genomic coordinates of all SNPs. This file is used to estimate enrichments that can happen by chance. The SNPs of interest should be a subset of this file.

The files containing categories of genome annotations can be created with
```bash
DIR=/home/genomerunner/db_5.00_07-22-2015
for file in `find $DIR/grsnp_db/hg19/ENCODE -maxdepth 1 -mindepth 1 -type d`; do
	gf=`basename $file`; echo $gf;
	find `$file -type f -name "*.bed.gz"` | sort > gf_enc$gf.txt;
done
```
The logic here is that the first `find` command finds top folders containing categories of genomic annotation data. The second `find` command outputs paths for each genome annotation data file into a category-specific file name, e.g., `gf_encTfbs.txt`.

The analysis can be run as, e.g.,
```bash
gr -g hg19 -d /path/to/database -r /path/to/results/folder -a fois.txt gf_encTfbs.txt /path/to/background.bed 
```
The optional "-a" argument specify if annotation analysis should be run.

This release is manually sunchronized with [the developmental branch](https://github.com/mdozmorov/genome_runner/tree/shiny) of GenomeRunner using `clear; diff -r --brief  genome_runner/ genomerunner_web/ | grep differ | grep -v .git | sort` and `clear; diff -r --brief  R.GenomeRunner/ genomerunner_web/R.GenomeRunner/ | grep differ | grep -v .git | sort`
