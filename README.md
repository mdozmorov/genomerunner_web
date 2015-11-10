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

