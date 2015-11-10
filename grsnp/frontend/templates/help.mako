		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">How GenomeRunner Web can help me?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">GenomeRunner Web helps to interpret the regulatory effect of SNPs by identifying functional elements most statistically significantly co-localized with the SNPs of interest (see <a href="#enrichment">Enrichment analysis</a>). These regulatory elements may provide understanding which mechanisms (and in which cell type) may be affected by the SNPs of interest.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">If one analyzes three or more sets of SNPs, GenomeRunner Web visualizes their epigenomic similarity (see <a href="#episimilarity">Epigenomic Similarity analysis</a>). This information may be used to classify the sets of SNPs by similarity of their effect upon regulatory landscape.</span></div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">How my input data should look like?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Tab-separated text files with genomic coordinates of SNPs of interest in <a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1" target="_blank">.BED format</a>. As a bare minimun, chromosome, start, and end coordinates should be provided. One can upload a file, or copy-paste tab-separated coordinates. Multiple file upload/analysis is possible (e.g., patient-specific sets of SNPs). </span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;"><strong>Note:</strong> a set of SNPs should contain at least 5 SNPs to be eligible for the analyses</span></div>
		<div>
			<span style="font-family: arial, helvetica, sans-serif;"><strong>Note:</strong> coordinates should be 0-bases. The <em>end </em>coordinate should equal <em>start + 1</em>.</span></div>
		<div>
			&nbsp;</div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Several pre-defined sets of SNPs are available.&nbsp;</span><span style="font-family: arial, helvetica, sans-serif;">&nbsp;For&nbsp;</span><em style="font-family: arial, helvetica, sans-serif;">Homo Sapiens</em><span style="font-family: arial, helvetica, sans-serif;">&nbsp;these include:</span></div>
		<div>
			&nbsp;</div>
		<div>
			<table border="1" cellpadding="1" cellspacing="1" style="width: 1000px;">
				<tbody>
					<tr>
						<td style="text-align: center;">
							<strong><span style="font-family:arial,helvetica,sans-serif;">Pre-defined sets of SNPs</span></strong></td>
						<td style="text-align: center;">
							<strong><span style="font-family:arial,helvetica,sans-serif;">When to use</span></strong></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">GWAS_vs_DGV - selected disease-associated SNP sets from <a href="http://www.genome.gov/26525384" target="_blank">GWAScatalog</a>.</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">Used for <a href="./demo">Demo</a> run, to be run for enrichment vs. structural variants (dgvVariants)</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">GWAS_vs_tfbsEncode</span><span style="font-family: arial, helvetica, sans-serif;">&nbsp;- selected disease-associated SNP sets from&nbsp;</span><a href="http://www.genome.gov/26525384" style="font-family: arial, helvetica, sans-serif;" target="_blank">GWAScatalog</a><span style="font-family: arial, helvetica, sans-serif;">.</span></td>
						<td>
							<span style="font-family: arial, helvetica, sans-serif;">Used for&nbsp;</span><a href="./demo" style="font-family: arial, helvetica, sans-serif;">Demo</a><span style="font-family: arial, helvetica, sans-serif;">&nbsp;run, to be run for enrichment vs. 161 experimentally obtained transcription factor binding sites (tfbsEncode)</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">snp138Rand - randomly selected sets of SNPs from the <a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=snp138&amp;hgTracksConfigPage=configure" target="_blank">snp138 </a>database.</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">Use with snp138+ background and any (category of) genome annotation elements, to ensure random SNPs do not show significant associations</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family: arial, helvetica, sans-serif;">gwasRand - randomly selected sets of SNPs from&nbsp;&nbsp;</span><a href="http://www.genome.gov/26525384" style="font-family: arial, helvetica, sans-serif;" target="_blank">GWAScatalog</a><span style="font-family: arial, helvetica, sans-serif;">.</span></td>
						<td>
							<span style="font-family: arial, helvetica, sans-serif;">Use with gwascatalog background and any (category of) genome annotation elements, to ensure random SNPs do not show significant associations</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">GWASmore15 - all&nbsp;</span><span style="font-family: arial, helvetica, sans-serif;">disease- or trait-associated SNP sets from&nbsp;</span><a href="http://www.genome.gov/26525384" style="font-family: arial, helvetica, sans-serif;" target="_blank">GWAScatalog</a><font face="arial, helvetica, sans-serif">&nbsp;having more than 15 SNPs in a set</font></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">Used for explorartory analyses vs. any (category of) genome annotation elements</span></td>
					</tr>
				</tbody>
			</table>
		</div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">What is &ldquo;background&rdquo;?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">The background is the &ldquo;universe&rdquo; of all SNPs assessed in a study, from which SNPs of interest came from. Several pre-defined backgrounds are provided,&nbsp;</span><span style="font-family: arial, helvetica, sans-serif;">&nbsp;for&nbsp;</span><em style="font-family: arial, helvetica, sans-serif;">Homo Sapiens</em><span style="font-family: arial, helvetica, sans-serif;">&nbsp;these include:</span></div>
		<div>
			&nbsp;</div>
		<div>
			<table border="1" cellpadding="1" cellspacing="1" style="width: 1000px;">
				<tbody>
					<tr>
						<td style="text-align: center;">
							<strong><span style="font-family:arial,helvetica,sans-serif;">Pre-defined background</span></strong></td>
						<td style="text-align: center;">
							<strong><span style="font-family:arial,helvetica,sans-serif;">When to use</span></strong></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">snp138+ (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=snp138&amp;hgTracksConfigPage=configure" target="_blank">All&nbsp;Simple Nucleotide Polymorphisms (dbSNP 138)</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">For SNPs from whole-genome GWA studies</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">snp138Common+ (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=snp138Common&amp;hgTracksConfigPage=configure" target="_blank">Simple Nucleotide Polymorphisms (dbSNP 138) Found in &gt;= 1% of Samples</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">For SNPs from studies where rare variants were ignored</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">gwascatalog+ (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=gwasCatalog&amp;hgTracksConfigPage=configure" target="_blank">NHGRI Catalog of Published Genome-Wide Association Studies</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">For demo testing, to observe regulatory associations of disease-specific sets of SNPs, as compared with randomly selected SNPs from all GWAScatalog</span></td>
					</tr>
				</tbody>
			</table>
		</div>
		<div>
			&nbsp;</div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">For a GWAS, the background is likely to be all SNPs (snp138+). For a study using microarrays, the background would contain coordinates of all SNPs on the microarray - upload or copy/paste them.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;"><strong>Note:</strong> SNPs of interest should be a subset of the background SNPs. If some SNPs of interest do not overlap the background, a non-critical error is issued. Create a custom background that includes all SNPs of interest. An excellent set of <a href="https://code.google.com/p/bedtools/" target="_blank">bedtools </a>made by Aaron Quinlan is recommended.</span></div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">What are &ldquo;genome annotation features&rdquo;?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Genome annotations are discrete regions potentially having functional/regulatory properties.&nbsp;</span></div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">There are too many genome annotation features! What to choose?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">The genome annotation features are organized by categories mirrored from the UCSC genome browser. Using search box and/or checkboxes in the TreeView control, one can select one or more genome annotation categories. Clicking on a genome annotation&rsquo; name will bring up the description, if available.</span></div>
		<div>
			&nbsp;</div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">The ENCODE data are organized by the source/data type, tiers (quality), and by cell types.</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;"><strong>Hint:</strong> Several well-known/specially processes genome annotation features sets are brought forward as &ldquo;default genome annotation features&rdquo;. For <em>Homo Sapiens</em> these include:</span></div>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">&nbsp;</span></div>
		<div>
			<table border="1" cellpadding="1" cellspacing="1" style="width: 1000px;">
				<thead>
					<tr>
						<th scope="col">
							<span style="font-family:arial,helvetica,sans-serif;">Genome annotation category</span></th>
						<th scope="col">
							<span style="font-family:arial,helvetica,sans-serif;">Experimental question:&nbsp;Are the SNPs of interest...</span></th>
					</tr>
				</thead>
				<tbody>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">H3K4me3 (Cell/tissue-specific H3K4me3 marks, from <a href="http://immunogenomics.hms.harvard.edu/PDFs/2013%20-%20EpiGWAS.pdf">Trynka-Raychaudhuri Nature Genetics paper</a>)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... enriched in H3K4me3 active promoters and encancers mark? In which cell/tissue type?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">coriellVariants (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=coriellDelDup&amp;hgTracksConfigPage=configure" target="_blank">Coriell Cell Line Copy Number Variants</a>, split by cell types)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... enriched in CNVs, and in which cell type?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">dgvVariants (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=dgvPlus&amp;hgTracksConfigPage=configure" target="_blank">Database of Genomic Variants: Structural Variation (CNV, Inversion, In/del)</a>, split by variant type)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... enriched in CNVs, or other types of structural variations?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">repeats (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=rmsk&amp;hgTracksConfigPage=configure" target="_blank">Repeating Elements by RepeatMasker</a>, split by repeat class)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... enriched in regions of low complexity, and in which type?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">gwasCatalog (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=gwasCatalog&amp;hgTracksConfigPage=configure" target="_blank">NHGRI Catalog of Published Genome-Wide Association Studies</a>, split by disease/trait types)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... enriched in known disease-specific SNPs?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">altSplicing (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=knownAlt&amp;hgTracksConfigPage=configure" target="_blank">Alternative Splicing, Alternative Promoter and Similar Events in UCSC Genes</a>, split by splicing type)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... potentially disrupt a specific type of alternative spliced regions?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">ncRNAs (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=wgRna&amp;hgTracksConfigPage=configure" target="_blank">C/D and H/ACA Box snoRNAs, scaRNAs, and microRNAs from snoRNABase and miRBase</a>, split by ncRNA type)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... associated with a class of non-coding elements?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">tfbsEncode (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=wgEncodeRegTfbsClusteredV3&amp;hgTracksConfigPage=configure" target="_blank">Transcription Factor ChIP-seq Clusters V3 (161 targets, 189 antibodies) from ENCODE</a>, split by TFBS name)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... potentially disrupt a specific experimentally defined transcription factor binding site?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">tfbsConserved (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=tfbsConsSites&amp;hgTracksConfigPage=configure">HMR Conserved Transcription Factor Binding Sites</a>, split by TFBS name)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... potentially disrupt a specific&nbsp;computationally predicted&nbsp;transcription factor binding site?</span></td>
					</tr>
					<tr>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">chromStates (<a href="http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=354750109&amp;g=wgEncodeBroadHmm&amp;hgTracksConfigPage=configure" target="_blank">Chromatin State Segmentation by HMM from ENCODE/Broad</a>, Gm12878 cell line, split by chromatin state type)</span></td>
						<td>
							<span style="font-family:arial,helvetica,sans-serif;">... preferentially located in certain chromatin states?</span></td>
					</tr>
				</tbody>
			</table>
			<p>
				<span style="font-family: arial, helvetica, sans-serif;">Examples of what to choose:</span></p>
		</div>
		<ul>
			<li>
				<span style="font-family:arial,helvetica,sans-serif;">Select &ldquo;tfbsEncode&rdquo; category to get an answer whether the SNPs of interest are enriched in any of the 161 transcription factor binding sites identified by ChIP-seq.</span></li>
			<li>
				<span style="font-family:arial,helvetica,sans-serif;">Select &ldquo;ENCODE/BroadHistone/Tier1/Gm12878&rdquo; category to get an insight whether the SNPs of interest are enriched in histone marks in B lymphocytes.</span></li>
			<li>
				<span style="font-family:arial,helvetica,sans-serif;">Select &ldquo;genes&rdquo; category to answer a question, whether the SNPs of interest are enriched in genes/exons.</span></li>
		</ul>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;"><a name="enrichment"></a>What is Enrichment analysis?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Enrichment analysis answers the question whether a set of SNPs of interest collectively enriched or depleted in regulatory regions, as compared with randomly selected set of SNPs. GenomeRunner Web performs this analysis versus each selected (group of) genome annotation features, and prioritized over- and underrepresented associations by p-value. The p-values are calculated using Fisher&#39;s exact, or Chi-squared tests.</span></div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;"><a name="episimilarity"></a>What is Epigenomic similarity analysis?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Epigenomic similarity analysis visualizes similarity among enrichment profiles for three or more sets of SNPs of interest. It answers the question whether different sets of SNPs are enriched in similar epigenomic elements, hence, may affect similar regulatory networks. This analysis may be useful when comparing sets of SNPs from multiple GWA studies, or patient-specific variomes. It is particularly useful for the analysis of patient-specific rare variants, to observe whether patients can be classified by similar regulatory elements affected by their SNPs.</span></div>
		<h3>
			<span style="font-family: arial, helvetica, sans-serif;">What else can I do with the results?</span></h3>
		<div>
			<span style="font-family:arial,helvetica,sans-serif;">Head to&nbsp;<a href="https://github.com/mdozmorov/R.genomerunner" target="_blank">https://github.com/mdozmorov/R.genomerunner</a> for further ideas and the R code to explore and visualize GenomeRunner&#39;s results.&nbsp;</span></div>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">How to cite GenomeRunner?</span></h3>
		<p>
			<span style="font-family:arial,helvetica,sans-serif;">If you find GenomeRunner useful, please, cite:</span></p>
		<p>
			<span style="font-family:arial,helvetica,sans-serif;"><a href="http://www.ncbi.nlm.nih.gov/pubmed/22155868">Dozmorov MG, Cara LR, Giles CB, Wren JD. GenomeRunner: automating genome exploration.</a><br />
			Bioinformatics. 2012 Feb 1;28(3):419-20.&nbsp;<a href="http://bioinformatics.oxfordjournals.org/content/28/3/419">doi:10.1093/bioinformatics/btr666.</a><br />
			PubMed PMID: 22155868; PubMed Central PMCID: PMC3268239</span></p>
		<p>
			<span style="font-family:arial,helvetica,sans-serif;">GenomeRunner Web publication is currently submitted. Please, check back for the updates</span></p>
		<h3>
			<span style="font-family:arial,helvetica,sans-serif;">How GenomeRunner Web is different from previously published version?</span></h3>
		<p>
			<span style="font-family: arial, helvetica, sans-serif;">The original version of GenomeRunner, available at <a href="http://sourceforge.net/projects/genomerunner/" target="_blank">SourceForge</a>, was designed as an &quot;all purpose&quot; tool. It has several disadvantages, such as steep learning curve, download of large databases, restriction to Windows platform, lacking any visualization capabilities. GenomeRunner Web addresses these issues - it is:</span></p>
		<ul>
			<li>
				<span style="font-family: arial, helvetica, sans-serif;">Platform-independent, runs in a browser</span></li>
			<li>
				<span style="font-family: arial, helvetica, sans-serif;">Simple, due to careful interface design</span></li>
			<li>
				<span style="font-family: arial, helvetica, sans-serif;">Visual, the results are presented as interactive heatmaps</span></li>
			<li>
				<span style="font-family: arial, helvetica, sans-serif;">Unique, a novel idea of &#39;<a href="#episimilarity">epigenomic similarity</a>&#39; analysis is fully implemented and the results are visualized</span></li>
			<li>
				<span style="font-family: arial, helvetica, sans-serif;">Downloads-free, we host and maintain genome annotation databases</span></li>
			<li>
				<span style="font-family: arial, helvetica, sans-serif;">Supports local installation and command line usage in 3 easy steps (work in progress)</span></li>
		</ul>
		<h3>
			<span style="font-family: arial, helvetica, sans-serif;">Will the <a href="http://sourceforge.net/projects/genomerunner/" target="_blank">old </a>version be updated?</span></h3>
		<p>
			<span style="font-family: arial, helvetica, sans-serif;">The old version is in stable alpha status now. However, the databases will be updated periodically, to incorporate the growing amount of multi-organism genome annotation d</span><span style="font-family: arial, helvetica, sans-serif;">ata. Currently, the updated and reorganized databases for hg19 human genome assembly (December 2013) are available for download (Warning: large files):</span></p>
		<p>
			<span style="color:#0000ff;"><span style="font-family:arial,helvetica,sans-serif;">Files temporarily unavailable due to ongoing data migration to a new server. Contact the developer if you need the database files immediately</span></span></p>
		<p>
			<span style="font-family: arial, helvetica, sans-serif;"><a href="static/hg19encode.sqlite">hg19encode.sqlite</a>&nbsp;(87Gb) - the ENCODE project data from the UCSC genome browser (table names starting with <em>wgEncode...</em>). The data are organized as &quot;Source-type/Tier/Cell&quot; hierarchy, to better reflect the UCSC theme and make navigation/selection of groups of epigenomic elements easier;</span></p>
		<p>
			<span style="font-family: arial, helvetica, sans-serif;"><a href="static/hg19new.sqlite">hg19new.sqlite</a>&nbsp;(127Gb) - non-ENCODE genome annotation data from the UCSC genome browser. Manually organized into tiers.</span></p>
		<h3>
			<span style="font-family: arial, helvetica, sans-serif;">I still have questions/suggestions/bug report...</span></h3>
		<p>
			<font face="arial, helvetica, sans-serif">Please, contact&nbsp;</font><br />
			<img height="20%" src="static/images/e-mail.png" width="20%" /></p>
		<p>
			<span style="font-family:arial,helvetica,sans-serif;">with anything you may want to mention or discuss and expect a prompt response. If you have a run crashed, please, send the run ID (from the URL bar, like&nbsp;<em>result?id=example</em>) for troubleshooting.</span></p>
