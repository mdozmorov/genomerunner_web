		<h3 dir="ltr" style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
			<span style="font-family:arial,helvetica,sans-serif;">Quick start</span></h3>
		<h4 dir="ltr" style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
			<span style="font-family:arial,helvetica,sans-serif;">Analysis of one SNP set</span></h4>
		<ol dir="ltr">
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span id="docs-internal-guid-1dace67e-ba8a-c2ae-4c69-a40c66dc5114">In step 1, expand "Paste data" field - this SNP set (genomic coordinates) will be selected for the analyses. Alternatively, upload a file with genomic coordinates in BED format, or rsIDs of your SNPs of interest.</span></span></li>
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span id="docs-internal-guid-1dace67e-ba8a-c2ae-4c69-a40c66dc5114">In step 2, expand "Upload custom background"-"Paste data" fields - this "background" SNP set will be used to estimate associations that can happen by chance. If uploading a file, ensure all SNPs selected in Step 1 are included in the background SNP set.</span></span></li>
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span style="white-space: pre-wrap;">In step 3, expand "Choose regulatory datasets" field, expand "ENCODE" category and select "Histone" subcategory.</span></span></li>
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span style="white-space: pre-wrap; line-height: 1.15;">Hit &quot;Submit job&quot; - the analysis should take a minute. The results should look like <a href="http://www.integrativegenomics.org/results_shiny?id=example1">http://www.integrativegenomics.org/results_shiny?id=example1</a></span></span></li>
		</ol>

		<h4 dir="ltr" style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
			<span style="font-family:arial,helvetica,sans-serif;">Analysis of multiple SNP sets</span></h4>
		<ol dir="ltr">
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span id="docs-internal-guid-1dace67e-ba8a-c2ae-4c69-a40c66dc5114">In step 1, click "<strong>diseases</strong>" button in the "Demo SNP sets" section;. This will automatically load several disease-associated SNP sets for analysis.</span></span></li>
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span id="docs-internal-guid-1dace67e-ba8a-c2ae-4c69-a40c66dc5114">In step 2, expand "Upload custom background"-"Paste data" fields - this "background" SNP set will be used to estimate associations that can happen by chance.</span></span></li>
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span style="white-space: pre-wrap;">In step 3, expand "Choose regulatory datasets" field, expand "ROADMAP" category and select "Histone_bPk" subcategory.</span></span></li>
			<li style="line-height: 1.15; margin-top: 0pt; margin-bottom: 0pt;">
				<span style="font-family:arial,helvetica,sans-serif;"><span style="white-space: pre-wrap; line-height: 1.15;">Hit &quot;Submit job&quot; - the analysis should take a minute. The results should look like <a href="http://www.integrativegenomics.org/results_shiny?id=example2">http://www.integrativegenomics.org/results_shiny?id=example2</a></span></span></li>
		</ol>
        <br>
                <span style="font-family:arial,helvetica,sans-serif;">To get more individual disease- or trait-specific sets of SNPs, to be uploaded in Step 1, explore <a href="https://github.com/mdozmorov/gwas2bed" target="_blank">https://github.com/mdozmorov/gwas2bed</a> repository. See BED file example <a href="https://github.com/mdozmorov/gwas2bed/blob/master/gwasCatalog/bed/more15/systemic_lupus_erythematosus.bed" target="_blank">systemic_lupus_erythematosus.bed</a></span></p>
        <br>
                <span style="font-family:arial,helvetica,sans-serif;">To get more regulatory/epigenomic datasets, to be uploaded in Step 3, advanced features, explore <a href="https://github.com/mdozmorov/gwas2bed/tree/master/trynka-raychaudhuri" target="_blank">https://github.com/mdozmorov/gwas2bed/tree/master/trynka-raychaudhuri</a> repository. See BED file example <a href="https://github.com/mdozmorov/gwas2bed/blob/master/trynka-raychaudhuri/H3K4me3_CD34_Cultured_Cells.bed" target="-blank">H3K4me3_CD34_Cultured_Cells.bed</a></span></p>
