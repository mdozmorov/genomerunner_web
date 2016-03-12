<body bgcolor="#A8D5FF" style="width: 1000px">

	<div class="topbar">
		<div class="topbar-inner">
			<div class="container-fluid">
				<a class="brand" href="./" style="font-family:Futura,'Helvetica Neue', Arial">
					<img width="200px" src="static/images/GRLogo.png" />
				</a>
				<ul class="pull-right">
					<li><a href="./overview">Overview</a></li>
					<!--<li><a href="./news">News</a></li>-->
					<li><a href="./demo">Quick start</a></li>
					<!-- <li><a href="./cite">How to Cite</a></li>
					<li><a href="http://sourceforge.net/projects/genomerunner/">GenomeRunner on SourceForge</a></li>
					<li><a href="./roadmap">Roadmap</a></li> -->
					<li><a href="https://github.com/mdozmorov/genomerunner_web/wiki" target="_blank">Help</a></li>
<!-- <div>
					<li><a href="./google">Google Group</a></li>
					<li><img width="30px" src="static/new-icon.jpg" alt="New: GenomeRunnerSNP Google Groups" /></li>
				</div> -->
			</ul>
			<!-- <img height="40px" src="static/images/logo-reversed-small.jpg" align="right" /> -->
		</div>
	</div>
</div>
<form id="frmQuery" name="frmQuery" action="query" method="post" enctype="multipart/form-data">
	<div id="content">
		<div class="well" style="margin-top: -15px; padding: 0px">			
			<h3>Select Database Version:</h3>
			${database_versions}		
			<h3>GenomeRunner: Functional interpretation of SNPs within regulatory/epigenomic context</h3>
			<p>
				<span style="font-size: 16px;"><span style="font-family:arial,helvetica,sans-serif;">GenomeRunner is a tool for functional enrichment analysis of SNP sets within regulatory/epigenomic context. The philosophy behind GenomeRunner is that SNPs are not acting in isolation and may collectively alter regulatory/epigenomic features. Finding which regulatory features are affected may help to understand  mechanisms of complex diseases from a holistic perspective.</span></p>
				<p>
				<br />
				<span style="font-size:16px;"><span style="font-family:arial,helvetica,sans-serif;">GenomeRunner performs regulatory enrichment/annotation analyses, differential regulatory analysis, and cell type-specific enrichment analysis. The downloadable results are visualized as interactive heatmaps and tables <a href="results_shiny?id=example1">(Example 1, single SNP set analysis)</a>, <a href="results_shiny?id=example2">(Example 2, multiple SNP sets analysis)</a>.</span></span></p>
					<p>
					</p>
					</div>
					<div class="well">
						<h3>1. Select sets of SNPs of interest <img class="helptooltip" title="Upload file(s) with rsIDs of the SNPs of interest, or their genomic coordinates in BED format . At least 5 SNPs per file are required. Note: Avoid special characters and extra dots in the file names. Do not use SNP IDs other than rsIDs." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</h3>
						
<table border="0">
	<tr>
		<td>
				<h4 style="float:left;">Files:</h4><input type="file" id="inputbedfile" style="margin:5px" name="bed_file" multiple="multiple"/>
				<a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">What should the data in BED format look like?</a>
		</td>
		<td>
				<h4 style="float: left;margin-top: 4px;margin-right: 9px;">Organism:</h4>
				${paths.org_as_html(id="org")} 
				<img class="helptooltip" title="Select organism-specific genome assembly" style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/> 
		</td>
	</tr>
	<tr>
		<td>
			<h4 style="float: left;margin-top: 10px;margin-right: 9px;">Demo SNP sets (click to select): </h4>
			<div class="btn-group" id="btngroup_demo_fois" data-toggle="buttons-radio" data-toggle-name="demo_fois">
				<input type="hidden" style="margin-top: -3px;" name="demo_fois"/>
				${demo_snps}
			</div>
		<td>
	</tr>
</table>

						<div id="accfoi" class="accordion" style="padding-bottom: 0em;list-style:none;margin-top: 20px">       
							<h3  id="accordionheader"><a href="#" style="font-size:120%">Paste data</a><img class="helptooltip" title="Paste a list of rsIDs of SNPs of interest, or tab-separated genomic coordinates (recommended) in .BED format (no headers)." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/></h3>
							<div>          
								<table>
									<tr>
										<td>
											<textarea id="inputbeddata" rows=10 cols=95 style="margin:10px;margin-bottom:0px;" name="bed_data" wrap="off" disabled>
											</textarea>
										</td>
									</tr>
								</table>
							</div>
						</div>
					</div>	
					<div class="well">
						<h3 style="float:left">2. Define the background: ${default_background}<img class="helptooltip" title="By default, all common SNPs are used as a 'background' of all SNPs, used to calculate the probability of the SNPs of interest being enriched in regulatory datasets. The SNPs of interest should be a subset of the background, or the p-values may be incorrect. The default background, all common SNPs from the latest organism-specific database, is suitable when a genome-wide study was performed. When a microarray was used for SNPs profiling, it is advisable to upload all SNPs on that array as a background." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
						</h3>
						<!-- <input type="checkbox" style="font-size:120%;margin-top:1em" name="run_random">Run randomization test</input> -->
						<b style=" font-size:120%; margin-left: 10px"></b>
						<div class="accordion" style="padding-bottom: 1em;list-style:none;margin-top:2em">        
							<h3  id="accordionheader"><a href="#" style="font-size:120%">Upload/paste custom background
								<img class="helptooltip" title="Upload, or paste a list of rsIDs, or genomic coordinates in BED format, of all SNPs evaluated in a study." width=25 height=25 src="static/images/help-icon.png" alt="help"/>
							</a></h3>
							<div style="height=100px">       
								<h4 style="float:left;">File:</h4>
								<input type="file" id="inputbackgroundfile" style="margin:5px" name="background_file" />
								<a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">What should data in BED format look like?</a>
								<div id="accback" class="accordion" style="padding-bottom: 1em;list-style:none;margin-top: 20px">        
									<h3  id="accordionheader"><a href="#" style="font-size:120%">Paste data</a><img class="helptooltip" title="Paste a list of rsIDs of SNPs of interest, or tab-separated genomic coordinates (recommended) in .BED format (no headers)." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/></h3>					        
									<table style="margin-bottom:0px; padding-bottom:0px">
										<tr>
											<td>
												<textarea id="inputbackgrounddata" rows=5 cols=95 style="margin:10px" name="background_data" wrap="off" disabled>
												</textarea>
											</td>
										</tr>
									</table>
								</div>
							</div>
						</div>
					</div>
					<div class="well" style="padding-bottom:0em">
						<h3 style="margin-right: 10px;margin-top: 1px;">3. Select regulatory/epigenomic datasets<img class="helptooltip" title="Select and/or upload regulatory/epigenomic datasets (aka functional annotations) used to test for enrichments in the SNPs of interest. Each regulatory/epigenomic dataset contains genomic coordinates (BED format) of regions annotated as carrying functional/regulatory potential or having a biological property." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/></h3>
						${custom_gfs}
						<br>
						<div id="accordGFS" class="accordion" style="padding-bottom: 1em;list-style:none; margin-top: 0.5em;" onClick="renderCheckBoxTree()"> 
							<h3  id="accordionheader"><a id='gfselheader' href="#" style="font-size:120%; height: 100%">Choose regulatory/epigenomic datasets</a><img class="helptooltip" title="Use checkboxes to select any (category of) regulatory/epigenomic datasets from a Database/Category/Subcategory hierarchy." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/></h3>
							<div >								
								<div id="grfdroplist" style="display: table;">	
									<div id="divCheckBox" style="width: 70%; margin:15px; display: table-cell; verticle-align: top; visibility: hidden">
										<div style="display: table-row;">
											<div id="checkbuttons" style="margin-right: 1em;">
												<a class="btn" style="margin-top: 9px;"onClick="$('#jstree_gfs').jstree('open_all');">Expand All</a>
												<a class="btn" style="margin-top: 9px;" onClick="$('#jstree_gfs').jstree('close_all');">Collapse All</a>
												<a class="btn" style="margin-top: 9px;" onClick="$('#jstree_gfs').jstree('check_all');">Select All</a>
												<a class="btn" style="margin-top: 9px" onClick="$('#jstree_gfs').jstree('uncheck_all');">Deselect all</a>
												<!-- <a class="btn" style="margin-top: 9px;" id="descriptions">Track Descriptions</a> -->
											</div>
										</div>												
											<label>Search regulatory/epigenomic datasets<img class="helptooltip" title="Use fuzzy search to find regulatory/epigenomic datasets; use Select/Unselect buttons to mark them for the analysis." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/></label>
											<input style="margin-top:1.5em" id="txt_gfs_search" class='input' type="text"></input>
											<a class="btn" id="treeSelect" onClick="treeviewSelectSearchedClick()">Select</a>
											<a class="btn" id="treeSelect" onClick="treeviewDeselectSearchedClick()">Unselect</a>
											<div id="jstree_gfs" style="height:200px;overflow:auto"></div>
									</div>
								</div>
							</div>
						</div>
					</div>	

					<div class="well">
						<div>
							<button id="btnSubmit" class="btn btn-primary" type="submit" style="margin:1em">Submit job</button>									
							<input type="checkbox" id="disclaimer" checked="checked" style="margin:1em">The results of GenomeRunner analyses may be used for research purposes only.</input><br/>
							<h3 id="upmessage" style="margin-left:11em;margin-bottom:3em;visibility:hidden;">Uploading files. Please do not refresh the page.</h3>
							<br>
						</div>
						<div class="accordion" style="margin-top:-46px;margin-bottom:18px">
							<h3 id='lblAdvancedFeatures'><a id='gfselheader' href="#" style="font-size:120%; height: 100%">Advanced Features</a></h3>

							<div>
								<div id="accordGFSfile" class="accordion" style="padding-bottom: 0.5em;list-style:none; margin-top: 1em;" ><h3  id="accordionheader"><a id='gffileheader' href="#" style="font-size:120%; height: 100%">Upload custom regulatory/epigenomic datasets:<img class="helptooltip" title="If needed, provide custom regulatory/epigenomic datasets, to be used for the analysis. Each custom regulatory/epigenomic dataset should contains genomic coordinates (BED format) of regulatory regions of interest." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/></a></h3>
							<ul>
								<li id="list-bedbackground">					
									<h4 style="float:left;">Files:</h4>
									<input type="file" id="inputgenomicfeaturefile" style="margin:5px" name="genomicfeature_file" multiple="multiple"/>
									<a href="http://genome.ucsc.edu/FAQ/FAQformat.html#format1">What should the data in BED format look like?</a>
								</li>				
							</ul>
						</div>

								<label>Percent score threshold: ${pct_scores}</label>
								<img class="helptooltip" title="Increasing this number filters out more low-level signal in the regulatory datasets. If a regulatory dataset does not have a signal value, this setting is ignored." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								<input type="checkbox" style="font-size:150%;margin-top:1em;margin-left:3em" name="run_annot">Run annotation analysis</input><img class="helptooltip" title="Annotate each SNP in each file by the number of overlaps with the selected regulatory datasets. Increases run time." style="position: relative;top: 6px;" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								<br>
								<label style="margin-left:10px;">Statistical test: </label>
								<select name="stat_test">
									<option value="chisquare" selected>Chi-square</option>
									<option value="binomial">Binomial</option>
									<option value="montecarlo">Monte Carlo</option>
								</select>
								<img class="helptooltip" title="Select test to obtain enrichment p-values.  Chi-square test recommended." style="position: relative;top: 6px" width=25 height=25 src="static/images/help-icon.png" alt="help"/>
								<label style="margin-left:10px;visibility:hidden" id="lbl_num_mc">Number of Monte Carlo simulations: </label>
								<select style="visibility:hidden" name="num_mc">
									<option value="100" selected>100</option>
									<option value="1000">1000</option>
									<option value="10000">10000</option>
								</select>
								<img style="visibility:hidden" id="tooltip_mc" class="helptooltip" title="The number of Monte Carlo simulations defines precision of the p-values. 100 MC can calcupate p-values up to 0.01, 1000 MC - up to 0.001, etc. Note Monte carlo simulations are highly computationally intensive and time consumint - use with consideration." style="position: relative;top: 10px" width=25 height=25 src="static/images/help-icon.png" alt="help"/>

								<label style="margin-left:10px;visibility:hidden">Strand selection: </label>
								<select name="strand" style="visibility:hidden">
									<option value="both" selected>Both</option>
									<option value="plus">Plus</option>
									<option value="minus">Minus</option>
								</select>
								<img class="helptooltip" title="Sets whether or not to use strand-specific regulatory datasets, if available. If a regulatory dataset does not have a strand, this setting is ignored" style="position: relative;top: 6px;visibility:hidden" width=25 height=25 src="static/images/help-icon.png" alt="help"/>								

							</div>
						</div>
						<div class="ui-state-highlight ui-corner-all">
							<p><span class="ui-icon ui-icon-info" style="float: left; margin-right: .3em;">
							</span>GenomeRunner works best in Chrome (recommended) or Firefox. Mac users, use Safari.</p>
						</div>

					</div>
				</div>
			</form>
			<p style="text-align: center;">
				<span style="font-family:arial,helvetica,sans-serif;">You are the&nbsp;<img alt="stats counter" border="0" src="http://www.easycounter.com/counter.php?mdozmorov" />&nbsp;visitor</span>
			</p>
		</div>
	</body>
