<body bgcolor="#A8D5FF">	
	<div class="topbar">
		<div class="topbar-inner">
			<div class="container-fluid">
				<a class="brand" href="./" style="font-family:Futura,'Helvetica Neue', Arial">
					<img width="200px" src="static/images/GRLogo.png" />
				</a>
				<ul class="pull-right">
					<li><a href="./overview">Overview</a></li>
					<li><a href="./news">News</a></li>
					<li><a href="./demo">Demo</a></li>
					<!-- <li><a href="./cite">How to Cite</a></li>
					<li><a href="http://sourceforge.net/projects/genomerunner/">GenomeRunner on SourceForge</a></li>
					<li><a href="./roadmap">Roadmap</a></li> -->
					<li><a href="https://mdozmorov.github.io/grdocs/">Help</a></li>
<!-- <div>
					<li><a href="./google">Google Group</a></li>
					<li><img width="30px" src="static/new-icon.jpg" alt="New: GenomeRunnerSNP Google Groups" /></li>
</div> -->
				</ul>
			</div>
		</div>
	</div>
<div class="well">
	<div style="margin: 10px;margin-top:70px;">
	<div id="progressbar" ></div>
	<div class="ui-widget" style="">
		<div class="ui-state-highlight ui-corner-all">
			<p><span class="ui-icon ui-icon-info" style="float: left; margin-right: .3em;"></span><label id="status">Analysis job is in the queue.</label></p>			
		</div>
		<p><span class="ui-icon ui-icon-info" style="float: left; margin-right: .3em;"></span>Please bookmark this page to return to your results. The results will be deleted after three days.</p>
	</div>
	<!--<a href="./results_shiny?id=${run_id}" class="btn btn-info" role="button">Shiny Analysis</a>-->
</div>
</div>
<ul class="tabs" id="resultstab" data-tabs="tabs">
	<li ><a href="#detailedresults" data-toggle="tab" onclick="get_detailed()" >Detailed Results<img class="helptooltip" title="P-value calculations" style="position: relative;top: 6px;" width="25" height="25" src="static/images/help-icon.png" alt="help"></a></li>

	<li class="active"><a href="#debugresults" data-toggle="tab" onclick="get_log()">Analysis Log<img class="helptooltip" title="Log, warning, and error messages" style="position: relative;top: 6px;" width="25" height="25" src="static/images/help-icon.png" alt="help"></a></li>
	% if false:
 	<li><a href="#annotationresults" data-toggle="tab">Annotation Results<img class="helptooltip" title="Contains Annotation Results" style="position: relative;top: 6px;" width="25" height="25" src="static/images/help-icon.png" alt="help"></a></li> 
 	% endif
	% if false:
	<li><a href="#heatmapdownload" data-toggle="tab">Download Results<img class="helptooltip" title="Download tab-delimited results, and heatmaps" style="position: relative;top: 6px;" width="25" height="25" src="static/images/help-icon.png" alt="help"></a></li>
	% endif


</ul>

<div class="tab-content">	
	<div class="tab-pane" id="detailedresults">
		<div class="well">
			<h3>Detailed Results</h3>
			<textarea  id="txtDetailed" cols="150" rows="30"></textarea><br>			
		</div>
	</div>
	<div class="tab-pane active" id="debugresults">
		<div class="well">
			<h3>Analysis Log</h3>
			<textarea id="txtLog" cols="150" rows="30"></textarea><br>
		</div>
	</div>
	% if false: #run_annotation:
 	<div class="tab-pane" id="annotationresults">
			<ul class="tabs" id="annotationtab" data-tabs="tabs">
				% for i,item in enumerate(fois):
				  <li><a href="#ano_tab_${item}" data-toggle="tab" onclick="get_annotation('${item}')">${item}</a></li>
				% endfor
			</ul>
			<div class="tab-content">
				% for i,item in enumerate(fois):
				<div class="dt_annotations" ${"class='tab-pane active'" if i==0 else "class='tab-pane'"} id="ano_tab_${item}">
					<div class="well">
						<table id="dt_${item}"></table>
					</div>
				</div>
				% endfor
			</div>	
	</div>
	% endif
	% if false:
	<div class="tab-pane" id="heatmapdownload">
		<div id="wellDownloads" class="well">
			<div id="divDownload" style="margin: 20px">
				<a class="btn btn-primary disabled" style="margin-left: 230px;" type="button">Run files not ready for download</a>
			</div>						
		</div>
	</div>
	%endif
</div>

</body>