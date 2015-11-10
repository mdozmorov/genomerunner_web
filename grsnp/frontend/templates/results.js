var refresh_progress = ""
$(document).ready(function() {	
	$("#progressbar").progressbar({
						value: ${curprog},
						max:  ${progmax}
					});
	update_progress();
	refresh_progress = setInterval(update_progress, 10000);	
	$(".helptooltip").tooltip();	
});



function update_progress(){
	$.post('/get_progress?run_id=${run_id}',function(data){
		data = jQuery.parseJSON( data );
		// update progress bar
		$("#progressbar").progressbar({
						value: data['curprog'],
						max:  data['progmax']});
		if (data['status'] != "") $("#status").html(data['status']);
		if (data['status'] == "Analysis Completed" || data['status'].substring(0,18) == "Running Annotation"){ 
						
			$("#divDownload").html("<a class='btn btn-primary' style='margin-left: 230px;'  type='button' href='${zipfile}'>Download All Run Files</a>")
		}
		if(data['status'] == "Run crashed. See end of log for details."){
			$("#progressbar").progressbar({
						value: 1,
						max:  1});
			clearInterval(refresh_progress);
		} else if (data['status'] == "Analysis Completed"){
			window.location.replace("http://" + window.location.host + "/results_shiny?id=${run_id}");
		}
	});
	get_log();
	get_detailed();
}

function get_log(){
	$.post('/get_log?run_id=${run_id}',function(data){
		data = jQuery.parseJSON( data );
		$("#txtLog").html(data['log']);
	});
}

function get_detailed(){
	$.post('/get_detailed?run_id=${run_id}',function(data){
		data = jQuery.parseJSON( data );
		$("#txtDetailed").html(data['detailed']);
	});
}

// Sorting functions for the P-values
jQuery.extend( jQuery.fn.dataTableExt.oSort, {
    "scientific-pre": function ( a ) {
        return parseFloat(a);
    },
 
    "scientific-asc": function ( a, b ) {
        return ((a < b) ? -1 : ((a > b) ? 1 : 0));
    },
 
    "scientific-desc": function ( a, b ) {
        return ((a < b) ? 1 : ((a > b) ? -1 : 0));
    }
} );

// Sorting functions for annotation results
jQuery.extend( jQuery.fn.dataTableExt.oSort, {
    "annotation-pre": function ( a ) {
        return parseInt(a.split("|")[0]);
    },
 
    "annotation-asc": function ( a, b ) {
        return ((a < b) ? -1 : ((a > b) ? 1 : 0));
    },
 
    "annotation-desc": function ( a, b ) {
        return ((a < b) ? 1 : ((a > b) ? -1 : 0));
    }
} );

function get_annotation(foi_name){
	$.post("/get_annotation?run_id=${run_id}&foi_name="+foi_name, function(data){
		data = jQuery.parseJSON(data);
		if (data.length != 0) {
			// clear old data tables.
										
			$(".dt_annotations").each(function (){ 
				$(this).children().first().empty();
				tmp = $(this).attr("id").split("_")
				tmp = tmp.slice(2,tmp.length).join("_")
				 $(this).children().first().append("<table id='dt_" + tmp + "'></table>") 						

			});

			// Get column names
			var cols = []
			var cols = []
			var c = data.shift()
			for(var k in c){
				if (c[k] == "Region"){
					cols.push({"sTitle": c[k], "sType": "annotation"})
				}
				else {
					cols.push({"sTitle": c[k]})
				}
			}
			// Create data table
			$('#dt_'+foi_name).dataTable({
				"sDom": 'T<"clear">lfrtip',
				"bProcessing": true,
				"aaData": data,
				"aoColumns": cols,
				"asSorting": [],
				 "oTableTools": {
			            "aButtons": [
			                "copy",
			                "print",
			                {
			                    "sExtends":    "collection",
			                    "sButtonText": "Save",
			                    "aButtons":    [ "csv", "xls", "pdf" ]
			                }
			            ]
			        }
			});
		}
		else{
			$("#dt_"+foi_name).html("<h3>Table could not be created.  Does the data exist?</h3>") 
		}
	});
}

function get_enrichment(foi_name){
	$.post("/get_enrichment?run_id=${run_id}&foi_name="+foi_name, function(data){
		data = jQuery.parseJSON(data);
		if (data.length != 0) {
			// clear old data tables.												
			$(".dt_enrichments").each(function (){ 
				$(this).children().first().empty();
				tmp = $(this).attr("id").split("_")
				tmp = tmp.slice(2,tmp.length).join("_")
				 $(this).children().first().append("<table id='dter_" + tmp + "'></table>") 						

			});
			// Get column names
			var cols = []
			var c = data.shift()
			for(var k in c){
				if (c[k] == "P.value" || c[k] == 'P.adj'){
					cols.push({"sTitle": c[k], "sType": "scientific"})
				}
				else {
					cols.push({"sTitle": c[k]})
				}
			}

			// Create data table
			$('#dter_'+foi_name).dataTable({
				"sDom": 'T<"clear">lfrtip',
				"bProcessing": true,
				"aaData": data,
				"aoColumns": cols,
				 "oTableTools": {
			            "aButtons": [
			                "copy",
			                "print",
			                {
			                    "sExtends":    "collection",
			                    "sButtonText": "Save",
			                    "aButtons":    [ "csv", "xls", "pdf" ]
			                }
			            ]
			        }
			});
		}
		else{
			$("#dter_"+foi_name).html("<h3>Table could not be created.  Data does not exist?</h3>") 
		}
	});
}
