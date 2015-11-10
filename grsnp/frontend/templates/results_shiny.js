var refresh_progress = ""
$(document).ready(function() {	
	var run_id = "${run_id}";
});

	// generates both of the heatmap graphs
function add_heatmaps(){


	/*
	$("#heatmap").heatmap({
		ajaxuri: "/get_cluster?run_id=${run_id}",
		dwnl_link_id: "heatmap_download"		
	});
	
	$("#heatmap_cor").heatmap({
		ajaxuri: "/get_pcc?run_id=${run_id}",
		dwnl_link_id: "heatmap_cor_download"		
	});

	return false;
	*/	
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
