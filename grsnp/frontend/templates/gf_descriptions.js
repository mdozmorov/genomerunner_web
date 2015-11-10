$(document).ready(function() {	
	$.post("/get_gf_descriptions?db_version=${db_version}&organism=${organism}", function(data){
		data = jQuery.parseJSON(data);
		if (data.length != 0) {
			// insert blanks for GFs that do not have description
			for (var i = 0; i < data.length; i++) {
				if (data[i].length == 1){
					 data[i] = [data[i].toString(),"NA"];					
				}
			};
			// Create data table
			$('#descriptions').dataTable({
				"sDom": 'T<"clear">lfrtip',
				"bProcessing": true,
				"aaData": data,
				"aoColumns": [{ "sTitle": "Track Name" },
							{ "sTitle": "Description" }],
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
			$("#descriptions").html("<h3>Table could not be created.  Does the data exist?</h3>") 
		}
	});
});