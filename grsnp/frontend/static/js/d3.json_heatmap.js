  // jQuery Plugin Boilerplate
// A boilerplate for jumpstarting jQuery plugins development
// version 1.1, May 14th, 2011
// by Stefan Gabos

(function($) {

    $.heatmap = function(element, options) {

        var defaults = {
            ajaxuri: "/vis/get_data",
            dwnl_link_id: "",
            onFoo: function() {},
            cur_index: 0
        }

        var plugin = this;

        plugin.settings = {}

        var $element = $(element),
             element = element;

        plugin.init = function() {
            plugin.settings = $.extend({}, defaults, options);  
            $.post(plugin.settings.ajaxuri, function(res){

                data = jQuery.parseJSON(res);
                // check if the server returned an error
                if (data[plugin.settings.cur_index].matrix.substring(0,7) == "\"ERROR:" || 
                    data[plugin.settings.cur_index].matrix.substring(0,6) == "\"INFO:" ) {  
                  $(element).html("<h3>" + data[plugin.settings.cur_index].matrix.replace("\"","") + "</h3>"); 
                  return;
                }  
                 // Insert the row_name column row name
                $.each(data, function(i,k){
                   var lines_matrix = k.matrix.split("\n");
                   var header = "row_name\t" + lines_matrix[0]
                   k.matrix = header + "\n" + lines_matrix.slice(1,lines_matrix.length).join("\n")
                })
                // Parse the matrix data
                $.each(data,function(i,k){
                  k.matrix = d3.tsv.parse(k.matrix);
                });
                render_heatmap(element,data,0);  // Create heatmap

                // Create download link
                dwnl_link_id = plugin.settings.dwnl_link_id;
                if($('#heatmap').children('svg').length > 0 && $('#heatmap_download').length == 0){
                  $("#wellDownloads").append("<div style='margin:10px'><a style='font-size: 19px;margin-top:10px;' id='heatmap_download' >Download Clustered  Enrichment Matrix as SVG File (Right click - 'Save As' to save) </a><br></div>");
                }
                if($("#heatmap_cor").children('svg').length > 0 && $('#heatmap_cor_download').length == 0){
                  $("#wellDownloads").append("<div style='margin:10px'><a style='font-size: 19px' id='heatmap_cor_download' >Download Pearson's Correlation   Matrix as SVG File (Right click - 'Save As' to save) </a></div>");
                }
                if (dwnl_link_id != ""){
                  d3.selectAll("#" + dwnl_link_id)
                    .attr("href", "data:image/svg+xml;charset=utf-8;base64," + 
                      btoa(unescape(encodeURIComponent(
                        d3.selectAll("#"+element.id).selectAll("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg")
                       .node().parentNode.innerHTML)
                        )
                      )
                  );  
                }
            });

        }

        plugin.foo_public_method = function() {
            // code goes here
        }


        plugin.init();


        var sign = function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
         // ### helper functions
        var getMatrixValues = function(genes, colNames) {
              return genes.map(function(row) {
                return colNames.map(function(condition) {
                  return +row[condition];
                });
              });
            };

        var getColNames = function(matrix) {
              cols = [];
              for (k in matrix[0]){
                if (k != "row_name") { cols.push(k); }
              }
              return cols;
            };

        var isNumber = function(n) {
              return !isNaN(parseFloat(n)) && isFinite(n);
            };


         color_log10 = function(val) {
                if (Math.abs(val) <= 1) return 0
                var sign = 1
                if (val < 0) {sign = -1;}
                return sign * Math.log(Math.abs(val)) / Math.LN10;
              }

         // ### BACKBONE CODE ### 
         Heatmap = Backbone.View.extend({

            initialize: function() {      
              return this.render();
            },

            render: function() {   
              var cell_size, columns, colNames, conditionNamesMargin, extent, rowValues, rowNames, geneNamesMargin, getRow, heatmap, heatmapColor, height, margin, rows, textScaleFactor, width, x, y, legend;
              rowValues = this.model.get("rowValues");
              colNames = this.model.get("colNames");
              rowNames = this.model.get("rowNames");
              extent = this.model.get("extent");
              matrices = this.options.matrices;
              matrix_index = this.options.matrix_index;
              colNames = this.options.colNames;
              use_log = matrices[matrix_index].log

              // Sets the color range
              var num_range, col_range,legend;
              c_max_log = sign(color_range.max) * 1/(Math.pow(10,Math.pow(10,Math.abs(color_range.max))));
              c_min_log = sign(color_range.min) * 1/(Math.pow(10,Math.pow(10,Math.abs(color_range.min))));
              c_max = color_range.max
              c_min = color_range.min
              // case of complete over representation
              if (c_min >=0 && c_max >=0){
                num_range = [color_range.min,color_range.max];
                col_range = ["white","red"];
                diff= c_max - c_min; 
                legend = [c_min,diff*.25+c_min,diff*.50+c_min,diff*.75+c_min,c_max];
                
              } 
              // case of complete underrepresentation
              else if (c_min<=0 && c_max<=0){
                num_range = [color_range.min,color_range.max];
                col_range = ["green","white"];
                diff= c_min + c_max;   
                console.log("diff"+diff);
                legend = [c_min,diff*.75+c_max,diff*.50+c_max,diff*.25+c_max,c_max];
              }
              else if (c_min<=0 && c_max >= 0){
                num_range = [color_range.min,0,c_max];
                col_range = ["green","white","red"];
                diff= c_max - color_range.min;       
                legend = [c_min,c_min*.50,0,c_max*.5,c_max];
              }
              // if using log transform, then we want the labels to be 
              if (use_log) {
                legend_text = legend.map( function (d,i){ // the value used for the text of the legend
                  
                  return sign(d) * 1/(Math.pow(10,Math.pow(10,Math.abs(d)))) ;
                });
                legend_cell = legend;  // The value used for coloring the legend cells             
              }
              else{
                  legend_cell = legend;
                  legend_text = legend;
              }
              heatmapColor = d3.scale.linear().domain(num_range).range(col_range);
              

              textScaleFactor = 9;
              conditionNamesMargin = d3.max(colNames.map(function(conditionName) {
                return conditionName.length;
              }));
              geneNamesMargin = d3.max(rowNames.map(function(geneName) {
                return geneName.length;
              }));
              // If file names are too long, cap them to 20 characters
              if (conditionNamesMargin > 20) {
                conditionNamesMargin = 20;
              }
              if (geneNamesMargin > 30) {
                geneNamesMargin = 30;
              }


              // The buffer zone
              margin = {
                top: conditionNamesMargin * textScaleFactor,
                right: 150,
                bottom: 150,
                left: geneNamesMargin * textScaleFactor
              };
              cell_size = 40;
              width = cell_size * rowValues[0].length;
              height = cell_size * rowNames.length;
              $(this.el).html("") // This clears out old heatmap, if it already exists. Used to fix 'multiple heatmap bug'
              console.log("Creating heatmap at " + this.el.id)
              heatmap = d3.select(this.el).append("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg").attr("width", width + margin.right + margin.left)
                          .attr("height", height + margin.top + margin.bottom)
                          .append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
              x = d3.scale.ordinal().domain(d3.range(rowValues[0].length)).rangeBands([0, width]);
              y = d3.scale.ordinal().domain(d3.range(rowNames.length)).rangeBands([0, height]);
              columns = heatmap.selectAll(".column").data(colNames).enter().append("g").attr("class", "column").attr("transform", function(d, i) {
                return "translate(" + x(i) + ")rotate(-90)";
              });
              columns.append("text").attr("x", 6).attr("y", x.rangeBand() / 2).attr("dy", "-.5em").attr("dx", ".5em").attr("text-anchor", "start").attr("transform", "rotate(45)").text(function(d, i) {
                return colNames[i];
              });
             
             
              cur_row = -1;
              getRow = function(row,row_index) {
                var divtooltip = d3.select("body").append("div")   
                      .attr("class", "tooltip")               
                      .style("opacity", 0);

                    var cell;
                    cur_row += 1;
                    return cell = d3.select(this).selectAll(".cell").data(row).enter().append("rect").attr("class", "cell").attr("x", function(d, i) {
                      return x(i);
                    }).attr("width", x.rangeBand()).attr("height", x.rangeBand()).attr("row",matrix_data[cur_row].row_name)
                    .text(function(d,i) {
                      this.setAttribute("col",colNames[i]);
                      return d;
                    }).style("fill", function(d) {
                      if (matrices[matrix_index].log) return heatmapColor(color_log10(d));
                      else return heatmapColor(d);
                      
                    // ### Tool tips ###
                    }).attr("value",function(d){
                      return d
                    })
                    // Store the tooltip text
                    .attr("tt", function(d,iCol) {
                          alpha = matrices[matrix_index].alpha;
                          tip2 = ""   
                        if (alpha != "") {

                          d  = parseFloat(d);
                           if (d<0) { d = -1*(1/Math.pow(10,Math.abs(d))); }
                           else { d = 1/Math.pow(10,Math.abs(d)); }
                           
                           if (Math.abs(d) < alpha) { 
                              if (d<0) tip2 = "Underrepresented<br>";
                              else tip2 = "Overrepresented<br>";
                           }
                           else{
                             tip2 = "Not Significant<br>";
                           }
                           // Cycles through matrices       
                           $.each(matrices, function(i,cur_mat){
                                if (cur_mat.log){
                                   tip2 += cur_mat.name +":"
                                   tip2 += (1/Math.pow(10,Math.abs(cur_mat.matrix[row_index][colNames[iCol]]))).toExponential(2) + "</br>"
                                }                                  
                                else {
                                  tip2 += cur_mat.name + ": " + parseFloat(cur_mat.matrix[row_index][colNames[iCol]]).toExponential(2) + "</br>";
                                }
                           })              
                         }
                         else {
                           $.each(matrices, function(i,cur_mat){
                               tip2 += cur_mat.name + ": " + parseFloat(cur_mat.matrix[row_index][colNames[iCol]]).toFixed(2) + "</br>";
                           }) 
                         }
                         return  tip2;
                       }
                    )
                    .on("mouseover", function(d,i) {             
                
                        divtooltip.transition()        
                            .duration(50)      
                            .style("opacity", .9);

                       
                        divtooltip .html("<p align= 'left' style=\"color:#C1C1C1; margin-top: 4px; font-size: 16px;\">" + this.getAttribute("tt") + 
                                        "Row: " + this.getAttribute("row")  + "<br>" +
                                        "Column: " + this.getAttribute("col") + "</p>")  

                            .style("left", (d3.event.pageX+30) + "px")     
                            .style("top", (d3.event.pageY - 28) + "px");    
                        })                  
                      .on("mouseout", function(d) {
                        divtooltip.transition()        
                            .duration(50)      
                            .style("opacity", 0)});          
                    
                    };

                    var divtooltip = d3.select("body").append("div")   
                                      .attr("class", "tooltip")               
                                      .style("opacity", 0);   

                    rows = heatmap.selectAll(".row").data(rowValues).enter().append("g").attr("class", "row").attr("name", function(d, i) {
                      return "gene_" + i;
                    }).attr("transform", function(d, i) {
                      return "translate(0," + y(i) + ")";
                    }).each(getRow);
                    rows.append("text").attr("x", -6).attr("y", x.rangeBand() / 2).attr("dy", ".32em").attr("text-anchor", "end")
                        .text(function(d, i) {
                          return rowNames[i];
                    }); 

                    // ### The legend ###
                    var legend_y = parseInt($(this.el).find("rect").last().position().top) +  parseInt($(this.el).find("rect").last().attr("height")); // get the lowest point of the heatmap
                    var legend_c_size = 50, legend_margin = 50;
                    // Draw the boxes
                    heatmap.selectAll("rect.l_cell").data(legend_cell).enter().append("svg:rect").attr("class", "l_cell").attr("x", function(d, i) {
                           return i*legend_c_size;
                      }).attr("y", y(rowNames.length - 1) + cell_size + legend_margin)
                        .attr("width",legend_c_size).attr("height",legend_c_size).style("fill", function(d) {
                        return heatmapColor(d);
                      }).attr("value",function(d){
                        return d;
                      });

                    // Creates the value labels 
                     heatmap.selectAll(".legend").data(legend_text).enter().append("text").attr("class","legend").attr("x", function(d, i) {
                                  return i * legend_c_size;
                              })
                              .attr("y", y(rowNames.length-1) + legend_c_size + legend_margin + legend_c_size)
                              .attr("font-size", "13px")
                              .attr("text-anchor", "start")
                              .attr("transform", function(d,i){
                                return "rotate(45,"+String(i * legend_c_size) + "," + String( y(rowNames.length-1) + legend_c_size + legend_margin + legend_c_size)+")";
                              })
                              .text(function(d){
                                if(use_log){
                                  return d.toExponential(2);
                                }
                                else{
                                  return d.toPrecision(3);
                                }
                      });

              } // End Render

          }); // End Heatmap





        function render_heatmap(target,matrices,matrix_index) {             
              matrix_data  = matrices[matrix_index].matrix;


              matrix_range = function(arr,use_log) {

                arr = getMatrixValues(arr,getColNames(matrix_data));
                var max = -400;
                var min = 400;

                for (var i =0; i<arr.length;i++){
                  for (var o in arr[i]){
                    var val = arr[i][o];
                    if (isNumber(val)){
                      if(use_log){
                        val = color_log10(parseFloat(val));
                      } else { val = parseFloat(val);}
                      if (val > max){ max = val;}
                      if (val < min) { min = val; }
                    }
                  }
                }
                  return {"min": min, "max": max };
               }

             

               var geneExpressionModel, genes, heatmap, legend_c_size = 50, legend_y = $(target).attr("height") -100;

                color_range = matrix_range(matrix_data,matrices[matrix_index].log);
                create_heatmap(target,matrices,plugin.settings.cur_index);               
              

        }

          
          create_heatmap = function(target,matrices,matrix_index){
            matrix_data = matrices[matrix_index].matrix;
            // Modified to work with GR data
              matrixModel = new Backbone.Model;
              matrixModel.set({
                colNames: getColNames(matrix_data)
              });
              matrixModel.set({
                rowNames: matrix_data.map(function(row) {
                  return row.row_name;
                })
              });
              matrixModel.set({
                rowValues: getMatrixValues(matrix_data, matrixModel.get("colNames"))
              });
              matrixModel.set({
                extent: d3.extent($.map(matrixModel.get("rowValues"), function(item) {
                  return item;
                }))
              });
              return heatmap = new Heatmap({
                el: target,
                model: matrixModel,
                matrices: matrices,
                matrix_index: matrix_index,
                colNames:  matrixModel.get("colNames")
              });
            }
       




    }

    $.fn.heatmap = function(options) {

        return this.each(function() {
            if (undefined == $(this).data('heatmap')) {
                var plugin = new $.heatmap(this, options);               
                $(this).data('heatmap', plugin);
            }
        });

    }

})(jQuery);



 
 