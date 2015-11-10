 var Heatmap, getConditionNames, getGeneExpressions, isNumber,cur_heatmap,cur_tooltip_matrix, tooltip_matrices = {},log_transform_color = false;
  sign = function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
 // ### helper functions
   getGeneExpressions = function(genes, conditionNames) {
      return genes.map(function(gene) {
        return conditionNames.map(function(condition) {

          return +gene[condition];
        });
      });
    };

    getConditionNames = function(genes) {
      return Object.keys(genes[0]).filter(function(columnName) {
        return !columnName.match(/cluster/) && isNumber(genes[1][columnName]);
      });
    };

    isNumber = function(n) {
      return !isNaN(parseFloat(n)) && isFinite(n);
    };
 // ###
 Heatmap = Backbone.View.extend({

    initialize: function() {      
     
      return this.render();
    },

    render: function() {   
      var cell_size, clusterColor, clusters, columns, conditionNames, conditionNamesMargin, extent, geneExpressions, geneNames, geneNamesMargin, getRow, heatmap, heatmapColor, height, margin, rows, textScaleFactor, width, x, y, legend;
      geneExpressions = this.model.get("geneExpressions");
      conditionNames = this.model.get("conditionNames");
      geneNames = this.model.get("geneNames");
      tooltip_matrices[cur_heatmap] = getGeneExpressions(cur_tooltip_matrix, conditionNames); // for tooltips
      extent = this.model.get("extent");
      clusters = this.model.get("clusters");
      clusterColor = this.model.get("clusterColor");

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
        legend_log = [c_min_log,diff*.25+c_min_log,diff*.50+c_min_log,diff*.75+c_min_log,c_max_log];
        
      } 
      // case of complete underrepresentation
      else if (c_min<=0 && c_max<=0){
        num_range = [color_range.min,color_range.max];
        col_range = ["green","white"];
        diff= c_min + c_max;   
        console.log("diff"+diff);
        legend = [c_min,diff*.75+c_max,diff*.50+c_max,diff*.25+c_max,c_max];
        legend_log = [c_min_log,diff*.75+c_max_log,diff*.50+c_max_log,diff*.25+c_max_log,c_max_log];
      }
      else if (c_min<=0 && c_max >= 0){
        num_range = [color_range.min,0,c_max];
        col_range = ["green","white","red"];
        diff= c_max - color_range.min;       
        legend = [c_min,c_min*.50,0,c_max*.5,c_max];
        legend_log = [c_min_log,c_min_log*.50,0,c_max_log*.5,c_max_log];
      }
      heatmapColor = d3.scale.linear().domain(num_range).range(col_range);
      

      textScaleFactor = 9;
      conditionNamesMargin = d3.max(conditionNames.map(function(conditionName) {
        return conditionName.length;
      }));
      geneNamesMargin = d3.max(geneNames.map(function(geneName) {
        return geneName.length;
      }));

      // The buffer zone
      margin = {
        top: 0, //conditionNamesMargin * textScaleFactor,
        right: 150,
        bottom: conditionNamesMargin * textScaleFactor +150,
        left: geneNamesMargin * textScaleFactor
      };
      cell_size = 40;
      width = cell_size * geneExpressions[0].length;
      height = cell_size * geneNames.length;
      $(this.el).html("") // This clears out old heatmap, if it already exists. Used to fix 'multiple heatmap bug'
      console.log("Creating heatmap at " + this.el.id)
      heatmap = d3.select(this.el).append("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg").attr("width", width + margin.right + margin.left).attr("height", height + margin.top + margin.bottom).attr("id", "heatmap").append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
      x = d3.scale.ordinal().domain(d3.range(geneExpressions[0].length)).rangeBands([0, width]);
      y = d3.scale.ordinal().domain(d3.range(geneNames.length)).rangeBands([0, height]);
      columns = heatmap.selectAll(".column").data(conditionNames).enter().append("g").attr("class", "column").attr("transform", function(d, i) {
        return "translate(" + x(i) + ")rotate(-90)";
      });
      columns.append("text").attr("x", 6).attr("y", x.rangeBand() / 2).attr("dy", "-.5em").attr("dx", ".5em").attr("text-anchor", "start").attr("transform", "rotate(45)").text(function(d, i) {
        return conditionNames[i];
      });


     // Create the legend
     // NOTE: only works for the #heatmap NOT the #heatmap_cor
      var root = "svg#"+this.el.getAttribute("id");
      var legend_y = height + margin.bottom -  100
      var legend_c_size = 50
      // Draw the boxes
      d3.select(root).selectAll("rect").data(legend).enter().append("svg:rect").attr("class", "l_cell").attr("x", function(d, i) {
             return i*legend_c_size+40;
        }).attr("y", legend_y).attr("width",legend_c_size).attr("height",legend_c_size).style("fill", function(d) {
          return heatmapColor(d);
        }).attr("value",function(d){
          return d;
        });

      // Creates the value labels 
       d3.select(root).selectAll(".legend").data(legend_log).enter().append("text").attr("class","legend").attr("x", function(d, i) {
                      console.log(d);
                    return (i * legend_c_size) +45;
                })
                .attr("y", legend_y + legend_c_size + 30)
                .attr("font-size", "13px")
                .attr("text-anchor", "start")
                .attr("transform", function(d,i){
                  return "rotate(45,"+String((i * legend_c_size) +45) + "," + String(legend_y + legend_c_size + 30)+")";
                })
                .text(function(d){
                    return d.toExponential(2); // toPrecision(3);
                });
      

       
      color_log10 = function(val) {
        if (Math.abs(val) <= 1) return 0
        var sign = 1
        if (val < 0) {sign = -1;}
        return sign * Math.log(Math.abs(val)) / Math.LN10;
      }
     
     
      cur_row = -1;
      getRow = function(row) {
        var divtooltip = d3.select("body").append("div")   
          .attr("class", "tooltip")               
          .style("opacity", 0);

        var cell;
        cur_row += 1;
        return cell = d3.select(this).selectAll(".cell").data(row).enter().append("rect").attr("class", "cell").attr("x", function(d, i) {
          return x(i);
        }).attr("width", x.rangeBand()).attr("height", x.rangeBand()).attr("row",matrix[cur_row].gene_name)
        .text(function(d,i) {
          this.setAttribute("tt",tooltip_matrices[cur_heatmap][cur_row][i]);
          this.setAttribute("col",conditionNames[i]);
          return d;
        }).style("fill", function(d) {
          if (log_transform_color) return heatmapColor(color_log10(d));
          else return heatmapColor(d);
          
      // Tool tips
        }).attr("value",function(d){
          return d
        })
        .on("mouseover", function(d,i) {             
    
            divtooltip.transition()        
                .duration(50)      
                .style("opacity", .9);

            if (log10_tt_value == true) {

               d  = parseFloat(this.getAttribute("value"));
               if (d<0) { d = -1*(1/Math.pow(10,Math.abs(d))); }
               else { d = 1/Math.pow(10,Math.abs(d)); }
               
               if (Math.abs(d) < .05) { 
                  if (d<0) tip1 = "Underrepresented<br>";
                  else tip1 = "Overrepresented<br>";
               }
               else{
                 tip1 = "Not Significant<br>";
               }             
               tip2 = "P-value: "  + d.toExponential(4) + "<br>";               
             }
             else {
                tip1 = parseFloat(this.getAttribute("value")).toPrecision(2);
                tip1 = "Pearsons: " + tip1 + "<br>"
                tip2 = parseFloat(this.getAttribute("tt")).toExponential(4);
                tip2 = "P-value: " + tip2 + "<br>"
             }
            divtooltip .html("<p style=\"color:#C1C1C1; margin-top: 4px; font-size: 16px;\">" + tip1 + tip2 + 
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

      rows = heatmap.selectAll(".row").data(geneExpressions).enter().append("g").attr("class", "row").attr("name", function(d, i) {
        return "gene_" + i;
      }).attr("transform", function(d, i) {
        return "translate(0," + y(i) + ")";
      }).each(getRow);
      return rows.append("text").attr("x", -6).attr("y", x.rangeBand() / 2).attr("dy", ".32em").attr("text-anchor", "end").text(function(d, i) {
        desc = matrix_data_gf_description[i]
        if (desc == "") this.setAttribute("tt", "No Description");
        else this.setAttribute("tt", desc);
        return geneNames[i];
      }).on("mouseover", function(d,i) {
            if (cur_heatmap == "heatmap") {              
              divtooltip.transition()        
                  .duration(50)      
                  .style("opacity", .9);
              
              divtooltip .html("<p style=\"color:#C1C1C1; margin-top: 4px; font-size: 16px;\">" + this.getAttribute("tt") + "</p>")  

                  .style("left", (d3.event.pageX+30) + "px")     
                  .style("top", (d3.event.pageY - 28) + "px");    
            }
          })                  
          .on("mouseout", function(d) {  
           if (cur_heatmap == "heatmap") {          
            divtooltip.transition()        
                .duration(50)      
                .style("opacity", 0)}
              }); 



          }


  });

function generate_heatmaps() {

      matrix_range = function(arr,log_trans) {
        var max = -400;
        var min = 400;

        for (var i =0; i<arr.length;i++){
          for (var o in matrix[i]){
            var val = matrix[i][o];
            if (isNumber(val)){
              if(log_trans){
                val = color_log10(parseFloat(val));
              } else { val = parseFloat(val);}
              if (val > max){ max = val;}
              if (val < min) { min = val; }
            }
          }
        }
          return {"min": min, "max": max };
       }

      color_log10 = function(val) {
        if (Math.abs(val) <= 1) return 0
        var sign = 1
        if (val < 0) {sign = -1;}
        return sign * Math.log(Math.abs(val)) / Math.LN10;
      }


      // The code that create the heatmaps
      $(document).ready(function() {
        matrix_data_gf_description = matrix_data_gf_description.split("\t");
        var geneExpressionModel, genes, heatmap
        matrix_cor_pvalues = d3.tsv.parse(matrix_cor_pvalues)
        log10_tt_value = false;
        matrix_cor = d3.tsv.parse(matrix_cor)
        matrix = matrix_cor;
        color_range = {"min": -1,"max": 1}; // Sets the color max to the largest value in the matrix
        console.log("range cor: "+ color_range.min + " : "+ color_range.max);

        cur_heatmap = "heatmap_cor";
        cur_tooltip_matrix = matrix_cor_pvalues;    
        create_heatmap("#heatmap_cor", matrix,matrix_cor);
        // NOTE: Draw the legend for the PCC heatmap. In the future, this should be done in the same place as  the legend for the "#heatmap"
        var root = "#heatmap_cor svg";
        var legend_y = $(root).attr("height") -100
        var legend_c_size = 50
        // Draw the boxes
        var PCC_color =  d3.scale.linear().domain([-1,0,1]).range(["green","white","red"]);
        d3.select(root).selectAll(".l_cell").data([-1,-.5,0,.5,1]).enter().append("svg:rect").attr("class", "l_cell").attr("x", function(d, i) {
               return i*legend_c_size+40;
          }).attr("y", legend_y).attr("width",legend_c_size).attr("height",legend_c_size).style("fill", function(d) {
            return PCC_color(d);
          }).attr("value",function(d){
            return d;
          });

        // Creates the value labels 
         d3.select(root).selectAll(".legend").data([-1,-.5,0,.5,1]).enter().append("text").attr("class","legend").attr("x", function(d, i) {
                        console.log(d);
                      return (i * legend_c_size) +45;
                  })
                  .attr("y", legend_y + legend_c_size + 30)
                  .attr("font-size", "13px")
                  .text(function(d){
                      return d.toPrecision(3);
                  });

        matrix_enrich_data = d3.tsv.parse(matrix_data);
        matrix = matrix_enrich_data
        log10_tt_value = true;
        cur_heatmap = "heatmap";
        cur_tooltip_matrix = matrix_enrich_data;
        color_range = matrix_range(matrix,true);
        log_transform_color = true;
        console.log("range: "+ color_range.min + " : "+ color_range.max);
        create_heatmap("#heatmap",matrix);
        // Creates links to download the svg file
        // Commented below is the code from 08/23/2013 commit 243f676. These changes broke SVG download.
       //  var svg = d3.selectAll("#heatmap_cor").selectAll("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg"); // or whatever you call it
       //  var serializer = new XMLSerializer();
       //  var str = serializer.serializeToString(svg);
       //  console.log(str);
       // $("heatmap").remove()
       //  d3.selectAll("#heatmap_download")
       //      .attr("href", "data:image/svg+xml;charset=utf-8;base64," + 
       //        btoa(unescape(encodeURIComponent(
       //            d3.selectAll("#heatmap svg").node().parentNode.innerHTML
       //          )
       //        )
       //      ));
        // Code before 243f676 commit, working download
        d3.selectAll("#heatmap_download")
            .attr("href", "data:image/svg+xml;charset=utf-8;base64," + 
              btoa(unescape(encodeURIComponent(
                d3.selectAll("#heatmap").selectAll("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg")
               .node().parentNode.innerHTML)
                )
              )
            );
         // Code unchanged through the commits
        d3.selectAll("#heatmap_cor_download")
            .attr("href", "data:image/svg+xml;charset=utf-8;base64," + 
              btoa(unescape(encodeURIComponent(
                d3.selectAll("#heatmap_cor").selectAll("svg").attr("version", "1.1").attr("xmlns", "http://www.w3.org/2000/svg")
               .node().parentNode.innerHTML)
                )
              )
        );
      });
 

 

}

  create_heatmap = function(target,matrix){
    // Modified to work with GR data
    // matrix_gfs, matrix_fois, matrix are given values in results.html

      if (matrix.length == 0){ 
          heatmap = d3.select("#heatmap").append("svg:heatmap").attr("width", 1195).attr("height", 500);
          heatmap.append("svg:rect")
          .attr("width", 1195)
          .attr("height", 500)
          .attr("fill","white")
          heatmap.append("svg:text")
          .attr("x", 500)
          .attr("y", 200)
          .style("font-size","20px")
          .attr("text-anchor", "middle")
          .text(matrix_data.replace("gene_name",""));
          return; 
      }
      geneExpressionModel = new Backbone.Model;
      geneExpressionModel.set({
        conditionNames: getConditionNames(matrix)
      });
      geneExpressionModel.set({
        geneNames: matrix.map(function(gene) {
          return gene.gene_name;
        })
      });
      geneExpressionModel.set({
        geneExpressions: getGeneExpressions(matrix, geneExpressionModel.get("conditionNames"))
      });
      geneExpressionModel.set({
        extent: d3.extent($.map(geneExpressionModel.get("geneExpressions"), function(item) {
          return item;
        }))
      });

      geneExpressionModel.set({
        clusters: matrix.map(function(gene) {
          return gene.cluster;
        })
      });
      geneExpressionModel.set({
        clusterColor: d3.scale.category20()
      });
      return heatmap = new Heatmap({
        el: target,
        model: geneExpressionModel
      });
    }




   

  // changes the tooltip matrix so that the correct values can be displayed by the tooltips for the 
  // currently active matrix
  change_active_heatmap = function(heatmap_name){
    cur_heatmap = heatmap_name;
    console.log(cur_heatmap);
    if (heatmap_name == "heatmap"){      
      log10_tt_value = true;
    }
    else if (heatmap_name == "heatmap_cor"){
      log10_tt_value = false;
    }
  }

 