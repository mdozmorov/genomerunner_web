library(shiny)
library(d3heatmap)
library(DT)

shinyUI(
  tags$head(
  fluidPage(tags$style(HTML("
    .progress-striped .bar {
      background-color: #FF3300;
      background-image: -webkit-gradient(linear, 0 100%, 100% 0, color-stop(0.25, rgba(255, 255, 255, 0.6)), color-stop(0.25, transparent), color-stop(0.5, transparent), color-stop(0.5, rgba(255, 255, 255, 0.6)), color-stop(0.75, rgba(255, 255, 255, 0.6)), color-stop(0.75, transparent), to(transparent));
      background-image: -webkit-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
      background-image: -moz-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
      background-image: -o-linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
      background-image: linear-gradient(45deg, rgba(255, 255, 255, 0.6) 25%, transparent 25%, transparent 50%, rgba(255, 255, 255, 0.6) 50%, rgba(255, 255, 255, 0.6) 75%, transparent 75%, transparent);
      -webkit-background-size: 60px 60px;
         -moz-background-size: 60px 60px;
           -o-background-size: 60px 60px;
              background-size: 60px 60px;
    }
   .progress-text {
    background-color: #FFCC00;
   }
  "))),
    sidebarLayout(
      uiOutput("sidebar"),
      mainPanel(
        tags$head(tags$style("td {align:center;}")),
        uiOutput("mainpage")
      )
    )
))
