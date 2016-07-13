shinyUI(fluidPage(
  ## for different themes install shinythemes
  #theme = shinytheme("spacelab"),
  #theme = "united.min.css",
  #theme = "style.css",
  # tags$head(
  #   tags$link(rel = "stylesheet", type = "text/css", href = "united.min.css")
  # ),

  div(titlePanel("Analysis of Sensory and Consumer data within a mixed effects model framework"), style = "color:#191970"),
  helpText("This application is a user-friendly interface for
      the R-package SensMixed"),
  fluidRow(
    column(4,
           uiOutput("antypeUI"),
#          bsTooltip("analysis", "title", placement = "bottom", trigger = "hover"),
          uiOutput("AttrUI"),
          #submitButton("Run Analysis")
         # actionButton("goButton", "Run Analysis"),
        bsButton("goButton", label = "Run Analysis", type = "action", 
                 style = "primary")
    ),
    column(8,
      uiOutput("theTabset")
  )
  )
  ))
