# NL-Logestic with bound V 0.4.5
# Payam Mokhtarian

##---------------------- Loading packages
library(shiny)

##---------------------- User interface
shinyUI(fluidPage(
  
  # Title
  titlePanel("Boundery decison in non-linear logistic regression"),
  
  # Sidebar controls
  sidebarLayout(
    sidebarPanel(
      selectInput("pattern", "Fitting Pattern:",
                  c("Convex" = "Convex",
                    "Close" = "Close")),
      sliderInput("degree",
                  "Degree Polynomial:",
                  min = 1,
                  max = 20,
                  value = 1),
      sliderInput("lambda",
                  "Lambda:",
                  min = 1,
                  max = 10,
                  value = 1),
      selectInput("opt", "Optimization Method:",
                  c("BFGS Quasi-Newton" = "BFGS",
                    "Nelder-Mead" = "Nelder-Mead",
                    "Conjugate Gradient" = "CG"))
    ),
    
    mainPanel(h4("Effective Genotype and Soil Acidity modelling"),
      plotOutput("da.plot")
    )
  )
))
