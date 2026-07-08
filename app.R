# preparations; required libraries
library(shiny)
library(shinyWidgets)

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(bslib)

# the post url
post <- "https://github.com/CAPTURE-TB-HIV/tb-lmic-patientcosts-review"

# user interface elements and layout
ui <- page_navbar( id = "inTabset",

    nav_panel(title = "About", 
              
              headerPanel(
                  
                  # HTML('<p><img src="logoZoom.png" height="600" float: right /></p>')
                  HTML('<center> <img src="logoZoom.png" height="500"  />')
                  
              ),

              setBackgroundColor("#BE91BE"),
              
              p("The TB Patient Cost Calculator (TB PaCC)  estimates country-level TB patient costs for 129 low- and middle income countries. We chose the name TB PaCC to represent the significant planning and resources required among persons undergoing TB treatment. Furthermore, the economic burden of TB treatment can place financial strain on households, like carrying a heavy pack throughout a journey."),
              p("The underlying model is informed from data collected from 104 studies of TB patient costs, representing more than 55,000 persons treated for TB. The model uses routinely reported country- and year- matched data to estimate the cost of TB, stratified by three cost types (direct medical, direct non-medical, and indirect) and across the phases of the TB care cascade (pre-treatment, intensive treatment, and continuation treatment). All estimates represent TB patient costs in 2024 in 2024 USD."),
              p("Please note: TB patient cost data varies substantially across settings and study design. As such, our estimates should be used as a guide point for analyses."), 
              
              actionButton('jumpToPaCC', 'Get started!')),
    
    nav_panel(title = "PaCC" , 
              setBackgroundColor("#BE91BE"),
              
              titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
              
              
              sidebarLayout(
                  sidebarPanel(
                      
                      sliderInput(inputId = "n", label = "Pick a country:",
                                  min = 10, max = 1000, value = 500),
                      sliderInput(inputId = "prob", label = "Pick a study year:",
                                  min = 2009, max = 2024, value = 2024),
                      
                      tags$div(tags$p(HTML("<br><br><br><br>
                        Predictive model code available ")),
                        tags$a(href=post, "in this repository."))
                      
                  ),
                  mainPanel(
                      plotOutput(outputId = "bars"))
                      # tableOutput(outputId = "table")
                  )
              
              )
)

# server-side computations
server <- function(input, output, session) {
    
    observeEvent(input$jumpToPaCC, {
        updateTabsetPanel(session, "inTabset",
                          selected = "PaCC")
    })
    
    # output$table <- renderTable({     # most of this is for ggplot2; note the input$x syntax
    #     flips <- tibble(flips = rbinom(input$n, 1, input$prob)) %>% 
    #         mutate(flips = if_else(flips == 1, "Heads", "Tails"))  
    #     
    #     flips, striped = TRUE})
    
    
    # the bar plot
    output$bars <- renderPlot({

        # most of this is for ggplot2; note the input$x syntax
        flips <- tibble(flips = rbinom(input$n, 1, input$prob/10000)) %>%
            mutate(flips = if_else(flips == 1, "Heads", "Tails"))

        flips %>%
            count(flips) %>%
            ggplot(aes(flips, n, fill = flips)) +
            geom_col() +
            geom_label(aes(flips, n, label = n), size = 5) +
            theme(legend.position = "none",
                  axis.text = element_text(size = 15)) +
            labs(x = "", y = "") +
            ggtitle(str_c("Results of ", input$n,
                          " flips with Heads probability ",
                          sprintf("%.2f", input$prob/10000)))
    })
}

# run it all
shinyApp(ui = ui, server = server)

