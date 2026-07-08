# preparations; required libraries
library(shiny)
library(shinyWidgets)

library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)

# the post url
post <- "https://github.com/CAPTURE-TB-HIV/tb-lmic-patientcosts-review"

# user interface elements and layout
ui <- fluidPage(
    setBackgroundColor("#BE91BE"),
    
    # titlePanel("TB PaCC"),
    
    titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
    
    # headerPanel(
    # 
    #     HTML('<p><img src="logo.svg"/></p>')
    # ),
    
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
        mainPanel(plotOutput(outputId = "bars"))
    )
)

# server-side computations
server <- function(input, output) {
    
    # the bar plot
    output$bars <- renderPlot({
        
        # most of this is for ggplot2; note the input$x syntax
        flips <- tibble(flips = rbinom(input$n, 1, input$prob)) %>% 
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
                          sprintf("%.2f", input$prob)))
    })
}

# run it all
shinyApp(ui = ui, server = server)

