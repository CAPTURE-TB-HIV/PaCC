# preparations; required libraries
library(shiny)
library(shinyWidgets)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(bslib)

#### Repository URL 
post <- "https://github.com/CAPTURE-TB-HIV/tb-lmic-patientcosts-review"


##############################################################################|
##### Read in the required data ##
##### If need to reduce bloat, this could be optimized to only read in 
##### data for selected country -- however could make switching slower for 
##### user. 
##############################################################################|
##### External covariate data 
covariates <- readRDS("utilities/data/covariates.rds")
##### Model parameter estimates 
params <- readRDS("utilities/data/parameters.rds")[c(1:54, 57),2] %>% as.numeric()
names(params) <- readRDS("utilities/data/parameters.rds")[c(1:54, 57),1]

##### Model estimates of country random effects
cntryREs <- readRDS("utilities/data/cntryRndEff.rds")

##### Available countries
countries <- covariates$study.country


##############################################################################|
##### Function to create a bar plot of covariate inputs for selected country ##
##############################################################################|

plot_inputs <- function(selectedCntry){
    covariates_filt <- covariates %>% filter(study.country == selectedCntry) %>% 
                       mutate(scale.median.income = median.income /100,
                              scale.median.income.trunc = median.income.trunc /100)
    
    
    ggplot(data = covariates_filt %>% 
                  pivot_longer(cols= c(scale.median.income.trunc, 
                               uhc.index.trunc,
                               oop.perc.che.trunc),names_to = "Covariate", values_to = "Value")) + 
        geom_col(aes(x=Covariate, y=Value), fill = "#A8DCAB", color="black") + theme_minimal(base_size = 14) + 
        ggtitle(paste("Values of covariates used to predict patient costs in", selectedCntry)) + 
        coord_flip() +
        expand_limits(y=c(0,120))
    
}


##############################################################################|
##### Function to calculate costs for selected country and return dataframe  ##
##############################################################################|

calculateCosts <- function(selectedCntry, 
                           selectedYr = 2024) {
    
    covs <- covariates %>% 
                       filter(study.country == selectedCntry) %>% 
                       mutate(scaled.med.inc = (log.median.income.trunc - 7.667872) / 0.8402527, 
                              scaled.oop.che = (log.median.income.trunc - 39.6697) / 16.0616, 
                              scaled.uhc =     (log.median.income.trunc - 61.46697) / 15.68885)
    
    cntryRE <- cntryREs %>% filter(Country == selectedCntry) %>% select(Value) %>% 
               unlist() %>% as.numeric()
    
    allCosts = data.frame(Type = c("Pre-treatment\ndirect medical", 
                                   "Pre-treatment\ndirect non-medical",
                                   "Pre-treatment\nindirect",
                                   
                                   "Intensive treatment\ndirect medical", 
                                   "Intensive treatment\ndirect non-medical",
                                   "Intensive treatment\nindirect", 
                                   
                                   "Continuation treatment\ndirect medical", 
                                   "Continuation treatment\ndirect non-medical",
                                   "Continuation treatment\nindirect"), 
                          Estimate = c(
                                  exp(params["beta_dme_prediag"] + cntryRE + 
                                      params["beta_med_income_dme_prediag"] * covs$scaled.med.inc +
                                      params["beta_oop_che_dme_prediag"] * covs$scaled.oop.che +
                                      params["beta_uhc_dme_prediag"] * covs$scaled.uhc) *
                                  exp(params["beta_year_dme_prediag"] * 1.991117) *       
                                  (1-params["prob_cost_zero"]), 
                                  
                                  ### DNME PREDIAGNOSIS 
                                  exp(params["beta_dnme_prediag"] + cntryRE + 
                                      params["beta_med_income_dnme_prediag"] * covs$scaled.med.inc +
                                      params["beta_oop_che_dnme_prediag"] * covs$scaled.oop.che +
                                      params["beta_uhc_dnme_prediag"] * covs$scaled.uhc) *
                                  exp(params["beta_year_dnme_prediag"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### INDIRECT PREDIAGNOSIS 
                                  exp(params["beta_ie_prediag"] + cntryRE + 
                                          params["beta_med_income_ie_prediag"] * covs$scaled.med.inc +
                                          params["beta_oop_che_ie_prediag"] * covs$scaled.oop.che +
                                          params["beta_uhc_ie_prediag"] * covs$scaled.uhc) *
                                      exp(params["beta_year_ie_prediag"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### DME NTENSIVE TREATMENT  
                                  exp(params["beta_dme_intsvtx"] + cntryRE + 
                                      params["beta_med_income_dme_intsvtx"] * covs$scaled.med.inc +
                                      params["beta_oop_che_dme_intsvtx"] * covs$scaled.oop.che +
                                      params["beta_uhc_dme_intsvtx"] * covs$scaled.uhc) *
                                      exp(params["beta_year_dme_intsvtx"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### DNME NTENSIVE TREATMENT  
                                  exp(params["beta_dnme_intsvtx"] + cntryRE + 
                                          params["beta_med_income_dnme_intsvtx"] * covs$scaled.med.inc +
                                          params["beta_oop_che_dnme_intsvtx"] * covs$scaled.oop.che +
                                          params["beta_uhc_dnme_intsvtx"] * covs$scaled.uhc) *
                                      exp(params["beta_year_dnme_intsvtx"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### INDIRECT INTENSIVE TREATMENT 
                                  exp(params["beta_ie_intsvtx"] + cntryRE + 
                                          params["beta_med_income_ie_intsvtx"] * covs$scaled.med.inc +
                                          params["beta_oop_che_ie_intsvtx"] * covs$scaled.oop.che +
                                          params["beta_uhc_ie_intsvtx"] * covs$scaled.uhc) *
                                      exp(params["beta_year_ie_intsvtx"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### DME CONTINUATION TREATMENT  
                                  exp(params["beta_dme_contntx"] + cntryRE + 
                                          params["beta_med_income_dme_contntx"] * covs$scaled.med.inc +
                                          params["beta_oop_che_dme_contntx"] * covs$scaled.oop.che +
                                          params["beta_uhc_dme_contntx"] * covs$scaled.uhc) *
                                      exp(params["beta_year_dme_contntx"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### DNME CONTINUATION TREATMENT  
                                  exp(params["beta_dnme_contntx"] + cntryRE + 
                                          params["beta_med_income_dnme_contntx"] * covs$scaled.med.inc +
                                          params["beta_oop_che_dnme_contntx"] * covs$scaled.oop.che +
                                          params["beta_uhc_dnme_contntx"] * covs$scaled.uhc) *
                                      exp(params["beta_year_dnme_contntx"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]), 
                                  
                                  ### INDIRECT CONTINUATION TREATMENT 
                                  exp(params["beta_ie_contntx"] + cntryRE + 
                                          params["beta_med_income_ie_contntx"] * covs$scaled.med.inc +
                                          params["beta_oop_che_ie_contntx"] * covs$scaled.oop.che +
                                          params["beta_uhc_ie_contntx"] * covs$scaled.uhc) *
                                      exp(params["beta_year_ie_contntx"] * 1.991117) *       
                                      (1-params["prob_cost_zero"]) 
                          ))
    
    return(allCosts)
}


##############################################################################|
##### Function to create a bar plot of calculated costs for selected country ##
##############################################################################|
plot_costs <- function(selectedCntry){
    
              costs_df <- calculateCosts(selectedCntry = selectedCntry)
    
              ggplot(costs_df) + 
                geom_col(aes(x=Type, y=Estimate), fill = "#A8DCAB", color="black") + 
                ggtitle(paste("Per-episode household cost of TB in", selectedCntry, "(2024 USD)")) + 
                coord_flip() + theme_minimal(base_size = 14) + 
                expand_limits(y=c(0,max(costs_df$Estimate) * 1.25)) +
              xlab("") + ylab("Costs (2024 USD)")
    
}

###############################################################################
###############################################################################
###############################################################################

###############################################################################
##### user interface elements and layout
###############################################################################

ui <- page_navbar( id = "inTabset",

    nav_panel(title = "About", 
              
              headerPanel(
                  
                  HTML('<center> <img src="logoZoom.png" height="500"  />')
                  
              ),

              setBackgroundColor("#BE91BE"),
              
              p("The TB Patient Cost Calculator (TB PaCC)  estimates country-level TB patient costs for 129 low- and middle income countries. We chose the name TB PaCC to represent the significant planning and resources required among persons undergoing TB treatment. Furthermore, the economic burden of TB treatment can place financial strain on households, like carrying a heavy pack throughout a journey."),
              p("The underlying model is informed from data collected from 104 studies of TB patient costs, representing more than 55,000 persons treated for TB. The model uses routinely reported country- and year- matched data to estimates the cost of TB, stratified by three cost Types (direct medical, direct non-medical, and indirect) and across the phases of the TB care cascade (pre-treatment, intensive treatment, and continuation treatment). All estimates represent TB patient costs in 2024 in 2024 USD."),
              p("Please note: TB patient cost data varies substantially across settings and study design. As such, our estimates should be used as a guide point for analyses."), 
              
              actionButton('jumpToPaCC', 'Get started!', icon("paper-plane"),
                           style="color: #000000; background-color: #FFFFFF")),
    
    nav_panel(title = "PaCC" , 
              setBackgroundColor("#BE91BE"),
              
              titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
              
              ##### Sidebar which includes the country dropdown selection #####
              sidebarLayout(
                  sidebarPanel(
                      
                      selectInput(
                          "country",
                          label = "Select a country",
                          choices = countries
                      ),
                      
                      tags$div(tags$p(HTML("<br><br>
                        Predictive model code available "), tags$a(href=post, "in this repository.")))
                      
                  ),
                  ##### Main panel of results on right of screen
                  mainPanel(
                      ##### Create an inset of tabs for the user to click through #####
                      navset_card_underline(
                          title = "",
                          ##### Plot of costs #####
                          nav_panel("Cost plots", 
                          plotOutput(outputId = "costbar")),
                          ##### Table of costs #####
                          nav_panel("Cost table",
                                    tableOutput(outputId = "costTable")),
                          ##### Plot of covariates #####
                          nav_panel("Covariate plots", 
                                    plotOutput(outputId = "covariates"))
                  )
              
              )
)
)
)

###############################################################################
##### server side 
###############################################################################

server <- function(input, output, session) {
    
    #### Observe the country selected 
    observeEvent(input$jumpToPaCC, {
        updateTabsetPanel(session, "inTabset",
                          selected = "PaCC")
    })
    
    
        #### Plot country specific costs as bar plot
        output$costbar <- renderPlot({
            plot_costs(input$country)
        })
        
        #### Present country specific costs as table
        output$costTable <- renderTable({
            calculateCosts(input$country)
        })
        
        #### Plot country specific covariates as bar plot
        output$covariates <- renderPlot({
            plot_inputs(input$country)
        })
    
}

###############################################################################
##### Run everything 
###############################################################################

shinyApp(ui = ui, server = server)

