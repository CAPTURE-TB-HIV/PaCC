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

##### Relevant studies 
studies <- readRDS("utilities/data/tbPatientCostStudies.rds")


##############################################################################|
##### Function to create a bar plot of covariate inputs for selected country ##
##############################################################################|

lookup_inputs <- function(selectedCntry){
    covariates_filt <- covariates %>% filter(study.country == selectedCntry) %>% 
        mutate(scale.median.income = median.income /100,
               scale.median.income.trunc = median.income.trunc /100) %>% 
        rename("Median individual\nannual income\n(in 100s of 2024 USD)" = scale.median.income.trunc, 
               "Universal healthcare\ncoverage index" = uhc.index.trunc, 
               "Out-of-pocket percentage\nof current health\nexpenditures" = oop.perc.che.trunc)
    
    return(covariates_filt)
}

plot_inputs <- function(selectedCntry){

    covariates_filt <- lookup_inputs(selectedCntry)
    
    ggplot(data = covariates_filt %>% 
                  pivot_longer(cols= c(
                       `Median individual\nannual income\n(in 100s of 2024 USD)`, 
                       `Universal healthcare\ncoverage index`,
                       `Out-of-pocket percentage\nof current health\nexpenditures`),names_to = "Covariate", values_to = "Value")) + 
        geom_col(aes(x=Covariate, y=Value), fill = "#A8DCAB", color="black") + theme_minimal(base_size = 14) + 
        ggtitle(paste("Values of covariates used to predict patient costs in", selectedCntry)) + 
        coord_flip() +
        expand_limits(y=c(0,120)) + xlab("") + 
        geom_text(aes(x=Covariate, y=Value, label=round(Value)), hjust = -.75)
    
}


##############################################################################|
##### Function to calculate costs for selected country and return dataframe  ##
##############################################################################|

calculateCosts <- function(selectedCntry, 
                           medIncInput, 
                           oopInput,
                           uhcInput,
                           selectedYr = 2024) {
    
    covs <- covariates %>% 
                       filter(study.country == selectedCntry) %>% 
                       mutate(scaled.med.inc = (log.median.income.trunc - 7.667872) / 0.8402527, 
                              scaled.oop.che = (oop.perc.che.trunc - 39.6697) / 16.0616, 
                              scaled.uhc =     (uhc.index.trunc - 61.46697) / 15.68885)
    
    covs[, "scaled.med.inc"] <- (log(medIncInput) - 7.667872) / 0.8402527
    covs[, "scaled.oop.che"] <- (oopInput - 39.6697) / 16.0616
    covs[, "scaled.uhc"]     <- (uhcInput - 61.46697) / 15.68885
    
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
plot_costs <- function(selectedCntry,
                       medIncInput,
                       oopInput,
                       uhcInput){
    
  costs_df <- calculateCosts(selectedCntry = selectedCntry, 
                             medIncInput = medIncInput, 
                             oopInput = oopInput, 
                             uhcInput = uhcInput)
  # costs_df[10,] <- c("UHC Input", uhcInput)
  # costs_df[11,] <- c("OOP Input", oopInput)
  # costs_df[12,] <- c("MI Input", medIncInput)
  # costs_df[,2] <- as.numeric(costs_df[,2])
  
  ggplot(costs_df) + 
    geom_col(aes(x=Type, y=Estimate), fill = "#A8DCAB", color="black") + 
    ggtitle(paste("Per-episode household cost of TB in", selectedCntry, "(2024 USD)")) + 
    coord_flip() + theme_minimal(base_size = 14) + 
    expand_limits(y=c(0,max(costs_df$Estimate) * 1.25)) +
  xlab("") + ylab("Costs (2024 USD)") +
    geom_text(aes(x=Type, y=Estimate, label=round(Estimate)), hjust = -.75)
             
}


##############################################################################|
##### Function to lookup relevant studies                                    ##
##############################################################################|

lookup_studies <- function(selectedCntry){
    studies_filt <- studies %>% filter(Country == selectedCntry) %>%
                    arrange(Year) %>% 
                    select(! c(Year, Country)) %>% 
                    rename(`Studies identified in systematic review` = Formatted)
    
    if (nrow(studies_filt) < 1) {
        studies_filt[1,"Studies identified in systematic review"] <- "No studies identified."
    }
    
    return(studies_filt)
}

###############################################################################
###############################################################################
###############################################################################

###############################################################################
##### user interface elements and layout
###############################################################################

ui <- page_navbar( 
        shinyjs:: useShinyjs(),
    
        id = "inTabset",
        tags$head(tags$style(HTML("a {color: black}"))),
                   
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
                      
                      tags$div(tags$p(HTML("<br><b>
                        Default covariate values for selected country: </b>"))),
                      
                      numericInput( 
                          "medInc", 
                          "Median annual household income (in 2025 INTD):", 
                          value = NULL, 
                          min = 4540, 
                          max = 81454 
                      ),
                      
                      numericInput( 
                          "oopPerc", 
                          "Out-of-pocket percentage of current health expenditures:", 
                          value = NULL, 
                          min = 4, 
                          max = 74 
                      ),
                      
                      numericInput( 
                          "uhcIndex", 
                          "Universal healthcare coverage index:", 
                          value = NULL, 
                          min = 19, 
                          max = 84 
                      ), 
                      
                      tags$div(tags$p(HTML("<br><b>
                        To calculate costs with custom covariate values, click the 'Customize covariates' button below. 
                        Update the values, then click the 'Calculate costs' button. </b>"))),
                      
                      actionButton('custCovs', 'Customize covariates',
                                   style="color: #000000; background-color: #FFFFFF"),
                  
                      actionButton('reRun',    '   Calculate costs  ',
                                   style="color: #000000; background-color: #FFFFFF"),
                      
                      
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
                                    plotOutput(outputId = "covariates")), 
                          
                          nav_panel("Relevant studies", 
                                    tableOutput(outputId = "studyTable")), 
                  )
              
              ),
              
              # HTML('<footer>
              #         CAPTURE funded by Gates
              #       </footer>')
),

tags$div(tags$p(HTML("<br>
                        Methodological overview provided in the Documentation tab.\nPredictive model code available "), tags$a(href=post, "in this repository.")))
), ##### Close PaCC nav panel 

nav_panel(title = "Documentation" , 
          setBackgroundColor("#BE91BE"),
          titlePanel(title = span(img(src = "logoZoom.png", height = 100), "Documentation and methods behind PaCC")),
          
          fluidRow(
              column(8,
                     includeMarkdown("inst/md/docAndMethods.md"))
          )
), 

nav_panel(title = "Guidance" , 
          setBackgroundColor("#BE91BE"),
          
          titlePanel(title = span(img(src = "logoZoom.png", height = 100), "Guidance and recommendations on TB patient costs")),

          fluidRow(
              column(8,
                     includeMarkdown("inst/md/guidance.md"))
          )
          )#,
# )

# tags$footer(strong("This tool was developed by CAPTURE study team."), # strong() = bold
#             align = "right",
#             style = "
#                  position:absolute;
#                  bottom:0px;
#                  width:100%;
#                  height:40px;
#                  color: black;
#                  padding: 10px;
#                  background-color: white;
#                  z-index: 0;
#                 "),
)

###############################################################################
##### server side 
###############################################################################

server <- function(input, output, session) {
    
    observeEvent(input$jumpToPaCC, {
        updateTabsetPanel(session, "inTabset",
                          selected = "PaCC")
        
        shinyjs::disable("reRun")
    })
    
    
    reactiveVar <- reactiveValues(uhcReact = NULL)
    
    #### Observe the country selected 
    observeEvent(input$country, {
        
        covariateDF <- lookup_inputs(input$country)
        
        reactiveVar$medIncReact <- round(as.numeric(covariateDF[,"median.income.trunc"]))
        reactiveVar$oopReact <- as.numeric(covariateDF[,"oop.perc.che"])
        reactiveVar$uhcReact <- as.numeric(covariateDF[,"uhc.index"])
                                           
        updateNumericInput(session,
                           inputId = "medInc",
                           value = round(as.numeric(covariateDF[,"median.income.trunc"])))
        
        shinyjs::disable("medInc")
        
        updateNumericInput(session,
                           inputId = "oopPerc",
                           value = as.numeric(covariateDF[,"oop.perc.che"]))
        
        shinyjs::disable("oopPerc")
        
        updateNumericInput(session,
                           inputId = "uhcIndex",
                           value = as.numeric(covariateDF[,"uhc.index"]))
        
        shinyjs::disable("uhcIndex")
    })

    
    ##### If the user clicks the customize covariates button, unlock the covariate boxes 
    observeEvent(input$custCovs, { 
    
        shinyjs::enable("medInc")
        shinyjs::enable("oopPerc")
        shinyjs::enable("uhcIndex")
        
        shinyjs::enable("reRun")
    
    })
    
    ##### Once the covariates are updated and the user clicks calculate costs, 
    ##### relock the input boxes. 
    
    observeEvent(input$reRun, {

        shinyjs::disable("medInc")
        shinyjs::disable("oopPerc")
        shinyjs::disable("uhcIndex")
        
        reactiveVar$medIncReact <- input$medInc
        reactiveVar$oopReact <- input$oopPerc
        reactiveVar$uhcReact <- input$uhcIndex
    })

        #### Plot country specific costs as bar plot
        output$costbar <- renderPlot({
            plot_costs(input$country, 
                       # uhcIndexCalc
                       reactiveVar$medIncReact,
                       reactiveVar$oopReact,
                       reactiveVar$uhcReact 
                       )
        })
        
        #### Present country specific costs as table
        output$costTable <- renderTable({
            calculateCosts(input$country)
        })
        
        #### Plot country specific covariates as bar plot
        output$covariates <- renderPlot({
            plot_inputs(input$country)
        })
        
        #### Present country specific studies as a table
        output$studyTable <- renderTable({
            lookup_studies(input$country)
        })
        
    
}

###############################################################################
##### Run everything 
###############################################################################

shinyApp(ui = ui, server = server)

