# preparations; required libraries
library(shiny)
library(shinyWidgets)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(bslib)


################# REMOVE THIS WHEN UNCERTAINTY VALUES ADDED 

tempUC <- rnorm(1000)

################# REMOVE THIS WHEN UNCERTAINTY VALUES ADDED

#### Repository URL 
post <- "https://github.com/CAPTURE-TB-HIV/tb-lmic-patientcosts-review"

##############################################################################|
##### Read in the required data ##
##### If need to reduce bloat, this could be optimized to only read in 
##### data for selected country -- however could make switching slower for 
##### user. 
##############################################################################|
##### External covariate data 
covariates <- readRDS("utilities/data/covariatesTime.rds")
##### Model parameter estimates 
params <- readRDS("utilities/data/parameters.rds")
params$prob_cost_zero <- 0

##### Scaling factors (mean and SD of model data) 
scalars <- readRDS("utilities/data/scalars.rds")

##### Model estimates of country random effects
cntryREs <- readRDS("utilities/data/cntryRndEff.rds") 
##### Available countries
countries <- sort(unique(covariates$study.country))

##### Relevant studies 
studies <- readRDS("utilities/data/tbPatientCostStudies.rds")


##############################################################################|
##### Function to create a bar plot of covariate inputs for selected country ##
##############################################################################|

lookup_inputs <- function(selectedCntry, 
                          selectedYear = 2024){
    covariates_filt <- covariates %>% filter(study.country == selectedCntry, 
                                             study.end.year == selectedYear) %>% 
        # mutate(scale.median.income = median.income.hh /100,
        #        scale.median.income.trunc = median.income.hh.trunc /100) %>% 
        mutate(#"Median household\nannual income\n(in 100s of 2025 INTD)" = scale.median.income.trunc, 
                "Median household\nannual income\n(2025 INTD)" = median.income.hh.trunc, 
               # "Universal healthcare\ncoverage index" = uhc.index.trunc, 
               "TB case\ndetection rate" = cdr.trunc, 
               
               "Out-of-pocket percentage\nof current health\nexpenditures" = oop.perc.che.trunc)
    
    return(covariates_filt)
}

plot_inputs <- function(selectedCntry){
    
    covariates_filt <- lookup_inputs(selectedCntry)
    
    ggplot(data = covariates_filt %>% 
               pivot_longer(cols= c(
                   `Median individual\nannual income\n(in 100s of 2025 INTD)`, 
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
                           yearInput = 2024,
                           medIncInput, 
                           oopInput,
                           cdrInput) {
    
    covs <- covariates %>% 
        filter(study.country == selectedCntry,
               study.end.year == yearInput)# %>% 
        # mutate(scaled.med.inc = (log.median.income.hh.trunc - scalars["mean", "log.median.income.hh"]) / scalars["sd", "log.median.income.hh"], 
        #        scaled.oop.che = (oop.perc.che.trunc - scalars["mean", "oop.perc.che"]) / scalars["sd", "oop.perc.che"], 
        #        scaled.cdr =     (cdr.trunc - scalars["mean", "cdr.rate"]) / scalars["sd", "cdr.rate"])
    
    covs[, "scaled.med.inc"] <- (log(medIncInput) - scalars["mean", "log.median.income.hh"]) / scalars["sd", "log.median.income.hh"] 
    covs[, "scaled.oop.che"] <- (oopInput - scalars["mean", "oop.perc.che"]) / scalars["sd", "oop.perc.che"]
    covs[, "scaled.cdr"]     <- (cdrInput - scalars["mean", "cdr.rate"]) / scalars["sd", "cdr.rate"]

    # covs[, "scaled.uhc"]     <- (uhcInput - 61.46697) / 15.68885
    
    cntryRE2 <- cntryREs %>% filter(Country == selectedCntry) %>% select(!Country) %>% 
        unlist() %>% as.numeric()
    
    scaledYear <- (yearInput - scalars["mean", "study.end.year"]) / scalars["sd", "study.end.year"]
    
    estimates <- data.frame(dme.prediag = 
                                exp(params[,"beta_dme_prediag"] + cntryRE2 + 
                                    params[,"mult_med_income_dme_prediag"] * covs$scaled.med.inc +
                                    params[,"mult_oop_perc_che_dme_prediag"] * covs$scaled.oop.che +
                                    params[,"mult_cdr_dme_prediag"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_dme_prediag"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]),
                            ### DNME PREDIAGNOSIS 
                            dnme.prediag =   
                                exp(params[,"beta_dnme_prediag"] + cntryRE2 + 
                                        params[,"mult_med_income_dnme_prediag"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_dnme_prediag"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_dnme_prediag"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_dnme_prediag"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]), 
                            ie.prediag =   
                            ### INDIRECT PREDIAGNOSIS 
                            exp(params[,"beta_ie_prediag"] + cntryRE2 + 
                                    params[,"mult_med_income_ie_prediag"] * covs$scaled.med.inc +
                                    params[,"mult_oop_perc_che_ie_prediag"] * covs$scaled.oop.che +
                                    params[,"mult_cdr_ie_prediag"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_ie_prediag"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]), 
    
                            dme.intsvtx = 
                                exp(params[,"beta_dme_intsvtx"] + cntryRE2 + 
                                        params[,"mult_med_income_dme_intsvtx"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_dme_intsvtx"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_dme_intsvtx"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_dme_intsvtx"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]),
                            ### DNME INTENSIVE TX 
                            dnme.intsvtx =   
                                exp(params[,"beta_dnme_intsvtx"] + cntryRE2 + 
                                        params[,"mult_med_income_dnme_intsvtx"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_dnme_intsvtx"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_dnme_intsvtx"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_dnme_intsvtx"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]), 
                            ### INDIRECT INTENSIVE TX 
                            ie.intsvtx =   
                                exp(params[,"beta_ie_intsvtx"] + cntryRE2 + 
                                        params[,"mult_med_income_ie_intsvtx"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_ie_intsvtx"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_ie_intsvtx"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_ie_intsvtx"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]),
                            
                            
                            dme.contntx = 
                                exp(params[,"beta_dme_contntx"] + cntryRE2 + 
                                        params[,"mult_med_income_dme_contntx"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_dme_contntx"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_dme_contntx"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_dme_contntx"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]),
                            ### DNME contntxNOSIS 
                            dnme.contntx =   
                                exp(params[,"beta_dnme_contntx"] + cntryRE2 + 
                                        params[,"mult_med_income_dnme_contntx"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_dnme_contntx"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_dnme_contntx"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_dnme_contntx"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]), 
                            ### INDIRECT CONTINUATION TX 
                            ie.contntx =   
                                exp(params[,"beta_ie_contntx"] + cntryRE2 + 
                                        params[,"mult_med_income_ie_contntx"] * covs$scaled.med.inc +
                                        params[,"mult_oop_perc_che_ie_contntx"] * covs$scaled.oop.che +
                                        params[,"mult_cdr_ie_contntx"] * covs$scaled.cdr) *
                                exp(params[,"mult_year_ie_contntx"] * scaledYear) *       
                                (1-params[,"prob_cost_zero"]))
    
    allCosts <- data.frame(Category = c(rep("Pre-treatment", 3), 
                                        rep("Intensive\ntreatment", 3),
                                        rep("Continuation\ntreatment", 3)), 
                           Phase =   rep(c("Direct\nmedical", "Direct\nnon-medical", "Indirect"), 3), 
                           Estimate = round(apply(estimates, MARGIN = 2, mean)),
                           LB = round(apply(estimates, MARGIN = 2, quantile, probs=0.025)),
                           UB = round(apply(estimates, MARGIN = 2, quantile, probs=0.975))) %>%
                 mutate(`95% uncertainty interval` = paste0("(", LB, ", ", UB, ")"))
    
    return(allCosts)
}


##############################################################################|
##### Function to create a bar plot of calculated costs for selected country ##
##############################################################################|

plot_stackP_costs <- function(selectedCntry,
                              costs_df
){
    ggplot(costs_df %>% 
               dplyr::mutate(Phase = factor(Phase, levels = c( "Indirect", 
                                                               "Direct\nnon-medical", 
                                                               "Direct\nmedical"
               ), ordered=TRUE))) + 
        geom_col(aes(x=Phase, y=Estimate, fill = Category), color="black", position = "stack") + 
        ggtitle(paste("Per-episode household cost of TB in", selectedCntry, "by cost category")) + 
        coord_flip() +
        theme_minimal(base_size = 18) + 
        expand_limits(y=c(0,max(costs_df$Estimate) * 1.25)) +
        scale_fill_manual(values = c("#60697C", "#94AEE9","#B5C7EA"), 
                          guide = guide_legend(reverse = TRUE)) + 
        xlab("") + ylab("Costs (2025 INTD)") +
        theme(legend.position = "bottom",
              plot.title = element_text(size = 18), axis.text = element_text(size = 16)) #+ 
    # geom_text(aes(x=Phase, y=Estimate, label=round(Estimate)), position = position_stack(0.5))
}

plot_stackC_costs <- function(selectedCntry,
                              costs_df
){
    ggplot(costs_df %>% 
               dplyr::mutate(Phase = factor(Phase, levels = c( "Indirect", 
                                                               "Direct\nnon-medical", 
                                                               "Direct\nmedical"
               ), ordered=TRUE))) + 
        geom_col(aes(x=Category, y=Estimate, fill = Phase), color="black") + 
        ggtitle(paste("Per-episode household cost of TB in", selectedCntry, "by treatment phase")) + 
        coord_flip() + theme_minimal(base_size = 18) + 
        expand_limits(y=c(0,max(costs_df$Estimate) * 1.25)) +
        scale_fill_manual(values = c("#2F3D2F", "#4E664F", "#A8DCAB"), 
                          guide = guide_legend(reverse = TRUE)) + # "#8CB88E",
        xlab("") + ylab("Costs (2025 INTD)")+
        theme(legend.position = "bottom", 
              plot.title = element_text(size = 18), axis.text = element_text(size = 16))
}


plot_costs <- function(selectedCntry,
                       costs_df
){
    ### CREATE THE TYPE VARIABLE 
    costs_dfB <- costs_df %>% 
        mutate(Type = c("Pre-treatment\ndirect medical", 
                        "Pre-treatment\ndirect non-medical",
                        "Pre-treatment\nindirect",
                        
                        "Intensive treatment\ndirect medical", 
                        "Intensive treatment\ndirect non-medical",
                        "Intensive treatment\nindirect", 
                        
                        "Continuation treatment\ndirect medical", 
                        "Continuation treatment\ndirect non-medical",
                        "Continuation treatment\nindirect")) 
    
    
    ggplot(costs_dfB) + 
        geom_col(aes(x=Type, y=Estimate), fill = "#A8DCAB", color="black") + 
        ggtitle(paste("Per-episode household cost of TB in", selectedCntry, "(2025 INTD)")) + 
        coord_flip() + theme_minimal(base_size = 14) + 
        expand_limits(y=c(0,max(costs_df$Estimate) * 1.25)) +
        xlab("") + ylab("Costs (2025 INTD)") +
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
              p("The underlying model is informed from data collected from 104 studies of TB patient costs, representing more than 55,000 persons treated for TB. The model uses routinely reported country- and year- matched data to estimates the cost of TB, stratified by three cost Types (direct medical, direct non-medical, and indirect) and across the phases of the TB care cascade (pre-treatment, intensive treatment, and continuation treatment). All estimates represent TB patient costs in 2024 in 2025 INTD."),
              p("Please note: TB patient cost data varies substantially across settings and study design. As such, our estimates should be used as a guide point for analyses."), 
              
              actionButton('jumpToPaCC', 'Get started!', icon("paper-plane"),
                           style="color: #000000; background-color: #FFFFFF")),
    
    nav_panel(title = "PaCC" , 
              setBackgroundColor("#BE91BE"),
              
              # titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
              
              ##### Sidebar which includes the country dropdown selection #####
              sidebarLayout(
                  
                  sidebarPanel(
                      titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
                      
                      tags$div(tags$p(HTML("<b>
                        Select a country and year below and the per-episode household costs of TB will be automatically calculated.</b>" ))),
                      
                      selectInput(
                          "country",
                          label = "Select a country",
                          choices = countries
                      ),
                      
                      selectInput(
                          "year",
                          selected = 2024,
                          label = "Select a year",
                          choices = 2009:2024
                      ),
                      
                      # checkboxGroupInput( 
                      #     "drOpt", 
                      #     "Resistance status", 
                      #     c( 
                      #         "Drug-sensitive TB" = "dstb", 
                      #         "RR/MDR TB" = "drtb" 
                      #     ) 
                      # ), 
                      
                      hr(),
                      tags$div(tags$p(HTML("<b>
                        Default covariate values for selected country and year: </b>"))),
                      
                      # numericInput( 
                      #     "year", 
                      #     "Year:", 
                      #     value = NULL, 
                      #     min = 2009, 
                      #     max = 2024 
                      # ),
                      
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
                          "cdr", 
                          "TB case detection rate", 
                          value = NULL, 
                          min = 24, 
                          max = 99 
                      ), 
                      
                      tags$div(tags$p(HTML("<br><b>
                        To calculate costs with custom covariate values, click the 'Customize covariates' button below. 
                        Update the values, then click the 'Calculate costs' button to recalculate costs. </b>"))),
                      
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
                          nav_panel("Cost figures", 
                                    plotOutput(outputId = "costCstack", height = "490px"),
                                    plotOutput(outputId = "costPstack", height = "490px"),
                                    
                                    downloadButton(
                                        "downloadFigure",
                                        label = "Download cost figures",
                                        icon = shiny::icon("download"),
                                        enabled = TRUE
                                    )
                                    # plotOutput(outputId = "costbar")
                          ),
                          ##### Table of costs #####
                          nav_panel("Cost table",
                                    tableOutput(outputId = "costTable"),
                                    downloadButton(
                                        "downloadTable",
                                        label = "Download cost table",
                                        icon = shiny::icon("download"),
                                        enabled = TRUE
                                    )),
                          ##### Plot of covariates #####
                          # nav_panel("Covariate plots", 
                          #           plotOutput(outputId = "covariates")), 
                          
                          nav_panel("Country-specific cost studies", 
                                    tableOutput(outputId = "studyTable")), 
                      ), 
                      
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
    
    ### Always start at 2024 
    baseYear <- 2024
    
    reactiveVar <- reactiveValues(cdrReact = NULL, oopReact = NULL, 
                                  medIncReact = NULL, yearReact = NULL, 
                                  calcCosts = NULL, costPlot = NULL, 
                                  costPlotStckP = NULL, costPlotStckC = NULL)
    
    #### Observe the country selected 
    observeEvent(input$country, {
        
        covariateDF <- lookup_inputs(selectedCntry = input$country, 
                                     selectedYear = baseYear)
        
        reactiveVar$medIncReact <- round(as.numeric(covariateDF[,"median.income.hh.trunc"]))
        reactiveVar$oopReact <- as.numeric(covariateDF[,"oop.perc.che.trunc"])
        reactiveVar$cdrReact <- as.numeric(covariateDF[,"cdr.trunc"])
        # reactiveVar$cdrReact <- as.numeric(covariateDF[,"uhc.index"])
        reactiveVar$yearReact <- baseYear
        
        updateNumericInput(session,
                           inputId = "medInc",
                           value = round(as.numeric(covariateDF[,"median.income.hh.trunc"])))
        
        shinyjs::disable("medInc")
        
        updateNumericInput(session,
                           inputId = "oopPerc",
                           value = as.numeric(covariateDF[,"oop.perc.che.trunc"]))
        
        shinyjs::disable("oopPerc")
        
        updateNumericInput(session,
                           inputId = "cdr",
                           value = as.numeric(covariateDF[,"cdr.trunc"]))
        
        shinyjs::disable("cdr")
        
        
        updateSelectInput(
            session,
            inputId = "year",
            selected = baseYear
        )
        
        # print(input$country)
        # print(reactiveVar$yearReact)
        # print(reactiveVar$medIncReact)
        # print(reactiveVar$oopReact)
        # print(reactiveVar$cdrReact)
        
        reactiveVar$calcCosts <- calculateCosts(selectedCntry = input$country,
                                                yearInput = reactiveVar$yearReact,
                                                medIncInput = reactiveVar$medIncReact,
                                                oopInput = reactiveVar$oopReact,
                                                cdrInput = reactiveVar$cdrReact)
    })
    
    observeEvent(input$year, {
        reactiveVar$yearReact <- as.numeric(input$year) 
        
        covariateDF <- lookup_inputs(selectedCntry = input$country, 
                                     selectedYear = as.numeric(input$year))
        
        reactiveVar$medIncReact <- round(as.numeric(covariateDF[,"median.income.trunc"]))
        reactiveVar$oopReact <- as.numeric(covariateDF[,"oop.perc.che.trunc"])
        reactiveVar$cdrReact <- as.numeric(covariateDF[,"cdr.trunc"])
        
        updateNumericInput(session,
                           inputId = "medInc",
                           value = round(as.numeric(covariateDF[,"median.income.trunc"])))
        
        updateNumericInput(session,
                           inputId = "oopPerc",
                           value = as.numeric(covariateDF[,"oop.perc.che.trunc"]))
        
        updateNumericInput(session,
                           inputId = "cdr",
                           value = as.numeric(covariateDF[,"cdr.trunc"]))
        
        # print(input$country)
        # print(reactiveVar$yearReact)
        # print(reactiveVar$medIncReact)
        # print(reactiveVar$oopReact)
        # print(reactiveVar$cdrReact)
        
        reactiveVar$calcCosts <- calculateCosts(selectedCntry = input$country,
                                                yearInput = reactiveVar$yearReact,
                                                medIncInput = reactiveVar$medIncReact,
                                                oopInput = reactiveVar$oopReact,
                                                cdrInput = reactiveVar$cdrReact)
    })
    
    ##### If the user clicks the customize covariates button, unlock the covariate boxes 
    observeEvent(input$custCovs, { 
        
        shinyjs::enable("medInc")
        shinyjs::enable("oopPerc")
        shinyjs::enable("cdr")
        
        shinyjs::enable("reRun")
        
    })
    
    ##### Once the covariates are updated and the user clicks calculate costs, 
    ##### relock the input boxes. 
    
    observeEvent(input$reRun, {
        
        shinyjs::disable("medInc")
        shinyjs::disable("oopPerc")
        shinyjs::disable("cdr")
        
        reactiveVar$medIncReact <- input$medInc
        reactiveVar$oopReact <- input$oopPerc
        reactiveVar$cdrReact <- input$cdr
        
        reactiveVar$calcCosts <- calculateCosts(selectedCntry = input$country,
                                                yearInput = reactiveVar$yearReact,
                                                medIncInput = reactiveVar$medIncReact,
                                                oopInput = reactiveVar$oopReact,
                                                cdrInput = reactiveVar$cdrReact)  
        
        shinyjs::disable("reRun")
    })
    
    #### Plot country specific costs as bar plots
    output$costPstack <- renderPlot({
        reactiveVar$costPlotStckP <- plot_stackP_costs(input$country, 
                                                       reactiveVar$calcCosts
        )
        reactiveVar$costPlotStckP
    })
    
    output$costCstack <- renderPlot({
        reactiveVar$costPlotStckC <- plot_stackC_costs(input$country, 
                                                       reactiveVar$calcCosts
        )
        reactiveVar$costPlotStckC
    })
    
    # output$costbar <- renderPlot({
    #     reactiveVar$costPlot <- plot_costs(input$country, 
    #                            reactiveVar$calcCosts
    #     )
    #     reactiveVar$costPlot
    # })
    
    #### Present country specific costs as table
    output$costTable <- renderTable({
        reactiveVar$calcCosts %>% 
            select(! c(LB,UB))
    }, digits = 0)
    
    #### Plot country specific covariates as bar plot
    output$covariates <- renderPlot({
        plot_inputs(input$country)
    })
    
    #### Present country specific studies as a table
    output$studyTable <- renderTable({
        lookup_studies(input$country)
    })
    
    
    output$downloadFigure <- downloadHandler(
        filename = function() {
            # Use the selected dataset as the suggested file name
            paste0("PaCC_estimate_",  input$country, ".pdf")
        },
        content = function(file) {
            # Write the dataset to the `file` that will be downloaded
            pList <- list(reactiveVar$costPlotStckP, 
                          reactiveVar$costPlotStckC)
            
            pdf(file = file, onefile = TRUE, width=11, height=8.5)
            for(i in 1:length(pList)){
                print(pList[[i]])
            }
            dev.off()
        }
    )
    
    output$downloadTable <- downloadHandler(
        filename = function() {
            # Use the selected dataset as the suggested file name
            paste0("PaCC_estimate_", input$country, ".csv")
        },
        content = function(file) {
            # Write the dataset to the `file` that will be downloaded
            write.csv(reactiveVar$calcCosts, file)
        }
    )
    
    
}

###############################################################################
##### Run everything 
###############################################################################

shinyApp(ui = ui, server = server)
