# preparations; required libraries
library(shiny)
library(shinyWidgets)
library(bslib)
library(dplyr)
library(ggplot2)
##### This load is required for webR/shinylive to work 
if (FALSE) {
    library(munsell)
}
#### YOU MUST COMMENT THIS OUT IF YOU WANT TO RUN THE SHINY APP LIKE NORMAL
webr::install("markdown")

#### Repository URL 
post <- "https://github.com/CAPTURE-TB-HIV/tb-lmic-patientcosts-review"

##############################################################################|
##### Read in the required data ##
##### If need to reduce bloat, this could be optimized to only read in 
##### data for selected country -- however could make switching slower for 
##### user. 
##############################################################################|
##### External covariate data 
tb_covariates <- readRDS("utilities/data/TB/covariatesTimeSmall.rds")
# hiv_covariates <- readRDS("utilities/data/HIV/covariatesTime.rds")

##### Comparison data 
tb_comp_df <- readRDS("utilities/data/TB/compData.rds")
# hiv_comp_df <- readRDS("utilities/data/HIV/compData.rds")

##### Model parameter estimates 
tb_params <- readRDS("utilities/data/TB/parameters.rds")

# hiv_params <- readRDS("utilities/data/HIV/parameters.rds")

##### Scaling factors (mean and SD of model data) 
tb_scalars <- readRDS("utilities/data/TB/scalars.rds")
# hiv_scalars <- readRDS("utilities/data/HIV/scalars.rds")

##### Model estimates of country random effects
tb_cntryREs <- readRDS("utilities/data/TB/cntryRndEff.rds") 
# tb_cntryREs %>% filter(Country == "China") %>% select(!Country) %>% as.numeric() %>% mean()
# hiv_cntryREs <- readRDS("utilities/data/HIV/cntryRndEff.rds") 

##### Available countries
countries <- sort(unique(tb_covariates$study.country))

##### Relevant studies 
studies <- readRDS("utilities/data/tbPatientCostStudies.rds")


##############################################################################|
##### Function to create a bar plot of covariate inputs for selected country ##
##############################################################################|

lookup_inputs <- function(selectedCntry, 
                          selectedYear = 2024, 
                          infection = "TB"){
    
    covariates_filt <- switch(infection, 
                        "TB" = tb_covariates %>% filter(study.country == selectedCntry, 
                                                         study.end.year == selectedYear) %>% 
                            # mutate(scale.median.income = median.income.hh /100,
                            #        scale.median.income = median.income.hh /100) %>% 
                            mutate(#"Median household\nannual income\n(in 100s of 2025 INTD)" = scale.median.income, 
                                "Median household\nannual income\n(2025 INTD)" = median.income.hh, 
                                # "Universal healthcare\ncoverage index" = uhc.index, 
                                "TB case\ndetection rate" = cdr.rate, 
                                "Out-of-pocket percentage\nof current health\nexpenditures" = oop.perc.che), 
                        "HIV" =  hiv_covariates %>% filter(study.country == selectedCntry, 
                                                          study.end.year == selectedYear) %>% 
                            mutate("Median household\nannual income\n(2025 INTD)" = median.income, 
                                   "Universal healthcare\ncoverage index" = uhc.index,
                                   "Out-of-pocket percentage\nof current health\nexpenditures" = oop.perc.che)
    )

    return(covariates_filt)
}

##############################################################################|
##### Function to calculate TB costs for selected country & return dataframe ##
##############################################################################|

tb_calculateCosts <- function(selectedCntry, 
                           yearInput = 2024,
                           medIncInput, 
                           oopInput,
                           cdrInput) {
    
    covs <- tb_covariates %>% 
        filter(study.country == selectedCntry,
               study.end.year == yearInput)
    
    covs[, "scaled.med.inc"] <- (log(medIncInput) - tb_scalars["mean", "log.median.income.hh"]) / tb_scalars["sd", "log.median.income.hh"] 
    covs[, "scaled.oop.che"] <- (oopInput - tb_scalars["mean", "oop.perc.che"]) / tb_scalars["sd", "oop.perc.che"]
    covs[, "scaled.cdr"]     <- (cdrInput - tb_scalars["mean", "cdr.rate"]) / tb_scalars["sd", "cdr.rate"]
    
    cntryRE2 <- tb_cntryREs %>% filter(Country == selectedCntry) %>% select(!Country) %>% 
        unlist() %>% as.numeric()
    scaledYear <- (yearInput - tb_scalars["mean", "study.end.year"]) / tb_scalars["sd", "study.end.year"]
    
    estimates <- data.frame(dme.prediag = 
                                exp(tb_params[,"beta_dme_prediag"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_dme_prediag"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_dme_prediag"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_dme_prediag"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_dme_prediag"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_dme_prediag"])),
                            ### DNME PREDIAGNOSIS 
                            dnme.prediag =   
                                exp(tb_params[,"beta_dnme_prediag"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_dnme_prediag"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_dnme_prediag"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_dnme_prediag"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_dnme_prediag"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_dnme_prediag"])), 
                            ie.prediag =   
                                ### INDIRECT PREDIAGNOSIS 
                                exp(tb_params[,"beta_ie_prediag"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_ie_prediag"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_ie_prediag"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_ie_prediag"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_ie_prediag"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_ie_prediag"])), 
                            
                            dme.intsvtx = 
                                exp(tb_params[,"beta_dme_intsvtx"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_dme_intsvtx"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_dme_intsvtx"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_dme_intsvtx"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_dme_intsvtx"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_dme_intsvtx"])),
                            ### DNME INTENSIVE TX 
                            dnme.intsvtx =   
                                exp(tb_params[,"beta_dnme_intsvtx"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_dnme_intsvtx"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_dnme_intsvtx"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_dnme_intsvtx"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_dnme_intsvtx"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_dnme_intsvtx"])), 
                            ### INDIRECT INTENSIVE TX 
                            ie.intsvtx =   
                                exp(tb_params[,"beta_ie_intsvtx"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_ie_intsvtx"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_ie_intsvtx"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_ie_intsvtx"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_ie_intsvtx"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_ie_intsvtx"])),
                            
                            
                            dme.contntx = 
                                exp(tb_params[,"beta_dme_contntx"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_dme_contntx"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_dme_contntx"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_dme_contntx"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_dme_contntx"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_dme_contntx"])),
                            ### DNME contntxNOSIS 
                            dnme.contntx =   
                                exp(tb_params[,"beta_dnme_contntx"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_dnme_contntx"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_dnme_contntx"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_dnme_contntx"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_dnme_contntx"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_dnme_contntx"])), 
                            ### INDIRECT CONTINUATION TX 
                            ie.contntx =   
                                exp(tb_params[,"beta_ie_contntx"] + cntryRE2 + 
                                        tb_params[,"mult_med_income_ie_contntx"] * covs$scaled.med.inc +
                                        tb_params[,"mult_oop_perc_che_ie_contntx"] * covs$scaled.oop.che +
                                        tb_params[,"mult_cdr_ie_contntx"] * covs$scaled.cdr) *
                                exp(tb_params[,"mult_year_ie_contntx"] * scaledYear) *       
                                (1-exp(tb_params[,"beta_theta_ie_contntx"])))
    
    estimates$Total <- rowSums(estimates)
    
    allCosts <- data.frame(Category = c(rep("Pre-treatment", 3), 
                                        rep("Intensive\ntreatment", 3),
                                        rep("Continuation\ntreatment", 3), 
                                        ""), 
                           Phase =   c(rep(c("Direct\nmedical", "Direct\nnon-medical", "Indirect"), 3), ""), 
                           Description = c("Pre-treatment\ndirect medical", 
                                           "Pre-treatment\ndirect non-medical",
                                           "Pre-treatment\nindirect",
                                           
                                           "Intensive treatment\ndirect medical", 
                                           "Intensive treatment\ndirect non-medical",
                                           "Intensive treatment\nindirect", 
                                           
                                           "Continuation treatment\ndirect medical", 
                                           "Continuation treatment\ndirect non-medical",
                                           "Continuation treatment\nindirect", 
                                           "Total"),  
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
                              costs_df, 
                              infection = "TB"
){
    plotPal <- switch(infection, 
                  TB = c("#60697C", "#94AEE9","#B5C7EA"), 
                  HIV = "#B5C7EA")
    
    plot.opts <- switch(infection, 
                        TB = c("bottom", "Per-episode"), 
                        HIV = c("none", "Annual"))
    
    infection <- ifelse(infection == "HIV", 
                        "ART", "TB")
    

    ggplot(costs_df %>% 
               filter(Description != "Total") %>%
               dplyr::mutate(Phase = factor(Phase, levels = c( "Indirect", 
                                                               "Direct\nnon-medical", 
                                                               "Direct\nmedical"
               ), ordered=TRUE))) + 
        geom_col(aes(x=Phase, y=Estimate, fill = Category), color="black", position = "stack") + 
        ggtitle(paste(plot.opts[2], "shousehold cost of", infection, "in", selectedCntry, "by cost category")) + 
        coord_flip() +
        theme_minimal(base_size = 18) + 
        expand_limits(y=c(0,max(costs_df$Estimate[-10]) * 1.25)) +
        scale_fill_manual(values = plotPal, 
                          guide = guide_legend(reverse = TRUE)) + 
        xlab("") + ylab("Costs (2025 INTD)") +
        theme(legend.position = plot.opts[1],
              plot.title = element_text(size = 18), axis.text = element_text(size = 16)) #+ 
    # geom_text(aes(x=Phase, y=Estimate, label=round(Estimate)), position = position_stack(0.5))
}

plot_stackC_costs <- function(selectedCntry,
                              costs_df
){
    ggplot(costs_df %>% 
               filter(Description != "Total") %>%
               dplyr::mutate(Phase = factor(Phase, levels = c( "Indirect", 
                                                               "Direct\nnon-medical", 
                                                               "Direct\nmedical"
               ), ordered=TRUE))) + 
        geom_col(aes(x=Category, y=Estimate, fill = Phase), color="black") + 
        ggtitle(paste("Per-episode household cost of TB in", selectedCntry, "by treatment phase")) + 
        coord_flip() +
        theme_minimal(base_size = 18) + 
        expand_limits(y=c(0,max(costs_df$Estimate[-10]) * 1.25)) +
        scale_fill_manual(values = c("#2F3D2F", "#4E664F", "#A8DCAB"), 
                          guide = guide_legend(reverse = TRUE)) + # "#8CB88E",
        xlab("") + ylab("Costs (2025 INTD)")+
        theme(legend.position = "bottom", 
              plot.title = element_text(size = 18), axis.text = element_text(size = 16))
}

plot_comp_costs <- function(selectedCntry,
                            costs_df, 
                            infection = "TB"
){
    comp_df <- switch(infection, 
                      "TB" = tb_comp_df, 
                      "HIV" = hiv_comp_df)
    
    covariates <- switch(infection, 
                         "TB" = tb_covariates, 
                         "HIV" = hiv_covariates)
    
    infection <- ifelse(infection == "HIV", 
                        "ART", "TB")
    
    region <- covariates %>% ungroup() %>% 
        filter(study.country == selectedCntry) %>% 
        select(who.region) %>% unique() %>% unlist()
    
    summaryCosts <- comp_df %>% 
        filter(geography %in% c("LMICs", region)) %>% 
        rename(Estimate = est) %>% 
        mutate(Description = ifelse(geography == "LMICs", 
                                    "All LMICs\naverage", 
                                    "Regional\naverage")) %>% 
        select(Description, Estimate) %>% 
        rbind(costs_df %>% filter(Description == "Total") %>%
                  select(Description, Estimate) %>% 
                  mutate(Description = paste(selectedCntry))) %>% 
        mutate(color = ifelse(Description == paste(selectedCntry),
                              "dark", "light"), 
               Description = factor(Description, levels = c("All LMICs\naverage", "Regional\naverage", 
                                                            paste(selectedCntry))))
    ggplot(summaryCosts) + 
        geom_col(aes(x=Description, y=Estimate, fill=color), color="black") + 
        ggtitle(paste("Comparison of household cost of", infection, "in", selectedCntry, "with regional and all LMIC average")) + 
        coord_flip() + theme_minimal(base_size = 18) + 
        expand_limits(y=c(0,max(summaryCosts$Estimate[-10]) * 1.25)) +
        scale_fill_manual(values = c("#A567A5", "#D7BBD7")) + # "#8CB88E",
        xlab("") + ylab("Costs (2025 INTD)")+
        theme(legend.position = "none", 
              plot.title = element_text(size = 18), axis.text = element_text(size = 16))  +
        geom_text(aes(x=Description, y=Estimate, label=round(Estimate)), hjust = -.5, fontface="bold")
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
                  # p("TB Patient Cost Calculator (TB PaCC) and HIV Patient Cost Calculator (HIV PaCC)", style="color: #000000; font-size: 36px;"),
                  p("TB Patient Cost Calculator (TB PaCC)", style="color: #000000; font-size: 36px;"),
                  
                  # HTML('<center> <img src="logoZoom.png" height="500"  /><img src="HIV_logoZoom.png" height="500"  />')
                  
              ),
              
              setBackgroundColor("#BE91BE"),

              
              fluidRow(
                  # p("The TB Patient Cost Calculator (TB PaCC) and HIV Patient Cost Calculator (HIV PaCC) estimate country-level patient costs for 129 low- and middle income countries. We chose the name PaCC to represent the significant planning and resources required among persons undergoing TB and HIV treatment. Furthermore, the economic burden of TB treatment can place financial strain on households, like carrying a heavy pack throughout a journey.", style="color: #000000; font-size: 24px;"),
                  p("The TB Patient Cost Calculator (TB PaCC) estimates country-level TB household costs for 129 low- and middle income countries. We chose the name PaCC to represent the significant planning and resources required among persons undergoing TB treatment. Furthermore, the economic burden of TB treatment can place financial strain on households, like carrying a heavy pack throughout a journey.", style="color: #000000; font-size: 24px;"),
                  
              ),
              fluidRow(
                  # column(width=1.5),
                  column(
                      width = 5,
                      offset = 0,
                      HTML('<center> <img src="logoZoom.png" height="500"  />'),
                      
                      p("The TB PaCC model is informed from data collected from 104 studies of TB patient costs, representing more than 55,000 persons treated for TB. The model uses routinely reported country- and year- matched data to estimates the cost of TB, stratified by three cost types (direct medical, direct non-medical, and indirect) and across the phases of the TB care cascade (pre-treatment, intensive treatment, and continuation treatment). All estimates represent TB patient costs in constant 2025 International Dollars", style="color: #000000; font-size: 18px;"),
                      p("Please note: TB patient cost data varies substantially across settings and study design. As such, our estimates should be used as a guide point for analyses.", style="color: #000000; font-size: 18px;"),
                      
              actionButton('jumpToTBPaCC', 'Calculate TB costs!', icon("paper-plane"),
                           style="color: #000000; background-color: #FFFFFF; font-size: 20px; inline-block; margin-left:20%; margin-right:20%")
              ), 
              # column(width=1),
              # column(
              #     width = 5,
              #     offset = 0,
              #     HTML('<center> <img src="HIV_logoZoom.png" height="500"  />'),
              #     p("The HIV PaCC model is informed from data collected from 32 studies of HIV public ART costs. The model uses routinely reported country- and year- matched data to estimates the cost of HIV, stratified by three cost types (direct medical, direct non-medical, and indirect). All estimates represent ART patient costs in constant 2025 International Dollars.", style="color: #000000; font-size: 18px;"),
              #     
              # actionButton('jumpToHIVPaCC', 'Calculate HIV costs!', icon("paper-plane"),
              #              style="color: #000000; background-color: #FFFFFF; font-size: 20px; display: inline-block; margin-left:20%; margin-right:20%")
              # ), 
              # column(width=1),
              style= "margin-left:15%;  margin-right:15%;"
              
              ),

    ),
    ###########################################################################
    ###### TB PaCC webtool tab
    ###########################################################################
    nav_panel(title = "TBPaCC" , 
              setBackgroundColor("#BE91BE"),
              
              # titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
              
              ##### Sidebar which includes the country dropdown selection #####
              sidebarLayout(
                  
                  sidebarPanel(
                      titlePanel(title = span(img(src = "logoZoom.png", height = 100), "TB Patient Cost Calculator")),
                      
                      tags$div(tags$p(HTML("<b>
                        Select a country and year below and the per-episode household costs of TB will be automatically calculated.</b>" ))),
                      
                      selectInput(
                          "tb_country",
                          label = "Select a country",
                          choices = countries
                      ),
                      
                      selectInput(
                          "tb_year",
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
                      
                      numericInput( 
                          "tb_medInc", 
                          "Median annual household income (in 2025 INTD):", 
                          value = NULL, 
                          min = 4540, 
                          max = 81454 
                      ),
                      
                      numericInput( 
                          "tb_oopPerc", 
                          "Out-of-pocket percentage of current health expenditures:", 
                          value = NULL, 
                          min = 4, 
                          max = 74 
                      ),
                      
                      numericInput( 
                          "tb_cdr", 
                          "TB case detection rate", 
                          value = NULL, 
                          min = 24, 
                          max = 99 
                      ), 
                      
                      tags$div(tags$p(HTML("<br><b>
                        To calculate costs with custom covariate values, click the 'Customize tb_covariates' button below. 
                        Update the values, then click the 'Calculate costs' button to recalculate costs. </b>"))),
                      
                      actionButton('tb_custCovs', 'Customize covariates',
                                   style="color: #000000; background-color: #FFFFFF"),
                      
                      actionButton('tb_reRun',    '   Calculate costs  ',
                                   style="color: #000000; background-color: #FFFFFF"),
                      
                      
                  ),
                  ##### Main panel of results on right of screen
                  mainPanel(
                      
                      ##### Create an inset of tabs for the user to click through #####
                      navset_card_underline(
                      title = "",
                          ##### Plot of costs #####
                          nav_panel("Cost figures",
                                    # splitLayout(cellWidths = c("50%", "50%"), plotOutput("costCstack", height="600px"), plotOutput("costPstack", height="600px")),

                                    plotOutput(outputId = "tb_costCstack", height = "490px"),
                                    plotOutput(outputId = "tb_costPstack", height = "490px"),
                                    plotOutput(outputId = "tb_costCompPlot", height = "490px"),


                                    downloadButton(
                                        "tb_downloadFigure",
                                        label = "Download cost figures",
                                        icon = shiny::icon("download"),
                                        enabled = TRUE
                                    )
                          ),
                          ##### Table of costs #####
                          nav_panel("Cost table",
                                    tableOutput(outputId = "tb_costTable"),
                                    # downloadButton(
                                    #     "tb_downloadTable",
                                    #     label = "Download cost table",
                                    #     icon = shiny::icon("download"),
                                    #     enabled = TRUE)
                                    ),
                          
                          nav_panel("Country-specific cost studies", 
                                    tableOutput(outputId = "tb_studyTable")), 
                      ), ##### Close navset card 
                  ), ##### Close main panel 
              ), ##### Close sidebar layout 
              
              tags$div(tags$p(HTML("<br>
                        Methodological overview provided in the Documentation tab.\nPredictive model code available "), tags$a(href=post, "in this repository.")))
    ), ##### Close TBPaCC nav panel
    
    ###########################################################################
    ###### Documentation tab which renders a markdown document
    ###########################################################################
    nav_panel(title = "Documentation" , 
              setBackgroundColor("#BE91BE"),
              titlePanel(title = span(img(src = "logoZoom.png", height = 100), 
                                      "Documentation and methods behind PaCC")),
              
              fluidRow(
                  column(8,
                         includeMarkdown("inst/md/docAndMethods.md"))
              )
    ), 
    
    ###########################################################################
    ###### Guidance tab which renders a markdown document
    ###########################################################################
    nav_panel(title = "Guidance" , 
              setBackgroundColor("#BE91BE"),
              
              titlePanel(title = span(img(src = "logoZoom.png", height = 100), "Guidance and recommendations on TB patient costs")),
              
              fluidRow(
                  column(8,
                         includeMarkdown("inst/md/guidance.md"))
              )
    )#,
    # )
)

###############################################################################
##### server side 
###############################################################################

server <- function(input, output, session) {

        observeEvent(input$jumpToTBPaCC, {
            updateTabsetPanel(session, "inTabset",
                              selected = "TBPaCC")
            
            shinyjs::disable("tb_reRun")
        })              
        
        ### Always start at 2024 
        baseYear <- 2024
        
        TBreactiveVar <- reactiveValues(cdrReact = NULL, oopReact = NULL, 
                                        medIncReact = NULL, yearReact = NULL, 
                                        calcCosts = NULL, costPlot = NULL, 
                                        costPlotStckP = NULL, costPlotStckC = NULL, 
                                        costCompPlot = NULL)
        
        #### Observe the country selected 
        observeEvent(input$tb_country, {
            
            covariateDF <- lookup_inputs(selectedCntry = input$tb_country, 
                                         selectedYear = baseYear, 
                                         infection = "TB")
            
            TBreactiveVar$medIncReact <- round(as.numeric(covariateDF[,"median.income.hh"]))
            TBreactiveVar$oopReact <- as.numeric(covariateDF[,"oop.perc.che"])
            TBreactiveVar$cdrReact <- as.numeric(covariateDF[,"cdr.rate"])
            TBreactiveVar$yearReact <- baseYear
            
            updateNumericInput(session,
                               inputId = "tb_medInc",
                               value = round(as.numeric(covariateDF[,"median.income.hh"])))
            
            shinyjs::disable("tb_medInc")
            
            updateNumericInput(session,
                               inputId = "tb_oopPerc",
                               value = as.numeric(covariateDF[,"oop.perc.che"]))
            
            shinyjs::disable("tb_oopPerc")
            
            updateNumericInput(session,
                               inputId = "tb_cdr",
                               value = as.numeric(covariateDF[,"cdr.rate"]))
            
            shinyjs::disable("tb_cdr")
            
            
            updateSelectInput(
                session,
                inputId = "year",
                selected = baseYear
            )
            
            TBreactiveVar$calcCosts <- tb_calculateCosts(selectedCntry = input$tb_country,
                                                         yearInput = TBreactiveVar$yearReact,
                                                         medIncInput = TBreactiveVar$medIncReact,
                                                         oopInput = TBreactiveVar$oopReact,
                                                         cdrInput = TBreactiveVar$cdrReact)
        })
        
        observeEvent(input$tb_year, {
            TBreactiveVar$yearReact <- as.numeric(input$tb_year) 
            
            covariateDF <- lookup_inputs(selectedCntry = input$tb_country, 
                                         selectedYear = as.numeric(input$tb_year), 
                                         infection = "TB")
            
            TBreactiveVar$medIncReact <- round(as.numeric(covariateDF[,"median.income.hh"]))
            TBreactiveVar$oopReact <- as.numeric(covariateDF[,"oop.perc.che"])
            TBreactiveVar$cdrReact <- as.numeric(covariateDF[,"cdr.rate"])
            
            updateNumericInput(session,
                               inputId = "tb_medInc",
                               value = round(as.numeric(covariateDF[,"median.income.hh"])))
            
            updateNumericInput(session,
                               inputId = "tb_oopPerc",
                               value = as.numeric(covariateDF[,"oop.perc.che"]))
            
            updateNumericInput(session,
                               inputId = "tb_cdr",
                               value = as.numeric(covariateDF[,"cdr.rate"]))

            
            TBreactiveVar$calcCosts <- tb_calculateCosts(selectedCntry = input$tb_country,
                                                         yearInput = TBreactiveVar$yearReact,
                                                         medIncInput = TBreactiveVar$medIncReact,
                                                         oopInput = TBreactiveVar$oopReact,
                                                         cdrInput = TBreactiveVar$cdrReact)
        })
        
        ##### If the user clicks the customize tb_covariates button, unlock the covariate boxes 
        observeEvent(input$tb_custCovs, { 
            
            shinyjs::enable("tb_medInc")
            shinyjs::enable("tb_oopPerc")
            shinyjs::enable("tb_cdr")
            
            shinyjs::enable("tb_reRun")
            
        })
        
        ##### Once the tb_covariates are updated and the user clicks calculate costs, 
        ##### relock the input boxes. 
        
        observeEvent(input$tb_reRun, {
            
            shinyjs::disable("tb_medInc")
            shinyjs::disable("tb_oopPerc")
            shinyjs::disable("tb_cdr")
            
            TBreactiveVar$medIncReact <- input$tb_medInc
            TBreactiveVar$oopReact <- input$tb_oopPerc
            TBreactiveVar$cdrReact <- input$tb_cdr
            
            TBreactiveVar$calcCosts <- tb_calculateCosts(selectedCntry = input$tb_country,
                                                         yearInput = TBreactiveVar$yearReact,
                                                         medIncInput = TBreactiveVar$medIncReact,
                                                         oopInput = TBreactiveVar$oopReact,
                                                         cdrInput = TBreactiveVar$cdrReact)
            
            shinyjs::disable("tb_reRun")
        })
        
        #### Plot country specific costs as bar plots
        output$tb_costPstack <- renderPlot({
            TBreactiveVar$costPlotStckP <- plot_stackP_costs(input$tb_country,
                                                             TBreactiveVar$calcCosts,
                                                             infection = "TB"
            )
            TBreactiveVar$costPlotStckP
        })

        output$tb_costCstack <- renderPlot({
            TBreactiveVar$costPlotStckC <- plot_stackC_costs(input$tb_country,
                                                             TBreactiveVar$calcCosts
            )
            TBreactiveVar$costPlotStckC
        })
        #
        output$tb_costCompPlot <- renderPlot({
            TBreactiveVar$costCompPlot <- plot_comp_costs(input$tb_country,
                                                          TBreactiveVar$calcCosts,
                                                          "TB"
            )
            TBreactiveVar$costCompPlot
        })
        
        ### Present country specific costs as table
        output$tb_costTable <- renderTable({
            TBreactiveVar$calcCosts %>%
                select(! c(LB,UB,Phase,Category))
        }, digits = 0)
        
        #### Present country specific studies as a table
        output$tb_studyTable <- renderTable({
            lookup_studies(input$tb_country)
        })
        
        output$tb_downloadFigure <- downloadHandler(
            filename = function() {
                # Use the selected dataset as the suggested file name
                paste0("PaCC_estimate_",  input$tb_country, ".pdf")
            },
            content = function(file) {
                # Write the dataset to the `file` that will be downloaded
                pList <- list(TBreactiveVar$costPlotStckP,
                              TBreactiveVar$costPlotStckC,
                              TBreactiveVar$costCompPlot)
                
                pdf(file = file, onefile = TRUE, width=11, height=8.5)
                for(i in 1:length(pList)){
                    print(pList[[i]])
                }
                dev.off()
            }
        )
        
        output$tb_downloadTable <- downloadHandler(
            filename = function() {
                # Use the selected dataset as the suggested file name
                paste0("PaCC_estimate_", input$tb_country, ".csv")
            },
            content = function(file) {
                # Write the dataset to the `file` that will be downloaded
                write.csv(TBreactiveVar$calcCosts, file)
            }
        )
    
}

###############################################################################
##### Run everything 
###############################################################################

shinyApp(ui = ui, server = server)
