#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#defaultW <- getOption("warn")
#options(warn = -1)


library(shiny)

library(mvtnorm)
library(ggplot2)
library(tidyverse)
options(tibble.print_max = Inf, tibble.width = Inf)
library(MplusAutomation)

source('fun_fixC_datprep_Z.R')
source('sourcefun_fixC_sr_Z.R')
source('sourcefun_plot_fixC_sr_Z.R')
# source('sensiFuns_fixChybridC.R')
source('sensiFuns_fixC_MCSA.R')
source('fun_BayesSensiC_Z.R')
source('fun_rC_aCbC.R')

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
    # transformation between rC (zero-order confounder correlations) and aCbC (confounder path coeffs), and check if they are within admissible ranges
  output$rC_aCbC_UI <-  renderUI({
    if ( !input$rC2aCbC ) return(NULL)
    if( input$rC2aCbC ) {
      div(id = "rC_2_aCbC",
          h4("Under the sensitivity analysis model, the set of specified confounder correlations corresponds to the
             following set of confounder path coefficients:",
             style = "display: inline-block;"),
          # a("Show/Hide",
          #    href = "javascript:toggleVisibility('rC_2_aCbC_TableSection');",
          #    class = "left-space"),
          div(id = "rC_2_aCbC_TableSection",
              tableOutput("rC_aCbC_tab") , style = "font-size:100%"
              # ,downloadButton("downloadSummary", "Download results")
          )
      )
    }

  })
  
  
  rC_aCbC_tabeventReactive <- eventReactive(input$update,{
     if(!input$rC2aCbC){
       rC2aCbC_tab <-NULL
     }
     if(input$rC2aCbC){
       datprep_list = fun.fixC.datprep_Z( 
         dat = read.csv(input$datcsv$datapath, header = T),
         Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
         SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
         Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
         SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
         Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
         Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
         Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
         Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
         within_person_residual_cov = input$withinperson_residual_cov,                 alph=input$alpha  # alpha level per effect
       )
       if ( length(mplusOut$results$errors)  > 0 ) {
         rC2aCbC_tab=data.frame( error_originalmodel=c('Mplus fitting original model M0 returns errors and/or warnings.',
                                                       unlist(mplusOut$results$errors), unlist(mplusOut$results$warnings) ) 
         )
       }
       if ( length(mplusOut$results$errors) ==0 ){
         rC2aCbC_tab = rC.to.aCbC( datpreplist = datprep_list , rC_ImSmIySy = as.numeric(unlist(strsplit(input$rC_ImSmIySy,",")))  )
         
       } 
     }
     
    rC2aCbC_tab
   })
  
   output$rC_aCbC_tab <- renderTable({
     rC_aCbC_tabeventReactive()
   })
   
   ## aCbC to rC ----
   output$aCbC_rC_UI <-  renderUI({
     if ( !input$aCbC2rC ) return(NULL)
     if( input$aCbC2rC ) {
       div(
           h4("Under the sensitivity analysis model, the set of specified confounder path coefficients corresponds to the
             following set of confounder correlations:",
              style = "display: inline-block;"),
           # a("Show/Hide",
           #    href = "javascript:toggleVisibility('aCbC_rC_tabsec');",
           #    class = "left-space"),
           div(id = "aCbC_rC_tabsec",
               tableOutput("aCbC_rC_tab") , style = "font-size:100%"
               # ,downloadButton("downloadSummary", "Download results")
           )
       )
     }
     
   })
   aCbC_rC_tab_eventReactive <- eventReactive(input$update,{
     if(!input$aCbC2rC){
       aCbC_rCtab <-NULL
     }
     if(input$aCbC2rC){
       datprep_list = fun.fixC.datprep_Z( 
         dat = read.csv(input$datcsv$datapath, header = T),
         Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
         SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
         Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
         SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
         Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
         Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
         Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
         Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
         within_person_residual_cov = input$withinperson_residual_cov,                 alph=input$alpha  # alpha level per effect
       )
       if ( length(mplusOut$results$errors)  > 0 ) {
         aCbC_rCtab=data.frame( error_originalmodel=c('Mplus fitting original model M0 returns errors and/or warnings.',
                                                       unlist(mplusOut$results$errors), unlist(mplusOut$results$warnings) ) 
         )
       }
       if ( length(mplusOut$results$errors)  ==0 ) {
         aCbC_rCtab = aCbC.to.rC( datpreplist = datprep_list, aCbC = as.numeric(unlist(strsplit(input$aCbC,","))) )
         
       }
     }
     
     aCbC_rCtab
     
   })
   
   output$aCbC_rC_tab <- renderTable({ aCbC_rC_tab_eventReactive() })
   
    # MCSA ----
   output$MCSA_UI <- renderUI({
     if(!input$hybridC){  NULL }
     if( input$hybridC ){
        div( h4('Monte Carlo sensitivity analysis',style = "display: inline-block;" ), 
             helpText('In the "param" column, the first six parameters are the path coefficients for the mediator-outcome relations and input-mediator relations (e.g., "SY.ON_IM" denotes the path from the mediator intercept to outcome slope);
                        the last four parameters are the mediation effects (e.g., "X_IM_IY" denotes the mediation effect of the input variable on the outcome intercept through the mediator intercept)'),
         tableOutput('hybridCsensi')
         
         )
       
     }
   })
    hybridCsensitab <- eventReactive(input$update, {
        if(!input$hybridC){
            hybridC_sensitab=NULL
        }
        if(input$hybridC){
            dat=read.csv(input$datcsv$datapath, header = T)
            # Ycols = c(1:5)  # columns for repeated measures of outcome from the first to the last time points 
            # SY.loadings = c(-4,-3,-2,-1,0)  # loadings of the outcome slope
            # Mcols = c(6:10)  # columns for repeated measures of mediator from the first to the last time points
            # SM.loadings = c(0,1,2,3,4) # loadings of the mediator slope
            # Xcol = 11 # column for the independent variable 
            # Z1cols = c(12) # columns for the level-2 covariates in both the level-2 models of mediator and outcome  
            # Z2cols = NULL # columns for the level-2 covariates in the level-2 model of mediator only
            # Z3cols = NULL # columns for the level-2 covariates in the level-2 model of outcome only
            # alph=0.05
            
            datprep_list = fun.fixC.datprep_Z( 
                dat = read.csv(input$datcsv$datapath, header = T),
                Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
                SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
                Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
                SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
                Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
                Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
                Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
                Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
                within_person_residual_cov = input$withinperson_residual_cov,                 alph=input$alpha  # alpha level per effect
            )
            if ( length(mplusOut$results$errors)  > 0 ) {
              hybridC_sensitab=data.frame( error_originalmodel=c('Mplus fitting original model M0 returns errors and/or warnings.',
                                                           unlist(mplusOut$results$errors), unlist(mplusOut$results$warnings) ) 
              )
            }
            if( length(mplusOut$results$errors)  == 0 ) {
              if( input$hybrid_whichC== 'rC'){
                #hybrid.rC
                hybridrC_reslist=fun.hybridrC_sr(datpreplist = datprep_list,
                                                 K=input$hybridnumK
                                                 ,Min_rC= as.numeric(unlist(strsplit(input$MinsrC,","))) # c(0.3,0.3,0.3,0.3),
                                                 ,Max_rC= as.numeric(unlist(strsplit(input$MaxsrC,","))) #c(0.5,0.5,0.5,0.5)
                )
                hybridrC_summ=fun.summhybridC( hybridreslist = hybridrC_reslist )
                hybridC_sensitab = hybridrC_summ$MCSAout
              }
              if( input$hybrid_whichC=='muC'){
                #hybrid.muC
                hybridmuC_reslist=fun.hybridmuC_sr(datpreplist = datprep_list
                                                   , Mean_C= as.numeric(unlist(strsplit(input$MeansC,","))) # c(0.39, 0.39, -0.39, -0.39)
                                                   , Sigma_C=diag( as.numeric(unlist(strsplit(input$VarsC,","))), nrow = 4)
                                                   ,K=input$hybridnumK
                )
                hybridmuC_summ=fun.summhybridC( hybridreslist = hybridmuC_reslist )
                hybridC_sensitab = hybridmuC_summ$MCSAout
              }
            }
            
        }
        
        hybridC_sensitab
    })
    
    output$hybridCsensi <- renderTable({
        hybridCsensitab()
    })
    
    
    # BSA ----
    output$BSA_UI <- renderUI({
      if(!input$BayesC){
        NULL
      }
      if(input$BayesC){
        div( h4('Bayesian sensitivity analysis'), 
             helpText('In the "param" column, the first six parameters are the path coefficients for the mediator-outcome relations and input-mediator relations (e.g., "SY.ON_IM" denotes the path from the mediator intercept to outcome slope);
                        the last four parameters are the mediation effects (e.g., "X_IM_IY" denotes the mediation effect of the input variable on the outcome intercept through the mediator intercept)'),
             tableOutput('BayesCsensi')
        )
      }
    })
    
    
    BayesCsensitab <- eventReactive(input$update, {
      if(!input$BayesC){
        BayesC_sensitab=NULL
      }
      if(input$BayesC){
        dat=read.csv(input$datcsv$datapath, header = T)
        
        BayesC_reslist = fun.BayesC_Z(
          dat = read.csv(input$datcsv$datapath, header = T),
          Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
          SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
          Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
          SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
          Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
          Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
          Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
          Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
          within_person_residual_cov = input$withinperson_residual_cov,
          
          priors_muC=input$priors_muC,
          BITERATIONS=input$BITERATIONS,THIN=input$THIN,CHAINS=input$CHAINS,
          
          non_default_priors_original = input$non_default_priors_original,
          Bayes_originalmodel= input$priors_originalmodel,
          priors_originalmodel=input$priors_originalmodel
        
        )
        BayesC_sensitab = BayesC_reslist$BSAout
      } 
      data.frame(BayesC_sensitab)
      
      } )
    
    output$BayesCsensi <- renderTable({
      BayesCsensitab()
    })
    
    # FSA ----
    output$FSA_UI <- renderUI({
      if(!input$fixC){
        NULL
      }
      if(input$fixC){
        div(   h4('Frequentist sensitivity analysis'),
               # textOutput('original_errors'),
          plotOutput("fixC.sensiplots_IyonIm"),
          plotOutput("fixC.sensiplots_SyonIm"),
          plotOutput("fixC.sensiplots_IyonSm"),
          plotOutput("fixC.sensiplots_SyonSm")
        )
      }
    })
    
    # original_errorseventReactive <- eventReactive(input$update,{
    #   datprep_list = fun.fixC.datprep_Z( 
    #     dat = read.csv(input$datcsv$datapath, header = T),
    #     Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
    #     SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
    #     Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
    #     SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
    #     Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
    #     Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
    #     Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
    #     Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
    #     within_person_residual_cov = input$withinperson_residual_cov,
    #     alph=input$alpha  # alpha level per effect
    #   )
    #   if(length(datprep_list$mplusOut$results$errors) >0 ){
    #     print(c('Mplus fitting original model M0 returns errors and/or warnings.',
    #             unlist(mplusOut$results$errors), unlist(mplusOut$results$warnings)) )
    #   }
    #   if( length(datprep_list$mplusOut$results$errors)==0 ){
    #     NULL 
    #   }
    #   
    # })
    # output$original_errors <-renderText( {original_errorseventReactive() })
    
    ## IyonIm
    sensiplots_IyonIm <- eventReactive(input$update,{ 
        if(input$fixC){
            dat=read.csv(input$datcsv$datapath, header = T)
            # datprep_list = fun.fixC.datprep_Z( 
            #     dat = read.csv(input$datcsv$datapath, header = T),
            #     Ycols = c(1:5), # columns for repeated measures of outcome from the first to the last time points
            #     SY.loadings = c(-4,-3,-2,-1,0), # loadings of the outcome slope
            #     Mcols = c(6:10), # columns for repeated measures of mediator from the first to the last time points
            #     SM.loadings = c(0,1,2,3,4), # loadings of the mediator slope
            #     Xcol = 11, # column for the independent variable
            #     Z1cols = c(12), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
            #     Z2cols = NULL, # columns for the level-2 covariates in the level-2 model of mediator only
            #     Z3cols = NULL, # columns for the level-2 covariates in the level-2 model of outcome only
            #     alph=0.05 # alpha level per effect
            # )
            datprep_list = fun.fixC.datprep_Z( 
                dat = read.csv(input$datcsv$datapath, header = T),
                Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
                SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
                Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
                SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
                Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
                Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
                Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
                Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
                within_person_residual_cov = input$withinperson_residual_cov,
                         alph=input$alpha  # alpha level per effect
            )
            
             if( length(datprep_list$mplusOut$results$errors)==0 ){
              out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                     rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
              )
              outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
              
              outplot=outsensi_ResPlot$sensiplots_IyonIm;
             }
             if( length(datprep_list$mplusOut$results$errors) > 0 ){
               outplot=NULL
             }
            
        }
        if(!input$fixC){
            outplot=NULL
        }
        outplot
    })
    output$fixC.sensiplots_IyonIm <- renderPlot({
        sensiplots_IyonIm()
    })
    
    ## SyonIm
    sensiplots_SyonIm <- eventReactive(input$update,{ 
        if(input$fixC){
            dat=read.csv(input$datcsv$datapath, header = T)
            # datprep_list = fun.fixC.datprep_Z( 
            #     dat = read.csv(input$datcsv$datapath, header = T),
            #     Ycols = c(1:5), # columns for repeated measures of outcome from the first to the last time points
            #     SY.loadings = c(-4,-3,-2,-1,0), # loadings of the outcome slope
            #     Mcols = c(6:10), # columns for repeated measures of mediator from the first to the last time points
            #     SM.loadings = c(0,1,2,3,4), # loadings of the mediator slope
            #     Xcol = 11, # column for the independent variable
            #     Z1cols = c(12), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
            #     Z2cols = NULL, # columns for the level-2 covariates in the level-2 model of mediator only
            #     Z3cols = NULL, # columns for the level-2 covariates in the level-2 model of outcome only
            #     alph=0.05 # alpha level per effect
            # )
            datprep_list = fun.fixC.datprep_Z( 
                dat = read.csv(input$datcsv$datapath, header = T),
                Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
                SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
                Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
                SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
                Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
                Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
                Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
                Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
                within_person_residual_cov = input$withinperson_residual_cov,                 alph=input$alpha  # alpha level per effect
            )
            if( length(datprep_list$mplusOut$results$errors) == 0 ){
              out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                     rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
              )
              outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
              
              outplot=outsensi_ResPlot$sensiplots_SyonIm
            }
            if( length(datprep_list$mplusOut$results$errors) > 0 ){
              outplot=NULL
            }
            
        }
        if(!input$fixC){outplot=NULL}
        
        outplot
    })
    output$fixC.sensiplots_SyonIm <- renderPlot({
        sensiplots_SyonIm()
    })
    
    ## IyonSm
    sensiplots_IyonSm <- eventReactive(input$update,{ 
        if(input$fixC){
            dat=read.csv(input$datcsv$datapath, header = T)
            
            # datprep_list = fun.fixC.datprep_Z( 
            #     dat = read.csv(input$datcsv$datapath, header = T),
            #     Ycols = c(1:5), # columns for repeated measures of outcome from the first to the last time points
            #     SY.loadings = c(-4,-3,-2,-1,0), # loadings of the outcome slope
            #     Mcols = c(6:10), # columns for repeated measures of mediator from the first to the last time points
            #     SM.loadings = c(0,1,2,3,4), # loadings of the mediator slope
            #     Xcol = 11, # column for the independent variable
            #     Z1cols = c(12), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
            #     Z2cols = NULL, # columns for the level-2 covariates in the level-2 model of mediator only
            #     Z3cols = NULL, # columns for the level-2 covariates in the level-2 model of outcome only
            #     alph=0.05 # alpha level per effect
            # )
            datprep_list = fun.fixC.datprep_Z( 
                dat = read.csv(input$datcsv$datapath, header = T),
                Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
                SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
                Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
                SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
                Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
                Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
                Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
                Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
                within_person_residual_cov = input$withinperson_residual_cov,                 alph=input$alpha  # alpha level per effect
            )
            if( length(datprep_list$mplusOut$results$errors) > 0 ){
              outplot=NULL
            }
            if( length(datprep_list$mplusOut$results$errors) == 0 ){
              out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                     rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
              )
              outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
              
              outsensi_ResPlot$sensiplots_IyonSm
              outplot=outsensi_ResPlot$sensiplots_IyonSm
            }
            
        }
        if(!input$fixC){outplot=NULL}
        
        outplot
    })
    output$fixC.sensiplots_IyonSm <- renderPlot({
        sensiplots_IyonSm()
    })
    
    ## SyonSm
    sensiplots_SyonSm <- eventReactive(input$update,{ 
        if(input$fixC){
            dat=read.csv(input$datcsv$datapath, header = T)
            
            # datprep_list = fun.fixC.datprep_Z( 
            #     dat = read.csv(input$datcsv$datapath, header = T),
            #     Ycols = c(1:5), # columns for repeated measures of outcome from the first to the last time points
            #     SY.loadings = c(-4,-3,-2,-1,0), # loadings of the outcome slope
            #     Mcols = c(6:10), # columns for repeated measures of mediator from the first to the last time points
            #     SM.loadings = c(0,1,2,3,4), # loadings of the mediator slope
            #     Xcol = 11, # column for the independent variable
            #     Z1cols = c(12), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
            #     Z2cols = NULL, # columns for the level-2 covariates in the level-2 model of mediator only
            #     Z3cols = NULL, # columns for the level-2 covariates in the level-2 model of outcome only
            #     alph=0.05 # alpha level per effect
            # )
            datprep_list = fun.fixC.datprep_Z( 
                dat = read.csv(input$datcsv$datapath, header = T),
                Ycols = as.numeric(unlist(strsplit(input$Ycolumns,","))), # columns for repeated measures of outcome from the first to the last time points
                SY.loadings =  as.numeric(unlist(strsplit(input$Yslope_loadings,","))) , # loadings of the outcome slope
                Mcols = as.numeric(unlist(strsplit(input$Mcolumns,","))), # columns for repeated measures of mediator from the first to the last time points
                SM.loadings =  as.numeric(unlist(strsplit(input$Mslope_loadings,","))), # loadings of the mediator slope
                Xcol = as.numeric(unlist(strsplit(input$Xcolumn,","))), # column for the independent variable
                Z1cols = as.numeric(unlist(strsplit(input$covariates.MYboth,","))), # columns for the level-2 covariates in both the level-2 models of mediator and outcome
                Z2cols = as.numeric(unlist(strsplit(input$covariates.Monly,","))), # columns for the level-2 covariates in the level-2 model of mediator only
                Z3cols =  as.numeric(unlist(strsplit(input$covariates.Yonly,","))), # columns for the level-2 covariates in the level-2 model of outcome only
                within_person_residual_cov = input$withinperson_residual_cov,                 alph=input$alpha  # alpha level per effect
            )
            if( length(datprep_list$mplusOut$results$errors) > 0 ){
              outplot=NULL
            }
            if( length(datprep_list$mplusOut$results$errors) == 0 ){
              out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                     rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
              )
              outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
              
              outsensi_ResPlot$sensiplots_SyonSm
              outplot=outsensi_ResPlot$sensiplots_SyonSm
            }
            
        }
        if(!input$fixC){outplot=NULL}
        
        outplot
    })
    
    output$fixC.sensiplots_SyonSm <- renderPlot({
        sensiplots_SyonSm()
    })
    
})




