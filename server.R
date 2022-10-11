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

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    # hybridC ----
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
                alph=input$alpha # alpha level per effect
            )
            
            if( input$hybrid_whichC== 'rC'){
                #hybrid.rC
                hybridrC_reslist=fun.hybridrC_sr(datpreplist = datprep_list,
                                                 K=input$hybridnumK
                                                 ,Min_rC= as.numeric(unlist(strsplit(input$MinsrC,","))) # c(0.3,0.3,0.3,0.3),
                                                 ,Max_rC= as.numeric(unlist(strsplit(input$MaxsrC,","))) #c(0.5,0.5,0.5,0.5)
                )
                hybridrC_summ=fun.summhybridC( hybridreslist = hybridrC_reslist )
                hybridC_sensitab = hybridrC_summ$abmedsensi
            }
            if( input$hybrid_whichC=='muC'){
                #hybrid.muC
                hybridmuC_reslist=fun.hybridmuC_sr(datpreplist = datprep_list
                                                   , Mean_C= as.numeric(unlist(strsplit(input$MeansC,","))) # c(0.39, 0.39, -0.39, -0.39)
                                                   , Sigma_C=diag( as.numeric(unlist(strsplit(input$VarsC,","))), nrow = 4)
                                                   ,K=input$hybridnumK
                )
                hybridmuC_summ=fun.summhybridC( hybridreslist = hybridmuC_reslist )
                hybridC_sensitab = hybridmuC_summ$abmedsensi
            }
        }
        
        hybridC_sensitab
    })
    
    output$hybridCsensi <- renderTable({
        hybridCsensitab()
    })
    
    
    # fixC ----
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
                alph=input$alpha # alpha level per effect
            )
            
            out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                   rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
                                                   )
            outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
            
            outplot=outsensi_ResPlot$sensiplots_IyonIm;
            outplot
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
                alph=input$alpha # alpha level per effect
            )
            out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                   rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
                                                   )
            outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
            
            outplot=outsensi_ResPlot$sensiplots_SyonIm
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
                alph=input$alpha # alpha level per effect
            )
            out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                   rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
                                                   )
            outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
            
            outsensi_ResPlot$sensiplots_IyonSm
            outplot=outsensi_ResPlot$sensiplots_IyonSm
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
                alph=input$alpha # alpha level per effect
            )
            out_sensifixC_JS_Z=fun.sensifixC_JS_Z( datpreplist = datprep_list,
                                                   rCvalues= as.numeric(unlist(strsplit(input$rCvalues,",")))
                                                   )
            outsensi_ResPlot=fun.plotsensifixC_JS_Z(outsensifixC_JS_Z = out_sensifixC_JS_Z )
            
            outsensi_ResPlot$sensiplots_SyonSm
            outplot=outsensi_ResPlot$sensiplots_SyonSm
        }
        if(!input$fixC){outplot=NULL}
        
        outplot
    })
    
    output$fixC.sensiplots_SyonSm <- renderPlot({
        sensiplots_SyonSm()
    })
    
})




