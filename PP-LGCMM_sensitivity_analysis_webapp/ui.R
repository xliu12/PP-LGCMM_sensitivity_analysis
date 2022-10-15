#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinybusy)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # Application title
    titlePanel("Sensitivity analysis of mediation inference from the parallel latent growth curve mediation model to the omission of a potential between-person pretreatment confounder of the mediator-outcome relations"), 
    
    add_busy_spinner(spin = "fading-circle"),
    
    sidebarLayout(
        sidebarPanel(
            fileInput("datcsv",
                      "csv file (include header) for the data frame containing variables in the model",
                      accept = c("text")
            ),
            helpText("See an example dataset with 5 time points (and a baseline covariate Z included in both the latent mediator models and latent outcome models)", a("exampledata", href='https://osf.io/cv25h/')),
            helpText("The illustrative values specified below are for this example datset."),
            textAreaInput('Ycolumns', 
                      'Column numbers for repeated measures of outcome from the first to the last time points in the input datcsv file (the column numbers should be separated by comma.)', 
                      value = "1,2,3,4,5"
            ),
            textAreaInput('Yslope_loadings', 
                      'Outcome slope loadings on repeated measures of the outcome from the first to the last time points (the loading values should be separated by comma.)', 
                      value = "-4,-3,-2,-1,0"
            ),
            textAreaInput('Mcolumns', 
                      'Column numbers for repeated measures of mediator from the first to the last time points in the input datcsv file (the column numbers should be separated by comma.)', 
                      value = "6,7,8,9,10"
            ),
            textAreaInput('Mslope_loadings', 
                      'Mediator slope loadings on repeated measures of the mediator from the first to the last time points (the loading values should be separated by comma).', 
                      value = "0,1,2,3,4"
            ),
            textAreaInput('Xcolumn', 
                      'Column number corresponding to the treatment variable (time-invariant)  in the input datcsv file', 
                      value = "11"
            ),
            textAreaInput('covariates.MYboth', 
                      'Column numbers corresponding to the covariates (time-invariant) included in both the latent mediator models and latent outcome models  in the input datcsv file (separated by comma).', 
                      value = "12"
            ),
            textAreaInput('covariates.Monly', 
                      'Columns in the input datcsv file corresponding to the covariates (time-invariant) included in only the latent mediator models; the column numbers should be separated by comma.', 
                      value = ""
            ),
            textAreaInput('covariates.Yonly', 
                      'Columns in the input datcsv file corresponding to the covariates (time-invariant) included in only the latent outcome models; the column numbers should be separated by comma.', 
                      value = ""
            ),
            textAreaInput('withinperson_residual_cov',
                      'Within-person residual covariance structure specification',
                      value = "eM1-eM5 (sigma_eM);
 eY1-eY5 (sigma_eY);"),
            helpText("The specification will be used to fit the original model in Mplus.
                     Therefore, the specification here follows Mplus syntax.",br(),
                     "Denote the mediator's within-person residual at the t-th time point as: eMt. For example, let eM2 denote the mediator's within-person residual at the 2nd time point.", br(),
                     "Denote the outcome's within-person residual at the t-th time point as: eYt. For example, let eY2 denote the outcome's within-person residual at the 2nd time point.", br(),
                     "In the example input shown above, there are 5 time points, and the within-person residual covariance structure are the 
                     covariance structure for the vector (eM1, eM2, eM3, eM4, eM5, eY1, eY2, eY3, eY4, eY5)'. "),
            numericInput("alpha",
                         "significance level for testing a specific indirect effect",
                         min = 0,
                         value = .05)
            ##---FSA----
            ,checkboxInput("fixC",
                          "Frequentist sensitivity analysis",
                          value = FALSE
                          )
            ,conditionalPanel(
                condition="input.fixC==1",
                textAreaInput('rCvalues', 
                          'plausible values of the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma; these will be used as the panel names of sensitivity plots):', 
                          value = "0.2,0.2,0.2,0.2" )
            )
            
            ##---MCSA----
            ,checkboxInput("hybridC",
                           "Monte Carlo sensitivity analysis",
                           value = FALSE
            )
            ,conditionalPanel(
                condition="input.hybridC==1",
                numericInput("hybridnumK",
                             "how many sets of confounder sensitivity parameters to draw?",
                             value=1000),
                radioButtons("hybrid_whichC",
                            "Specify priors for which type of confounder sensitivity paramters?",
                            choices = c('Confounder correlations (uniform priors)'='rC'
                                        ,'Confounder path coefficients (normal priors)'='muC'
                            ),
                            selected = 'muC')
                ,conditionalPanel(
                    condition="input.hybrid_whichC=='rC'",
                    textAreaInput('MinsrC', 
                              'minimums of the uniform priors for the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.2,0.2,0.2,0.2" ),
                    textAreaInput('MaxsrC', 
                              'maximums of the uniform priors for the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.4,0.4,0.4,0.4" )
                )
                ,conditionalPanel(
                    condition="input.hybrid_whichC=='muC'",
                    textAreaInput('MeansC', 
                              'means of the nromal priors for the confounder path coefficients to the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.3,0.3,0.3,0.3" ),
                    textAreaInput('VarsC', 
                              'variances of the normal priors for the confounder path coefficients to the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.01,0.01,0.01,0.01" )
                )
            )
            
            ##---BSA----
            ,checkboxInput("BayesC",
                           "Bayesian sensitivity analysis",
                           value = FALSE
            )
            ,conditionalPanel(
              condition="input.BayesC==1",
              textAreaInput('priors_muC', 
                        'Informative prior distribution for the confounder path coefficients:', 
                        value = "
aIM_C ~ N(0.39, 0.25); aSM_C ~ N(0.39, 0.25);
bIY_C ~ N(-0.39, 0.25); bSY_C ~ N(-0.39, 0.25); "),
              helpText("The specified informative prior will be used to estimate the sensitivity analysis model with the Bayes estimator in Mplus.",br(),
                       "The Mplus syntax involving the potential confounder C is set as follows (following Harring et al., 2017)", br(),
                       " C by; [C@0]; C@1;
  IM on C (aIM_C); SM on C (aSM_C); IY on C (bIY_C);  SY on C (bSY_C);",br(),
                       "Please input informative priors for the four confounder path coefficients in the format of the MODELPRIORS section in Mplus." ) ,
              numericInput("CHAINS", "Number of chains for Bayesian posterior sampling", value = 2),
              numericInput("BITERATIONS", "Number of iterations for each chain", value = 20000),
              numericInput("THIN", "Thinning rate for each chain", value = 5),
              
              checkboxInput("non_default_priors_original",
                            "whether use non-default priors for the original model parameters", value = FALSE)   
             ,helpText("By default, priors for the original model parameters are the default priors in Mplus."),
             conditionalPanel(
               condition = "input.non_default_priors_original==1", 
               textAreaInput("Bayes_originalmodel",
                         "Input the syntax of the original model in the format of the MODEL section in Mplus",
                         value = "
IM by M1-M5@1;
SM by M1@0 M2@1 M3@2 M4@3 M5@4;
IY by Y1-Y5@1;
SY by Y1@-4 Y2@-3 Y3@-2 Y4@-1 Y5@0;
[M1-M5@0]; [Y1-Y5@0]; M1-M5; Y1-Y5;
[IM SM];[IY SY]; IY on X Z; SY on X Z;

IM on X  (aIM_X);
SM on X  (aSM_X);
IY on IM SM (bIY_IM bIY_SM);
SY on IM SM (bSY_IM bSY_SM);

IM (revarIM);SM (revarSM); IM with SM (IMwSM);
IY (revarIY);SY (revarSY); IY with SY (IYwSY);
IY SY with IM@0 SM@0;
  "
               ),
               helpText("For obtaining inferences for the mediation effects, which are product of the mediation path coefficients, ",br(),
               "please label the path coefficients of mediator intercept and slope on X as aIM_X and aSM_X, respectively; ",br(),
                        "label the path coefficients of outcome intercept on mediator intercept and slope as bIY_IM and bIY_SM, respectively; ",br(),
                        "label the path coefficients of outcome slope on mediator intercept and slope as bSY_IM and bSY_SM, respectively;"),
               textAreaInput("priors_originalmodel",
                         "Input the priors for the original model parameters in the format of the MODELPRIORS section in Mplus",
                         value = "revarIM ~ iw(1,3); revarSM ~iw(1,3);IMwSM ~ iw(0,3);revarIY ~ iw(1,3); revarSY ~ iw(1,3);IYwSY ~ iw(0,3);"
               )
             )
              
            ),
            
            
            ####----transform rC and aCbC and check admissibility
            checkboxInput("rC2aCbC",
                          "Transform a set of confounder correlations to the corresponding set of confounder path coefficients under the sensitivity analysis model, 
                          and check whether the set of confounder correlations (and confounder path coefficients) is in the admissible ranges", 
                          value = FALSE),
            conditionalPanel(condition = "input.rC2aCbC==1", 
                             textAreaInput('rC_ImSmIySy', 
                                       'Values of the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                                       value = "0.1,0.1,0.1,0.1" )
                             ),
            checkboxInput("aCbC2rC",
                          "Transform a set of confounder path coefficients to the corresponding set of confounder correlations under the sensitivity analysis model, 
                          and check whether the set of confounder path coefficients (and confounder correlations) is in the admissible ranges", 
                          value = FALSE),
            helpText("Admissible confounder correlations (confounder path coefficients) are those with which the model-implied covariance matrix of the potential confounder C, latent intercepts and slopes of mediator and outcome, treatment variable, and observed baseline covariates under the sensitivity analysis model positive definite"),
            conditionalPanel(condition = "input.aCbC2rC==1", 
                             textAreaInput('aCbC', 
                                       'Values of the confounder path coefficients from potential confounder C to the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                                       value = "0.14,0.14,0.14,0.14" )
            ),
            ## go----
            actionButton("update", "Go"),
            helpText('In the Frequentist sensitivity analysis and Monte Carlo sensitivity analysis, the original model is estimated with the maximum likelihood (ML) estimator in Mplus',br(),
                     'In the Bayesian sensitivity analysis, the original and sensitivity analysis models are estimated with the Bayes estimator in Mplus.', br(),
                     'Mplus is run in R via the package MplusAutomation (Hallquist & Wiley, 2018)' ),
            helpText('References: ', br(),
                     'Hallquist, M. N. & Wiley, J. F. (2018). MplusAutomation: An R Package for Facilitating Large-Scale Latent Variable Analyses in Mplus. Structural Equation Modeling, 25, 621-638. doi: 10.1080/10705511.2017.1402334.', br(),
                     'Harring, Jeffrey R., McNeish, Daniel M., and Hancock, Gregory R., Using phantom variables in structural equation modeling to assess model sensitivity to external misspecification. Psychological Methods 22, 4 (2017), pp. 616--631.'   
            ),
            tags$div(class="body", checked=NA,
                     tags$a(href="https://github.com/xliu12/PP-LGCMM_sensitivity_analysis", "See here for the R code and more details about the sensitivity analysis methods.")
            )
        ),
        
        # Show a plot of the generated plots
        mainPanel(
          h4('Diagrams of the original PP-LGCMM (Model M0) and sensitivity analysis PP-LGCMM (Model M1)'),
          img(src='fig_models0and1.png',height=400),
          helpText('The sensitivity analysis results provide information regarding how the original mediation inferences would be altered by a potential omitted between-person confounder of the mediator-outcome relations.
        Interpreting and gauging the robustness of original mediation inference should be based on substantive knowledge.'),
          
          #   h4('Monte Carlo sensitivity analysis'),
          #   tableOutput('hybridCsensi'),
          # 
          # h4('Bayesian sensitivity analysis'),
          # tableOutput('BayesCsensi'),
          # 
          # helpText('In the "param" column, the first six parameters are the path coefficients for the mediator-outcome relations and input-mediator relations (e.g., "SY.ON_IM" denotes the path from the mediator intercept to outcome slope);
          #            the last four parameters are the mediation effects (e.g., "X_IM_IY" denotes the mediation effect of the input variable on the outcome intercept through the mediator intercept)'),
          
            # h4('Frequentist sensitivity analysis'),
            # plotOutput("fixC.sensiplots_IyonIm"),
            # plotOutput("fixC.sensiplots_SyonIm"),
            # plotOutput("fixC.sensiplots_IyonSm"),
            # plotOutput("fixC.sensiplots_SyonSm"),
          
          uiOutput('BSA_UI'),
          uiOutput('MCSA_UI'),
          
          # h4('Transformation between confounder correlations and confounder path coeffs, and check if they are within admissible ranges') ,
           
          #tableOutput('rC_aCbC_tab'),
          uiOutput('FSA_UI'),
          uiOutput('aCbC_rC_UI'),
         uiOutput('rC_aCbC_UI')
        
        )
    )
))
