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
    titlePanel("Sensitivity analysis of the mediation inference from the parallel latent growth curve mediation model to the omission of a potential time-invariant confounder of the mediator-outcome relations"), 
    
    add_busy_spinner(spin = "fading-circle"),
    
    sidebarLayout(
        sidebarPanel(
            fileInput("datcsv",
                      "csv file (include header) for the data frame containing variables in the LGCM model only",
                      accept = c("text")
            ),
            helpText("See an example dataset with 5 time points (and a baseline covariate Z included in both the latent mediator models and latent outcome models)", a("exampledata", href='https://osf.io/cv25h/')),
            helpText("The illustrative values specified below are for this example datset."),
            textInput('Ycolumns', 
                      'Column numbers for repeated measures of outcome from the first to the last time points in the input datcsv file (the column numbers should be separated by comma.)', 
                      value = "1,2,3,4,5"
            ),
            textInput('Yslope_loadings', 
                      'Outcome slope loadings on repeated measures of the outcome from the first to the last time points (the loading values should be separated by comma.)', 
                      value = "-4,-3,-2,-1,0"
            ),
            textInput('Mcolumns', 
                      'Column numbers for repeated measures of mediator from the first to the last time points in the input datcsv file (the column numbers should be separated by comma.)', 
                      value = "6,7,8,9,10"
            ),
            textInput('Mslope_loadings', 
                      'Mediator slope loadings on repeated measures of the mediator from the first to the last time points (the loading values should be separated by comma).', 
                      value = "0,1,2,3,4"
            ),
            textInput('Xcolumn', 
                      'Column number corresponding to the independent variable (time-invariant)  in the input datcsv file', 
                      value = "11"
            ),
            textInput('covariates.MYboth', 
                      'Column numbers corresponding to the covariates (time-invariant) included in both the latent mediator models and latent outcome models  in the input datcsv file (separated by comma).', 
                      value = "12"
            ),
            textInput('covariates.Monly', 
                      'the columns in the input datcsv file corresponding to the covariates (time-invariant) included in only the latent mediator models; the column numbers should be separated by comma.', 
                      value = ""
            ),
            textInput('covariates.Yonly', 
                      'the columns in the input datcsv file corresponding to the covariates (time-invariant) included in only the latent outcome models; the column numbers should be separated by comma.', 
                      value = ""
            ),
            
            numericInput("alpha",
                         "significance level for testing a specific indirect effect",
                         min = 0,
                         value = .05)
            ,checkboxInput("fixC",
                          "Frequentist sensitivity analysis (where confounder sensitivity parameters are viewed as fixed values)",
                          value = FALSE
                          )
            ,conditionalPanel(
                condition="input.fixC==1",
                textInput('rCvalues', 
                          'plausible values of the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma; these will be used as the panel names of sensitivity plots):', 
                          value = "0.2,0.2,0.2,0.2" )
            )
            
            ,checkboxInput("hybridC",
                           "Monte Carlo sensitivity analysis (where multiple [K] sets of confounder sensitivity parameters are randomly drawn from user-specified prior distributions, and then used to make mediational inferences using frequentist maximum likelihood methods)",
                           value = TRUE
            )
            ,conditionalPanel(
                condition="input.hybridC==1",
                numericInput("hybridnumK",
                             "how many sets of confounder sensitivity parameters to draw? (the results will be summarized across all the draws)",
                             value=1000),
                radioButtons("hybrid_whichC",
                            "Specify priors for which type of confounder sensitivity paramters?",
                            choices = c('Confounder correlations (uniform priors)'='rC'
                                        ,'Confounder path coefficients (normal priors)'='muC'
                            ),
                            selected = 'muC')
                ,conditionalPanel(
                    condition="input.hybrid_whichC=='rC'",
                    textInput('MinsrC', 
                              'minimums of the uniform priors for the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.2,0.2,0.2,0.2" ),
                    textInput('MaxsrC', 
                              'maximums of the uniform priors for the confounder correlations with the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.4,0.4,0.4,0.4" )
                )
                ,conditionalPanel(
                    condition="input.hybrid_whichC=='muC'",
                    textInput('MeansC', 
                              'means of the nromal priors for the confounder path coefficients to the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.3,0.3,0.3,0.3" ),
                    textInput('VarsC', 
                              'variances of the normal priors for the confounder path coefficients to the latent mediator intercept, mediator slope, outcome intercept and outcome slope (separated by comma):', 
                              value = "0.01,0.01,0.01,0.01" )
                )
            )
            
            ,helpText('The original model is estimated with the maximum likelihood (ML) estimator in Mplus via the package MplusAutomation (Hallquist & Wiley, 2018)',br(),
                      'Hallquist, M. N. & Wiley, J. F. (2018). MplusAutomation: An R Package for Facilitating Large-Scale Latent Variable Analyses in Mplus. Structural Equation Modeling, 25, 621-638. doi: 10.1080/10705511.2017.1402334.')
            
            ,actionButton("update", "Go")
        ),
        
        # Show a plot of the generated plots
        mainPanel(
            h4('Monte Carlo sensitivity analysis'),
            tableOutput('hybridCsensi'),
            helpText('In the "param" column, the first six parameters are the path coefficients for the mediator-outcome relations and input-mediator relations (e.g., "SY.ON_IM" denotes the path from the mediator intercept to outcome slope);
                     the last six parameters are the specific and total indirect effects (e.g., "X_IM_IY" denotes the indirect effect of the input variable on the outcome intercept through the mediator intercept)'),
            
            h4('Frequentist sensitivity analysis'),
            plotOutput("fixC.sensiplots_IyonIm"),
            plotOutput("fixC.sensiplots_SyonIm"),
            plotOutput("fixC.sensiplots_IyonSm"),
            plotOutput("fixC.sensiplots_SyonSm"),
            
            h5('criterion of robustness should be based on substantive knowledge.')
        )
    )
))
