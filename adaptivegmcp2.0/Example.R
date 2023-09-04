#Loading Functions
library(usethis)
library(devtools)
load_all()
#####################
##  Design Inputs  ##

#Number of looks
K <- 2

#Number of hypothesis
D <- 4

#Type-1 Error
alpha <- 0.05

#Initial Weights

WI <- c(1/4,1/4,1/4,1/4)

#Transition Matrix

G <- matrix(c(0,1/2,1/2,0,

              1/2,0,0,1/2,

              0,1,0,0,

              1,0,0,0),

            nrow = 4, byrow = T)

#Correlation Matrix: Non-Parametric
Correlation <- matrix(c(1,0.5,NA,NA,

                 0.5,1,NA,NA,

                 NA,NA,1,NA,

                 NA,NA,NA,1),

               nrow = 4, byrow = T)

Threshold = c(0.025,0.05)

gMCPLite::hGraph(nHypotheses = D, nameHypotheses = c('H1','H2','H3','H4'),
       alphaHypotheses = alpha*WI, m = G)


adaptGMCP_PC(K=K,D=D,G=G,WI = WI,
             Correlation = Correlation,Selection = T,
             Threshold =Threshold)
#############################################################


#Number of looks
K <- 2

#Number of hypothesis
D <- 4

#Type-1 Error
alpha <- 0.025

#Weights : Bonferroni equal weights for all primary hypothesis
WI <- c(1/4,1/4,1/4,1/4)

#Transition Matrix: Bonferroni Holms Procedure
G <- matrix(c(0,1/3,1/3,1/3,
              1/3,0,1/3,1/3,
              1/3,1/3,0,1/3,
              1/3,1/3,1/3,0),
            nrow = 4,byrow = T)

#Correlation Matrix: Non-Parametric
Correlation <- matrix(c(1,NA,NA,NA,
                        NA,1,NA,NA,
                        NA,NA,1,NA,
                        NA,NA,NA,1),
                      nrow = 4,byrow = T)

gMCPLite::hGraph(nHypotheses = D, nameHypotheses = c('H1','H2','H3','H4'),
                 alphaHypotheses = alpha*WI, m = G)

Threshold = c(0.0000001,0.025)

adaptGMCP_PC(K=K,D=D,G=G,WI = WI,
             Correlation = Correlation,Selection = T,
             Threshold =Threshold)
















