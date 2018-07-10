## Import libraries
library(knitr)
library(readr)
library(rjags)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

## Import data
salaries <- read_csv("~/Desktop/salaries.csv", col_types = c(col_double(), col_character()))
salaries.df <- data.frame(salaries)

## Set the seed
set.seed(55)

## Initialize variables
salary <- salaries.df$Salary
dept <- salaries.df$Dept
deptNumeric <- as.numeric(as.factor(dept))
deptNames <- unique(dept) # department names
nFac <- length(salary) # number of faculty members in survey
D <- length(unique(dept)) # number of departments

## Create objects for JAGS
dataList <- list("salary" = salary,
                 "nFac" = nFac,
                 "D" = D,
                 "dept" = deptNumeric)

## Initialize unknown parameters
parameters <- c("theta", 
                "sigma2",
                "mu",
                "tau2")

## Set initial values for unknown parameters
initsValues <- list("theta" = rep(100, D), 
                    "prec_sigma2" = 1/100,
                    "mu" = 100,
                    "prec_tau2" = 1/100)

## Set the number of iteration for "tuning" 
adaptSteps <- 5000

## Set the number of iterations for "burn-in" 
burnInSteps <- 5000

## Set the number of chains to run
nChains <- 2

## Set the total number of iterations to save
numSavedSteps <- 5000

## Thinning to keep every iteration
thinSteps <- 1

## Iterations per chain
ITER <- ceiling((numSavedSteps*thinSteps)/nChains)
```

The starting values are listed in the code above. We should set both the number of adaptive steps and burn-in steps to 5,000, and we should run two chains an additional 2,500 iterations each.

```{r jags}
# Create, initialize, and adapt the JAGS model
jagsModel <- jags.model("~/Desktop/salariesmodel.txt", 
                        data = dataList, 
                        inits = initsValues, 
                        n.chains = nChains, 
                        n.adapt = adaptSteps)

## Burn-in the algorithm
update(jagsModel, n.iter = burnInSteps)

## Run the algorithm to get interations for inference
codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameters, 
                            n.iter = ITER, 
                            thin = thinSteps)
                            
## Make a dataframe with the posterior samples
mcmcChainDF <- data.frame(as.matrix(codaSamples, iters = T, chains = T))

## Create a vector with the variable names
varNames <- names(mcmcChainDF)[3:(dim(mcmcChainDF)[2])]

## Number of variables
nVars <- length(varNames)
mcmcChainDF$CHAIN <- as.factor(mcmcChainDF$CHAIN)

## Construct trace plots
p <- list()
for( k in 1:nVars ) {
  plot_frame <- mcmcChainDF
  plot_frame$dep_var <- mcmcChainDF[,varNames[k]]
  p[[k]] <- ggplot(plot_frame, aes(x = ITER, y = dep_var)) +
    geom_line(aes(color = CHAIN)) + 
    labs(y = varNames[k])
}

## Trace plots
do.call(grid.arrange, c(p, list("ncol" = 1)))

## Summary of the posterior distribution
postDFreshape <- melt(mcmcChainDF, id.vars = "ITER", measure.vars = c("theta.1.", "theta.2.", "theta.3.",
                                                                      "theta.4.", "theta.5.", "theta.6.",
                                                                      "theta.7.", "theta.8.", "mu"))

ggplot(postDFreshape, aes(x = variable, y = value)) +
  scale_x_discrete(labels = c(deptNames, "mu")) +
  labs(x = "Theta", y = "Posterior") +
  geom_boxplot()
  
## Calculate the credible interval
mcmcChainDF$varProp <- mcmcChainDF$sigma2 / (mcmcChainDF$sigma2+mcmcChainDF$tau2)

## Round the credible interval
round(quantile(mcmcChainDF$varProp, probs = c(.025, .975)), 2)

## Calculate the posterior probability
mcmcChainDF$predBM <- rnorm(nChains*ITER,
                            mcmcChainDF$theta.1.,     
                            mcmcChainDF$sigma)

## Posterior probability greater than 140
more <- mean(mcmcChainDF$predBM > 140)
more