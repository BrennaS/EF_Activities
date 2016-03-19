if(!require(PEcAn.data.land)) devtools::install_github("PecanProject/pecan/modules/data.land")
require(rjags)
ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
## 1. Read tree data
trees <- read.csv("data/H2012AdultFieldData.csv")

## 2. Read tree ring data
rings <- Read_Tucson("data/TUCSON/")

## 3. merge inventory and tree ring data, extract most recent nyears
combined <- matchInventoryRings(trees,rings,nyears=15)

## take a look at the first few rows of data to see the structure
knitr::kable(combined[1:5,])

## 4. organize data into a list
data <- buildJAGSdata_InventoryRings(combined)

# y = increment (tree x year)
# z = dbh (tree x year)
# make sure to take a look at all the priors!
str(data)


n.iter = 500

## this code fuses forest inventory data with tree growth data (tree ring or dendrometer band)
## for the same plots. Code is a rewrite of Clark et al 2007 Ecol Appl into JAGS
TreeDataFusionMV = "
model{

### Loop over all individuals
for(i in 1:ni){

#### Data Model: DBH
for(t in 1:nt){
z[i,t] ~ dnorm(x[i,t],tau_dbh)
}

#### Data Model: growth
for(t in 2:nt){
inc[i,t] <- x[i,t]-x[i,t-1]
y[i,t] ~ dnorm(inc[i,t],tau_inc)
}

#### Process Model
for(t in 2:nt){
Dnew[i,t] <- x[i,t-1] + mu
x[i,t]~dnorm(Dnew[i,t],tau_add)
}

x[i,1] ~ dnorm(x_ic,tau_ic)
}  ## end loop over individuals

#### Priors
tau_dbh ~ dgamma(a_dbh,r_dbh)
tau_inc ~ dgamma(a_inc,r_inc)
tau_add ~ dgamma(a_add,r_add)
mu ~ dnorm(0.5,0.5)
}"

## state variable initial condition
z0 = t(apply(data$y,1,function(y){-rev(cumsum(rev(y)))})) + data$z[,ncol(data$z)] 

## JAGS initial conditions
nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(data$y,length(data$y),replace=TRUE)
  init[[i]] <- list(x = z0,tau_add=runif(1,1,5)/var(diff(y.samp),na.rm=TRUE),
                    tau_dbh=1,tau_inc=500,tau_ind=50,tau_yr=100,ind=rep(0,data$ni),year=rep(0,data$nt))
}

## compile JAGS model
j.model   <- jags.model (file = textConnection(TreeDataFusionMV),
                         data = data,
                         inits = init,
                         n.chains = 3)
## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_dbh","tau_inc","mu","tau_ind","tau_yr"),
                            n.iter = min(n.iter,2000))
plot(jags.out)

## run MCMC
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_dbh","tau_inc","mu",
                                               "tau_ind","tau_yr","ind","year"),
                            n.iter = n.iter)