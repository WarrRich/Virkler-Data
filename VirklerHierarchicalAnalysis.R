library(R2jags); library(tidyverse); library(MASS); library(mvtnorm)

### Data Prep
data <- read.csv('VirklerData.csv')
names(data) = c('item','cycles','length')

crack.length <- data[,3]; items <- data[,1]; cycles <- data[,2]; n <- length(crack.length)
itemNum <- rep(0,68)
for (i in 1:68) {
  itemNum[i] <- length(which(data[,1]==i))
}
cumItemNum <- cumsum(c(0,itemNum[-68]))


### EDA
#pdf(file="VirklerDataPlot.pdf",width=5.5,height=3.667)
ggplot(data, aes(x = cycles, y = length, group = as.factor(items))) +
  geom_line(size=0.4,color = "gray") + 
  labs(x = "Thousands of Cycles", y = "1/2 Crack Length in mm") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(0, 250)) +
  geom_hline(yintercept=30)+
  geom_point()
#dev.off()

##################
# Metal Fatigue Crack Growth JAGS Model

crack.growth.model <- "model {
  D0 <- 9
  pi <- 4*atan(1)
  for (i in 1:68) {
    for (j in 1:itemNum[i]) {
      crack.length[cumItemNum[i]+j] ~ dnorm(mu[cumItemNum[i]+j],1/sigma.eps^2)
      stuff[cumItemNum[i]+j] <- D0^(1-theta[i,2]/2)+(1-theta[i,2]/2)*exp(theta[i,1])*(pi)^(theta[i,2]/2)*cycles[cumItemNum[i]+j]
      mu[cumItemNum[i]+j] <- ifelse(abs(theta[i,2]) < 1e-12, 
        D0*exp(exp(theta[i,1]*pi)*cycles[cumItemNum[i]+j]), 
        ifelse(stuff[cumItemNum[i]+j] < 0,-100,stuff[cumItemNum[i]+j]^(2/(2-theta[i,2])))
      )
    }
    theta[i,1:2] ~ dmnorm.vcov(mu.theta[1:2], sigma.theta[1:2,1:2])
  }
  sigma.eps ~ dexp(1)
  mu.theta[1] ~ dnorm(-9,1) 
  mu.theta[2] ~ dnorm(2,0.1)
  rho ~ dunif(-1,1)
  sigma.theta[1,1] ~ dexp(1) 
  sigma.theta[1,2] <- rho*sigma.theta[1,1]*sigma.theta[2,2] 
  sigma.theta[2,1] <- sigma.theta[1,2]
  sigma.theta[2,2] ~ dexp(1) 
}
"

#Finding suitable initial conditions
mean.test <- function(theta1,theta2,cycles=250,D0=9) {
  (D0^(1-theta2/2)+(1-theta2/2)*exp(theta1)*(pi)^(theta2/2)*cycles)^(2/(2-theta2))
}

mean.test(-9.3,3.5,cycles=225)
inits1 <- list(theta = matrix(c(rep(-9.3,68),rep(3.5,68)),ncol=2), sigma.eps=0.2)

# Run the MCMC
set.seed(123)
CrackLength.sim <- jags(
  data=c('crack.length','cycles','itemNum','cumItemNum'),
  parameters.to.save=c('theta','mu.theta','sigma.theta','sigma.eps'),
  model.file=textConnection(crack.growth.model),
  inits=list(inits1,inits1),
  n.iter=120000,
  n.burnin=20000,
  n.chains=2,
  n.thin=10
)

# Some diagnostics
gelman.diag(CrackLength.sim$BUGSoutput,multivariate = F)
post <- CrackLength.sim$BUGSoutput$sims.matrix
post.means <- apply(post,2,mean)

#################
#Plot fitted curves
for (i in 1:68) {
  eval(parse(text=paste('f',i,' <- function(x) {mean.test(post.means[8+',i,'],post.means[76+',i,'],x)}',sep='')))
}
gPlot <- ggplot(data, aes(x = cycles, y = length)) +
  labs(x = "Thousands of Cycles", y = "1/2 Crack Length in mm") +
  theme(legend.position="none") +
  coord_cartesian(xlim = c(0, 250), ylim = c(9.394, 36.55))
for (i in 1:68) {
  gPlot <- gPlot + geom_function(fun = eval(parse(text=paste('f',i,sep=''))), color = "gray")
}
gPlot <- gPlot + geom_point() + annotate("point", x = 0, y = 9) + geom_hline(yintercept=30)
#pdf(file="VirklerFittedDataPlot.pdf",width=5.5,height=3.667)
gPlot
#dev.off()

################
# Predict when crossing 30mm

# 1st draw a quantile for the error term
quantile.norm.draws <- runif(length(post[,1]))

# 2nd for the ith draw from the posterior: sample a random theta[1] & theta[2]
samp.theta <- function(i) {
  mu <- post[i,2:3]
  sigma <- matrix(post[i,5:8],ncol=2)
  theta <- rmvnorm(1,mean=mu,sigma=sigma)
  theta
}

# 3rd solve for the quantile at which the function will cross the 30mm threshold
solve.for.cycles <- function(cycles,quantile,theta,sigma.eps) {
  qnorm(quantile,mean.test(theta[1],theta[2],cycles=cycles),sd=sigma.eps) - 30
}
first.pass <- rep(NA,length(post[,1]))
for (i in 1:length(post[,1])) {
  theta <- samp.theta(i)
  upper.lim <- (-9^(1-theta[2]/2))/((1-theta[2]/2)*exp(theta[1])*pi^(theta[2]/2)) - 1
  if (theta[2] <= 2) {upper.lim <- 1000}
  first.pass[i] <- uniroot(solve.for.cycles,interval=c(10,upper.lim),
                           quantile=quantile.norm.draws[i],theta=theta,
                           sigma.eps=post[i,4])$root
}

#pdf(file="Virkler30FirstPassHist.pdf",width=5.5,height=3.667)
hist(first.pass,br=20,main='',xlab='Thousand of Cycles',ylab='',xlim=c(150,275),freq=F,ylim=c(0,.027))
#dev.off()

quantile(first.pass,c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995))


