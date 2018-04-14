## all data
dpath = "./data/"
fs = dir(path = dpath,pattern = "barData")
d <- read.csv(paste0(dpath,fs[1]))
for (f in fs[2:40]) {
dprime <- read.csv(paste0(dpath, f))
d <- rbind(d,dprime)
rm(dprime)
}
d <- d[,1:(ncol(d)-2)]

## number of trials to consider
n <- 200

## choose a participant
x <- 3

## sequence of observations for that participant
obs <- c(d[d$Participant==x & d$Slot.Number==1,'Ball.Position'])[1:n]

## sequence of distributions for that participant
## (consistency checked by for-loop)
distr <- matrix(d[d$Participant==x & d$Trial <= n, c('Participant.Slot.Estimate')],n,40,byrow=T)



