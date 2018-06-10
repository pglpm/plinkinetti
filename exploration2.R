## This code studies the behaviour of a robot based on a piece-wise Johnson-Dirichlet model. The pieces are given by change-points determined with the Adams-MacKay algorithm (adamsetal2007).

## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
dev.off()
mmtoin <- 0.0393701
plotsdir <- './comparisons2/'

## load all data
dpath = "./_data/"
fs = dir(path = dpath,pattern = "barData")
d <- read.csv(paste0(dpath,fs[1]))
for (f in fs[2:40]) {
dprime <- read.csv(paste0(dpath, f))
d <- rbind(d,dprime)
rm(dprime)
}
d <- d[,1:(ncol(d)-2)]
person <- readPNG('face2.png')
robot <- readPNG('robot1.png')
robot2 <- readPNG('robot2.png')

## number of slots and sequence of zeros for later use
N <- 40
slots <- integer(N)
maxstd <- sqrt((1+N^2)*0.5-(N+1)*0.5)

## Definitions of various functions

## sequence of observations for a participant
observations <- function(participant, maxtrials=200){c(d[d$Participant==participant & d$Slot.Number==1,'Ball.Position'])[1:maxtrials]}


## sequence of distributions for a participant
## (consistency checked by for-loop)
## and find minimal non-zero value
distribution <- function(participant, maxtrials=200){matrix(d[d$Participant == participant & d$Trial <= maxtrials+1, c('Participant.Slot.Estimate')], maxtrials+1, N, byrow=T)}


## sequence of robot distributions: uses the changepoint & Johnson-Dirichlet model

## first define the probability of new changepoint, h(s)
#probchangepoint <- function(s){20^s*exp(-20)/factorial(s)}
probchangepoint <- function(s,m){1-dnorm(m-s-1,0,75)*sqrt(2*pi)*75*0.9}
#probchangepoint <- function(s,m){0.9-s/m}
#probchangepoint <- function(s){1}

robotdistribution <- function(priorfrequencies,stubbornness,obs,fprobchangepoint=probchangepoint){
    n <- length(obs)

    ## sequence of cumulative frequencies:
    ## column i are the frequencies up to observation i inclusive
    cumfreqs <- sapply(1:(n+1), function(i) sapply(1:N, function(j) sum(c(0,obs)[1:i] == j)))

    ## Johnson-Dirichlet model:
    ## stubbornness parameter and ensuing sequence
    L <- stubbornness # + (0:(n+1))
    ## initialize
    ldistr <- matrix(NA, n+1, N)
    probR <- list()

    ## m=0
    ldistr[1,] <- priorfrequencies
    AA <- ldistr[1,obs[1]]
    probR[[1]] <- 1

    ## m=1
    h <- fprobchangepoint(0,1)
    CC <- AA * c(h,1-h)
    tempB <- t(t(L*priorfrequencies +
                   (cumfreqs[,2]-cumfreqs[,2:1]))/(L+(0:1)))
    ldistr[2,] <- (tempB %*% CC)/AA
    probR[[2]] <- CC/AA
    AA1 <- AA
    AA <- ldistr[2,obs[2]]
    BB <- tempB[obs[2],]

    for(m in 2:n){
        ## calculate h(s)
        h <- sapply(0:(m-1),function(i){fprobchangepoint(i,m)})
        
        ## calculate C_{m+1}(r)
        CC <- c(sum(h * BB * CC), # r=0
        (1-h) * BB * CC # r>0
        )/AA1

        ## calculate B_{m+1}(r, d)
        tempB  <- t(t(L*priorfrequencies +
                   (cumfreqs[,m+1]-cumfreqs[,(m+1):1]))/(L+(0:m)))

        ## calculate A_{m+1}(d)
        ldistr[m+1,] <- (tempB %*% CC)/AA

        probR[[m+1]] <- CC/AA
        AA1 <- AA
        AA <- ldistr[m+1,obs[m+1]]
        BB <- tempB[obs[m+1],]
    }

    list(ldistr,probR)}

## sequences of means and stds
allmeans <- function(distrib){distrib %*% (1:N)}
allstds <- function(distrib){sqrt(distrib %*% (1:N)^2 - allmeans(distrib)^2)}


## Define a function that plots:
## - sequences of overlap, relative entropies, means, stds
## - histograms for a selected set of trials
## and that outputs the final distributions, rel-entropies, overlaps, means, stds
comparison <- function(participant,maxtrials=200,trialstoshow=c(1:10, 95:105, 190:200),stubbornness=0.01,label='',params=rep(0.5,6),nregions=3){
    n <- maxtrials ## number of trials to consider
    
    obs <- observations(participant,n)

    ## sequence of frequency parameters of the changepoint-JD model. At
    ## every "reset" we set the first equal to the participant's initial
    ## distribution plus a value equal to 1/1000 of the minimum nonzero value
    ## ever assigned
    tdistr <- distribution(participant,n)
    mini <- tdistr + min(c(tdistr[tdistr>0]))/1000
    pdistr <- t(sapply(1:(n+1),function(i){mini[i,]/sum(mini[i,])}))

    rdistr <- robotdistributiondiscr(pdistr[1,],stubbornness,obs,params,nregions)

    ## calculate the overlap and relative entropy of the participant's distr.
    temp <-  rdistr * log(rdistr/pdistr)
    temp[is.nan(temp)] <- 0
    rentropy <- apply(temp,1,sum)
    overlap <- diag(rdistr %*% t(pdistr))

    ## sequence of means and stds
    meanperson <- allmeans(pdistr)
    stdperson <- allstds(pdistr)
    meanrobot <- allmeans(rdistr)
    stdrobot <- allstds(rdistr)

## ## plot obzervations
## df <- data.frame(x=2:(n+1), y1=obs)
##     maxy <- N
##     miny <- 1
## g <- ggplot() + theme_classic() 
## g <- g + geom_point(data=df, aes(x,y1), colour='black') +
##     geom_line(data=df, aes(x,y1), colour='black') +
##     xlim(1,n+1) + ylim(miny,maxy) +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='observations',
##          title=paste0('participant ', participant))
## pdfname <- paste0(plotsdir,'observ_',participant,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()
## ##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)

## ## plot the overlap
## df <- data.frame(x=1:(n+1), y=overlap)
## g <- ggplot() + theme_classic() 
## g <- g + #geom_point(data=df, aes(x,y), colour=mygreen) +
##     geom_line(data=df, aes(x,y), colour=mygreen) +
##     xlim(1,n+1) + # ylim(0,1) +
##     theme_classic() +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='overlap',
##          title=paste0('participant ', participant,', stub ',stubbornness))
## pdfname <- paste0(plotsdir,'overlap_',participant,'-stub_',stubbornness,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()

## ## plot the relative entropy
## df <- data.frame(x=1:(n+1), y=rentropy)
## g <- ggplot() + theme_classic() 
## g <- g + #geom_point(data=df, aes(x,y), colour=myyellow) +
##     geom_line(data=df, aes(x,y), colour=myyellow) +
##     xlim(1,n+1) +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='relative entropy',
##          title=paste0('participant ', participant,', stub ',stubbornness))
## pdfname <- paste0(plotsdir,'rentropy_',participant,'-stub_',stubbornness,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()
## ##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)

## plot sequence of means
df <- data.frame(x=1:(n+1), y1=meanperson, y2=meanrobot, y3=c(NA,obs))
    maxy <- N #max(df$y1, df$y2)
    miny <- 1 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (n+1)/20/300*271
    personwidth <- (n+1)/20/300*223
g <- ggplot() + theme_classic() 
g <- g + annotation_raster(person, xmax = n+1, xmin = n+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = n+1, xmin = n+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
g <- g + geom_line(data=df, aes(x,y3), colour='black', alpha=0.33)
    g <- g + #geom_point(data=df, aes(x,y1), colour=myred, alpha=0.33) +
    geom_line(data=df, aes(x,y1), colour=myred, alpha=0.75) +
    #geom_point(data=df, aes(x,y2), colour=myblue, alpha=0.33) +
    geom_line(data=df, aes(x,y2), colour=myblue, alpha=0.75) +
    xlim(1,n+1) + ylim(miny,maxy) +
    theme(aspect.ratio=0.5) +
    labs(x='trial',y='mean',
         title=paste0('participant #', participant,', divisions = ',nregions))
g <- g + geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
pdfname <- paste0(plotsdir,label,'means_partc',participant,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height=148*0.6, units='mm', dpi = 300)
#ggsave(pdfname, width = 148, height = 148*0.2, units='mm', dpi = 300)
dev.off()

## plot sequence of stds
df <- data.frame(x=1:(n+1), y1=stdperson, y2=stdrobot)
    maxy <- maxstd #max(df$y1, df$y2)
    miny <- 0 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (n+1)/20/300*271
    personwidth <- (n+1)/20/300*223
g <- ggplot() + theme_classic() 
g <- g + annotation_raster(person, xmax = n+1, xmin = n+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = n+1, xmin = n+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
g <- g + #geom_point(data=df, aes(x,y1), colour=myred, alpha=0.33) +
    geom_line(data=df, aes(x,y1), colour=myred, alpha=0.75) +
    #geom_point(data=df, aes(x,y2), colour=myblue, alpha=0.33) +
    geom_line(data=df, aes(x,y2), colour=myblue, alpha=0.75) +
    xlim(1,n+1) + ylim(miny,maxy) +
    theme(aspect.ratio=0.5) +
    labs(x='trial',y='std',
         title=paste0('participant #', participant,', divisions = ',nregions))
g <- g + geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
pdfname <- paste0(plotsdir,label,'stds_partc',participant,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
dev.off()
##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)


## plot histograms for a range of trials
if(!is.null(trialstoshow)){rangehist <- trialstoshow
maxhist <- max(c(pdistr,rdistr))
pdfname <- paste0(plotsdir,label,'histogram_',rangehist[1],'-',rangehist[length(rangehist)],'_partc',participant,'.pdf')
pdf(pdfname,width = 148*mmtoin, height = 148*0.6*mmtoin)
for(j in 1:length(rangehist)){
    i <- rangehist[j]
    df <- data.frame(x=rep(1:N, 2),
                     who=rep(c('person','robot'), each=N),
                     y=c(pdistr[i,], rdistr[i,]))
    maxy <- max(df$y)
    miny <- 0
    iconheight <- (maxy-miny)/10
    robotwidth <- N/20/300*271
    personwidth <- N/20/300*223
    g <- ggplot() + theme_classic() 
g <- g + annotation_raster(person, xmax = N, xmin = N-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = N, xmin = N-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
    g <- g + geom_bar(data=df, aes(x=x,y=y,fill=who), alpha=0.5,stat='identity', position='identity') +
         scale_fill_manual(values=c(myred,myblue)) +
        scale_x_continuous(breaks=c(1,seq(5,N,length.out=8))) + 
        ylim(miny,maxy) +
        theme(aspect.ratio=0.5,
              legend.title=element_blank(),
              legend.background=element_blank(),
              legend.justification=c(1,1),
              legend.position='none',
              ) +
        labs(x='slot',y='probability',
             title=paste0('participant #', participant,', divisions = ',nregions,', trial ',i,', overlap ',signif(overlap[i],2),', rel-entr. ',signif(rentropy[i],2)))
g <- g + geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
print(g)
}
dev.off()}

    return(list(distr=pdistr,robotdistr=rdistr,rentropy=rentropy,overlap=overlap,meanperson=meanperson,meanrobot=meanrobot,stdperson=stdperson,stdrobot=stdrobot))
}

## function to assess speed of robot re-learning
comparerobots <- function(index,obs,changepoint,priorfrequencies=rep(1/N,N),stubbornness=0.1,label=''){
    n <- length(obs)
    distr1 <- robotdistribution(priorfrequencies,stubbornness,obs)[[1]]
    distr2 <- robotdistribution(distr1[changepoint-1,],stubbornness,obs[changepoint:n])[[1]]
    pdfname <- paste0(plotsdir,label,'compare_robots_',index,'-stub_',stubbornness,'.pdf')
    df0 <- data.frame(x=2:(n+1),y=obs)
    df1 <- data.frame(x=1:(n+1),y=allstds(distr1))
    df2 <- data.frame(x=changepoint:(n+1),y=allstds(distr2))
    ##
    maxy <- 15#maxstd #max(df$y1, df$y2)
    miny <- 0 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (n+1)/20/300*271
    g <- ggplot() + theme_bw()
    g <- g + annotation_raster(robot, xmax = n+1, xmin = n+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot2, xmax = n+1, xmin = n+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
#    g <- g + geom_line(data=df0, aes(x,y), colour='black', alpha=0.33)
    g <- g + geom_line(data=df1, aes(x,y), colour=myblue, alpha=0.75)
    g <- g + geom_line(data=df2, aes(x,y), colour=myredpurple, alpha=0.75)
    g <- g + xlim(1,n+1) + ylim(miny,maxy) +
        theme(aspect.ratio=0.5) +
        labs(x='trial', y='std', title=index)
    g <- g + geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myblue, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myredpurple, alpha=0.5)
    save_plot(pdfname, g, base_width = 148, base_height=148*0.6, units='mm', dpi = 300)
#ggsave(pdfname, width = 148, height = 148*0.2, units='mm', dpi = 300)
    dev.off()
    ## relative entropy and overlap
    distr1t <- distr1[changepoint:(n+1),]
    temp <-  distr2 * log(distr2/distr1t)
    temp[is.nan(temp)] <- 0
    rentropy <- apply(temp,1,sum)
    temp <-  distr1t * log(distr1t/distr2)
    temp[is.nan(temp)] <- 0
    rentropy2 <- apply(temp,1,sum)
    overlap <- diag(distr2 %*% t(distr1t))/sqrt(diag(distr2 %*% t(distr2)) * diag(distr1t %*% t(distr1t)))
    return(list(rentropy,
             rentropy2,
             overlap))
}

#######################################################################
## minimization of discrepancy between robot and participant

logit <- function(a){log(a/(1-a))}
ilogit <- function(a){exp(a)/(1+exp(a))}

## changepoint probability
## probchangepointdiscr <- function(s,m,n,params){
##     if(s>m | m<0 | s<0){return(NA)}
##     if(s+m>=2*n/sqrt(3)){return(params[3])}
##     if(s+m>=sqrt(2/3)*n){return(params[2])}
##     return(params[1])
## }

## probchangepointdiscr <- function(s,m,n,params){
##     ##if(s>m | m<0 | s<0){return(NA)}
##     if(s>=n/2){return(params[3])}
##     if(m>=n/2){return(params[2])}
##     return(params[1])
## }

## probchangepointdiscr <- function(s,m,n,params){
##     if(s>m | m<0 | s<0){return(NA)}
##     if(m-s >=2*n/3){return(params[3])}
##     if(m-s >= n/3){return(params[2])}
##     return(params[1])
## }

probchangepointdiscr <- function(s,m,n,nregions,params){
    ##if(s>m | m<0 | s<0){return(NA)}
    params[choose(floor(nregions*m/(n+1))+1,2)+1+nregions*s/(n+1)]
}

## probchangepointdiscr <- function(s,m,n,params){1-dnorm(m-s-1,0,75)*sqrt(2*pi)*75*0.9}


## simplified robotdistribution function
robotdistributiondiscr <- function(priorfrequencies,stubbornness,obs,params,nregions=3){
    n <- length(obs)

    ## sequence of cumulative frequencies:
    ## column i are the frequencies up to observation i inclusive
    cumfreqs <- sapply(1:(n+1), function(i) sapply(1:N, function(j) sum(c(0,obs)[1:i] == j)))

    ## Johnson-Dirichlet model:
    ## stubbornness parameter and ensuing sequence
    L <- stubbornness # + (0:(n+1))
    ## initialize
    ldistr <- matrix(NA, n+1, N)

    ## m=0
    ldistr[1,] <- priorfrequencies
    AA <- ldistr[1,obs[1]]

    ## m=1
    h <- probchangepointdiscr(0,1,n,nregions,params)
    CC <- AA * c(h,1-h)
    tempB <- t(t(L*priorfrequencies +
                   (cumfreqs[,2]-cumfreqs[,2:1]))/(L+(0:1)))
    ldistr[2,] <- (tempB %*% CC)/AA
    AA1 <- AA
    AA <- ldistr[2,obs[2]]
    BB <- tempB[obs[2],]

    for(m in 2:n){
        ## calculate h(s)
        h <- sapply(0:(m-1),function(i){probchangepointdiscr(i,m,n,nregions,params)})
        
        ## calculate C_{m+1}(r)
        CC <- c(sum(h * BB * CC), # r=0
        (1-h) * BB * CC # r>0
        )/AA1

        ## calculate B_{m+1}(r, d)
        tempB  <- t(t(L*priorfrequencies +
                   (cumfreqs[,m+1]-cumfreqs[,(m+1):1]))/(L+(0:m)))

        ## calculate A_{m+1}(d)
        ldistr[m+1,] <- (tempB %*% CC)/AA

        AA1 <- AA
        AA <- ldistr[m+1,obs[m+1]]
        BB <- tempB[obs[m+1],]
    }
    ldistr}




## probability-discrepancy measures
hellingerd <- function(a,b){sqrt(1-sum(sqrt(a*b)))}
kld <- function(a,b){
    temp <-  a * log(a/b)
    temp[is.nan(temp)] <- 0
    sum(temp)}
kl2d <- function(a,b){
    temp <-  b * log(b/a)
    temp[is.nan(temp)] <- 0
    apply(temp,1,sum)}


## function to calculate "discrepancy" between participant's and robot's sequences of distributions
discrepancypart <- function(params,participant,nregions=3,maxtrials=200,stubbornness=0.01){
    iparams <- ilogit(params) ## convert from (-inf,+inf) to (0,1)
    n <- maxtrials ## number of trials to consider
    obs <- observations(participant,n)
    distr <- distribution(participant,n)

    ## sequence of frequency parameters of the JD model. At every "reset"
    ## we set the first equal to the participant's initial distribution
    ## plus a value equal to 1/100 of the minimum nonzero value ever
    ## assigned (to avoid zero probabilities in the robot)
    mini <- min(c(distr[distr>0]))
    temp <- distr[1,] + mini/100
    temp <- temp/sum(temp)
    ldistr <- robotdistributiondiscr(temp,stubbornness,obs,iparams,nregions)

    ## calculate total discrepancy
    mean(sapply(1:n,function(i){kld(ldistr[i,],distr[i,])}))
}

discrepancy <- function(params,pdistr,obs,nregions=3,maxtrials=200,stubbornness=0.01){
    iparams <- ilogit(params) ## convert from (-inf,+inf) to (0,1)

    rdistr <- robotdistributiondiscr(pdistr[1,],stubbornness,obs,iparams,nregions)

    ## calculate total discrepancy
    mean(sapply(1:(maxtrials+1),function(i){kld(rdistr[i,],pdistr[i,])}))
}

reducediscrepancy <- function(participant,maxtrials,nregions=3,startpoints=10,seed=999){
    set.seed(seed)
    n <- maxtrials
    obs <- observations(participant,n)
    tdistr <- distribution(participant,n)

    ## sequence of frequency parameters of the JD model. At every "reset"
    ## we set the first equal to the participant's initial distribution
    ## plus a value equal to 1/100 of the minimum nonzero value ever
    ## assigned (to avoid zero probabilities in the robot)
    mini <- tdistr + min(c(tdistr[tdistr>0]))/1000
    pdistr <- t(sapply(1:(n+1),function(i){mini[i,]/sum(mini[i,])}))

    maxval <- Inf
    for(i in 1:startpoints){
        startpar <- logit(runif((nregions^2+nregions)/2))
        optrobot <- optim(startpar,discrepancy,gr=NULL, pdistr=pdistr,obs=obs,nregions=nregions,maxtrials=maxtrials,stubbornness=0.01)
        if(optrobot$value < maxval){maxval <- optrobot$value
            maxpars <- optrobot$par
            region <- startpar}
    }
    list(par=ilogit(maxpars),value=maxval,region=ilogit(region))}

reducediscrepancy3 <- function(participant,maxtrials){
    maxval=Inf
    prange <- logit(2/3)*c(-1,1)

    for(a1 in prange){
        for(a2 in prange){
            for(a3 in prange){
                optrobot <- optim(c(a1,a2,a3),discrepancypart,gr=NULL, participant=participant,maxtrials=maxtrials,stubbornness=0.01)
                if(optrobot$value < maxval){maxval <- optrobot$value
                    maxpars <- optrobot$par
                    region <- c(a1,a2,a3)}
            }}}
list(par=ilogit(maxpars),value=maxval,region=region)}


## as previous function, but without plots
comparerobotssimple <- function(obs,changepoint,priorfrequencies=rep(1/N,N),stubbornness=0.1){
    n <- length(obs)
    distr1 <- robotdistribution(priorfrequencies,stubbornness,obs)[[1]]
    distr2 <- robotdistribution(distr1[changepoint-1,],stubbornness,obs[changepoint:n])[[1]]
    ## std, relative entropy, and overlap
    std1 <- allstds(distr1)
    std2 <- allstds(distr2)
    distr1t <- distr1[changepoint:(n+1),]
    temp <-  distr2 * log(distr2/distr1t)
    temp[is.nan(temp)] <- 0
    rentropy <- apply(temp,1,sum)
    temp <-  distr1t * log(distr1t/distr2)
    temp[is.nan(temp)] <- 0
    rentropy2 <- apply(temp,1,sum)
    overlap <- diag(distr2 %*% t(distr1t))/sqrt(diag(distr2 %*% t(distr2)) * diag(distr1t %*% t(distr1t)))
    return(list(std1,std2,rentropy,
             rentropy2,
             overlap))
}

## average results of previous function over many experiments
averagecomparerobots <- function(index,numexperiments,meanstd1,meanstd2,numobs,changepoint,priorfrequencies=rep(1/N,N),stubbornness=0.1,seed=666,label=''){
    set.seed(seed)
    obs1 <- matrix(round(rnorm(numexperiments*(changepoint-1),meanstd1[1],meanstd1[2])),numexperiments,changepoint-1)
    obs1[obs1<1] <- 1
    obs1[obs1>N] <- N
    obs2 <- matrix(round(rnorm(numexperiments*(numobs-changepoint+1),meanstd2[1],meanstd2[2])),numexperiments,numobs-changepoint+1)
    obs2[obs2<1] <- 1
    obs2[obs2>N] <- N
    obsall <- cbind(obs1,obs2)
    ##initialize matrices
    mstd1 <- matrix(NA,numexperiments,numobs+1)
    mstd2 <- matrix(NA,numexperiments,numobs+1-changepoint+1)
    mrentropy <- mstd2
    mrentropy2 <- mstd2
    moverlap <- mstd2
    for(i in 1:numexperiments){
        temp <- comparerobotssimple(obsall[i,],changepoint)
        mstd1[i,] <- temp[[1]]
        mstd2[i,] <- temp[[2]]
        mrentropy[i,] <- temp[[3]]
        mrentropy2[i,] <- temp[[4]]
        moverlap[i,] <- temp[[5]]
    }
    ##plots
    n <- numobs
    pdfname <- paste0(plotsdir,label,'avg_compare_robots_',index,'-stub_',stubbornness,'.pdf')
    df1 <- data.frame(x=1:(n+1),y=apply(mstd1,2,mean))
    df2 <- data.frame(x=changepoint:(n+1),y=apply(mstd2,2,mean))
    ##
    maxy <- 15#maxstd #max(df$y1, df$y2)
    miny <- 0 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (n+1)/20/300*271
    g <- ggplot() + theme_bw()
    g <- g + annotation_raster(robot, xmax = n+1, xmin = n+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot2, xmax = n+1, xmin = n+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
#    g <- g + geom_line(data=df0, aes(x,y), colour='black', alpha=0.33)
    g <- g + geom_line(data=df1, aes(x,y), colour=myblue, alpha=0.75)
    g <- g + geom_line(data=df2, aes(x,y), colour=myredpurple, alpha=0.75)
    g <- g + xlim(1,n+1) + ylim(miny,maxy) +
        theme(aspect.ratio=0.5) +
        labs(x='trial', y='std', title=index)
    g <- g + geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myblue, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myredpurple, alpha=0.5)
    save_plot(pdfname, g, base_width = 148, base_height=148*0.6, units='mm', dpi = 300)
#ggsave(pdfname, width = 148, height = 148*0.2, units='mm', dpi = 300)
    dev.off()
    return(list(mstd1,mstd2,mrentropy,mrentropy2,moverlap))
}
