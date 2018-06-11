## This code studies the behaviour of a robot based on a piece-wise Johnson-Dirichlet model. The pieces are given by change-points determined with the Adams-MacKay algorithm (adamsetal2007).

## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
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
logit <- function(a){log(a/(1-a))}
ilogit <- function(a){exp(a)/(1+exp(a))}

## probability-discrepancy measures
hellingerd <- function(a,b){sqrt(1-sum(sqrt(a*b)))}
kld <- function(a,b){
    temp <-  a * log(a/b)
    temp[is.nan(temp)] <- 0
    sum(temp)}
lkd <- function(a,b){
    temp <-  b * log(b/a)
    temp[is.nan(temp)] <- 0
    apply(temp,1,sum)}

## sequence of observations for a participant
observations <- function(participant, maxtrials=200){c(d[d$Participant==participant & d$Slot.Number==1,'Ball.Position'])[1:maxtrials]}


## sequence of distributions for a participant
## (consistency checked by for-loop)
## and find minimal non-zero value
distribution <- function(participant, maxtrials=200){matrix(d[d$Participant == participant & d$Trial <= maxtrials+1, c('Participant.Slot.Estimate')], maxtrials+1, N, byrow=T)}


## sequence of robot distributions: uses the changepoint & Johnson-Dirichlet model

## Define the probability of new changepoint, h(s,m)
#probchangepoint <- function(s){20^s*exp(-20)/factorial(s)}
#probchangepoint <- function(s,m){0.9-s/m}
#probchangepoint <- function(s){1}

## This gave results very close to participant 12
#probchangepoint <- function(s,m,n,params){1-dnorm(m-s-1,0,75)*sqrt(2*pi)*75*0.9}


## probchangepoint <- function(s,m,n,params){
##     if(s>m | m<0 | s<0){return(NA)}
##     if(s+m>=2*n/sqrt(3)){return(params[3])}
##     if(s+m>=sqrt(2/3)*n){return(params[2])}
##     return(params[1])
## }

## probchangepoint <- function(s,m,n,params){
##     ##if(s>m | m<0 | s<0){return(NA)}
##     if(s>=n/2){return(params[3])}
##     if(m>=n/2){return(params[2])}
##     return(params[1])
## }

## probchangepoint <- function(s,m,n,params){
##     if(s>m | m<0 | s<0){return(NA)}
##     if(m-s >=2*n/3){return(params[3])}
##     if(m-s >= n/3){return(params[2])}
##     return(params[1])
## }

## This defines h region-wise
probchangepoint <- function(s,m,n,nregions,params){
    ##if(s>m | m<0 | s<0){return(NA)}
    params[choose(floor(nregions*m/(n+1))+1,2)+1+nregions*s/(n+1)]
}

probchangepoint <- function(s,m,n,nregions,params){
    ##if(s>m | m<0 | s<0){return(NA)}
    ilogit(sum(unlist(lapply(0:nregions, function(i){lapply(0:i, function(j){params[choose(i+1,2)+1+j]*s^j*m^(i-j)})}))))
}


## sequences of means and stds
allmeans <- function(distrib){distrib %*% (1:N)}
allstds <- function(distrib){sqrt(distrib %*% (1:N)^2 - allmeans(distrib)^2)}

## simplified robotdistribution function
robotdistribution <- function(priorfrequencies,stubbornness,obs,params,nregions=3){
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
    h <- probchangepoint(0,1,n,nregions,params)
    CC <- AA * c(h,1-h)
    tempB <- t(t(L*priorfrequencies +
                   (cumfreqs[,2]-cumfreqs[,2:1]))/(L+(0:1)))
    ldistr[2,] <- (tempB %*% CC)/AA
    AA1 <- AA
    AA <- ldistr[2,obs[2]]
    BB <- tempB[obs[2],]

    for(m in 2:n){
        ## calculate h(s)
        h <- sapply(0:(m-1),function(i){probchangepoint(i,m,n,nregions,params)})
        
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

#######################################################################
## minimization of discrepancy between robot and participant

## Discrepancy between robot's and participant's sequences of distributions
discrepancy <- function(params,pdistr,obs,nregions=3,maxtrials=200,stubbornness=0.01){
    ##iparams <- ilogit(params) ## convert from (-inf,+inf) to (0,1)

    rdistr <- robotdistribution(pdistr[1,],stubbornness,obs,params,nregions)

    ## calculate total discrepancy
    mean(sapply(1:(maxtrials+1),function(i){kld(rdistr[i,],pdistr[i,])}))
}

## Algorithm to seek discrepancy minimum
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
        startpar <- rnorm((nregions^2+nregions)/2,0,100)
        optrobot <- optim(startpar,discrepancy,gr=NULL, pdistr=pdistr,obs=obs,nregions=nregions,maxtrials=maxtrials,stubbornness=0.01)
        if(optrobot$value < maxval){maxval <- optrobot$value
            maxpars <- optrobot$par
            region <- startpar}
    }
    list(par=maxpars,value=maxval,region=region)}

## Function to plot the resulting h(s,m) from optimization
plotminparams <- function(nregions,params,label='',maxtrials=200){
	nobs <- 0:maxtrials
	probd <- function(s,m,n,nregions,params){
	if(s>m | m<0 | s<0){return(NA)}
	probchangepoint(s,m,n,nregions,params)
	}
	
	pmatrix <- sapply(nobs,function(i){sapply(nobs,function(j){probd(i,j,	maxtrials,nregions,params)})})
	pdf(paste0(plotsdir,label,'optimh.pdf'),width = 148*mmtoin, height = 148*1*mmtoin)
	image2D(pmatrix,x=nobs,y=nobs,xlab='m',ylab='s',zlim=c(0,1))
	dev.off()
}



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

    rdistr <- robotdistribution(pdistr[1,],stubbornness,obs,params,nregions)

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


    ## plot h(s,m)
    plotminparams(nregions,params,label=paste0(plotsdir,label,'_partc',participant,'.pdf'),maxtrials=maxtrials)


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


