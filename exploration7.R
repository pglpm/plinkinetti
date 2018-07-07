## Johnson;Dirichlet robot with changepoints:
## optimization with 4 variables:
## stubbornness and logistic-linear changepoint function


## libraries and colour-blind palette from http://www.sron.nl/~pault/
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('dfoptim')
mypurpleblue <- '#4477AA'
myblue <- '#66CCEE'
mygreen <- '#228833'
myyellow <- '#CCBB44'
myred <- '#EE6677'
myredpurple <- '#AA3377'
mygrey <- '#BBBBBB'
mycolours <- c(myblue, myred, mygreen, myyellow, myredpurple, mypurpleblue, mygrey, 'black')
palette(mycolours)
barpalette <- colorRampPalette(c(mypurpleblue,'white',myredpurple),space='Lab')
barpalettepos <- colorRampPalette(c('white','black'),space='Lab')
dev.off()
mmtoin <- 0.0393701
plotsdir <- './comparisons7c/'

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
logistic <- function(a){exp(a)/(1+exp(a))}

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
jsd <- function(a,b){
    c <- (a+b)/2
    temp1 <-  a * log(a/c)
    temp1[is.nan(temp1)] <- 0
    temp2 <-  b * log(b/c)
	temp2[is.nan(temp2)] <- 0
    sum(temp1+temp2)/2}

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
## probchangepoint <- function(s,m,n,nregions,params){
##     ##if(s>m | m<0 | s<0){return(NA)}
##     params[choose(floor(nregions*m/(n+1))+1,2)+1+nregions*s/(n+1)]
## }

## note for the case nregions=1:
## = logistic(param[1] + m*param[2] + s*param[3]),
## where m and s range from -1 to 1 inclusive
## probchangepoint <- function(s,m,n,nregions,params){
##     ##if(s>m | m<0 | s<0){return(NA)}
##     logistic(sum(unlist(lapply(0:nregions, function(i){lapply(0:i, function(j){params[choose(i+1,2)+1+j]*(((2*s-n)/n)^j)*(((2*m-n)/n)^(i-j))})}))))
## }

probchangepoint <- function(s,m,n,params){
    ##if(s>m | m<0 | s<0){return(NA)}
    logistic(params[1] + params[2]*(2*m-n)/n + params[3]*(2*s-n)/n)
}

## sequences of means and stds
allmeans <- function(distrib){distrib %*% (1:N)}
allstds <- function(distrib){sqrt(distrib %*% (1:N)^2 - allmeans(distrib)^2)}

## simplified robotdistribution function
robotdistribution <- function(priorfrequencies,allparams,obs){
    n <- length(obs)
    stubbornness <- exp(allparams[1])
    params <- allparams[-1]
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
    h <- probchangepoint(0,1,n,params)
    CC <- AA * c(h,1-h)
    tempB <- t(t(L*priorfrequencies +
                   (cumfreqs[,2]-cumfreqs[,2:1]))/(L+(0:1)))
    ldistr[2,] <- (tempB %*% CC)/AA
    AA1 <- AA
    AA <- ldistr[2,obs[2]]
    BB <- tempB[obs[2],]

    for(m in 2:n){
        ## calculate h(s)
        h <- sapply(0:(m-1),function(i){probchangepoint(i,m,n,params)})
        
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

## eliminate zero-probabilities from distributions
regularizedistr <- function(tdistr,n=NULL,fraction=1000000){
    if(length(tdistr)==1){tdistr <- distribution(tdistr,n)}
    if(is.vector(tdistr)){dim(tdistr) <- c(1,length(tdistr))}
    n <- dim(tdistr)[1]-1
    ## we set the distributions plus a value equal to 1/100 of the minimum
    ## nonzero value ever assigned
    toadd <- min(c(tdistr[tdistr>0]))/fraction
    pdistr <- t(sapply(1:(n+1),function(i){
        mini <- tdistr[i,]
        if(min(mini)==0){
            mini <- mini + toadd}
            mini/sum(mini)}))
    pdistr}

## Discrepancy between robot's and participant's sequences of distributions
discrepancy <- function(allparams,pdistr,obs,maxtrials=200,initial.distr=NULL){
    ##iparams <- logistic(params) ## convert from (-inf,+inf) to (0,1)
    ## check if pdistr is a participant's ID
    if(length(pdistr)==1){pdistr <- regularizedistr(pdistr,maxtrials)}
    ## check if obs is a participant's ID
    if(length(obs)==1){obs <- observations(obs,maxtrials)}

    ## robot's reset distribution
    if(length(initial.distr)==0){## use the first provided
        initial.distr <- pdistr[1,]
    }
    else{## use another participant's or custom
        initial.distr <- regularizedistr(initial.distr,maxtrials)[1,]
    }

    rdistr <- robotdistribution(initial.distr,allparams,obs)
    if(anyNA(rdistr)){return(NA)}

    ## calculate total discrepancy
    mean(sapply(1:(maxtrials+1),function(i){jsd(rdistr[i,],pdistr[i,])}))
}

## fast version, assumes pdistr, obs, and initial.distr are vectors
discrepancyfast <- function(allparams,pdistr,obs){
    
    rdistr <- robotdistribution(pdistr[1,],allparams,obs)
    if(anyNA(rdistr)){return(NA)}
    
    ## calculate total discrepancy
    mean(sapply(1:(length(obs)+1),function(i){jsd(rdistr[i,],pdistr[i,])}))
}

## Algorithm to seek discrepancy minimum using optim
## max + kl: 600
## max + js: 500
## mean + kl: 370
## mean + js: 360
reducediscrepancy <- function(participant,maxtrials,startpoints=4,seed=999){
    set.seed(seed)
    obs <- observations(participant,maxtrials)
    pdistr <- regularizedistr(participant,maxtrials)
    startstub <- c(log(0.1),log(10),log(1),rnorm(startpoints,0,2)) 

    maxval <- Inf
    maxpars <- NA
    region <- NA
    for(i in 1:startpoints){
        message("iteration ",i)
        testNA <- TRUE
        counti <- 0
        while(testNA){
            cat('looking for acceptable starting point...')
            startpar <- c(startstub[i],rnorm(3,0,1/sqrt(3)))
            testNA <- anyNA(discrepancyfast(startpar,pdistr,obs))
            counti <- counti + 1
            if(counti>10){message("can't find starting point")
            return(NA)}
        }
        cat('found (',exp(startpar[1]),") ",startpar,' in ',counti, ' trials.\n')
        optrobot <- optim(startpar,discrepancyfast,gr=NULL, pdistr=pdistr,obs=obs,control=list(maxit=2500))
        conv <- optrobot$convergence
        if(conv > 0){message("iteration ",i," didn't converge: ",conv)
        details <- optrobot}
        if(conv == 0 & optrobot$value < maxval){
            message('iteration ',i,' accepted:')
            maxval <- optrobot$value
            maxpars <- optrobot$par
            cat('discr = ',maxval,',  values = (',exp(maxpars[1]),') ',maxpars,'\n')
            region <- startpar
            details <- optrobot}
    }
    list(par=maxpars,value=maxval,region=region,details=details,convergence=conv)}

## Algorithm to seek discrepancy minimum using nlm
## time for nregions=2: 1107.021    0.101 1106.673
## time region 1, jsd: 1137.590    1.277 1138.407
reducediscrepancyhjk <- function(participant,maxtrials,nregions=1,startpoints=10,seed=999){
    set.seed(seed)
    n <- maxtrials
    nparams <- ((nregions+1)^2+nregions+1)/2
    obs <- observations(participant,n)
    pdistr <- regularizedistr(participant,n)

    maxval <- Inf
    maxpars <- NA
    region <- NA
    for(i in 1:startpoints){
        message("iteration ",i)
        testNA <- TRUE
        counti <- 0
        while(testNA){
            message('looking for acceptable starting point...')
            startpar <- rnorm(nparams,0,1/sqrt(nparams))
            testNA <- anyNA(discrepancy(startpar,pdistr,obs,nregions,maxtrials,0.01))
            counti <- counti + 1
        }
        cat('found (',exp(startpar[1]),") ",startpar,' in ',counti, ' trials.\n')
        optrobot <- hjk(startpar,discrepancy, pdistr=pdistr,obs=obs,nregions=nregions,maxtrials=maxtrials,stubbornness=0.01)##,control=list(maxit=2500))
        if(optrobot$convergence > 0){message("iteration ",i," didn't converge: ",optrobot$convergence)
        details <- optrobot}
        if(optrobot$convergence == 0 & optrobot$value < maxval){
            message('iteration ',i,' accepted:')
            maxval <- optrobot$value
            maxpars <- optrobot$par
            print(maxval)
            print(maxpars)
            region <- startpar
            details <- optrobot}
    }
    list(par=maxpars,value=maxval,region=region,details=details)}

## Function to plot the resulting h(s,m) from optimization
plotminparams <- function(allparams,label='',maxtrials=200){
    params <- allparams[-1]
    nobs <- 0:maxtrials
    probd <- function(s,m,n,params){
	if(s>m | m<0 | s<0){return(NA)}
	probchangepoint(s,m,n,params)
    }
	
    pmatrix <- sapply(nobs,function(i){sapply(nobs,function(j){
        probd(i,j,maxtrials,params)})})
    png(paste0(plotsdir,label,'_optimh.png'))
    image2D(pmatrix,x=nobs,y=nobs,xlab='m',ylab='s',zlim=c(0,1),
                col=barpalettepos(25))
    dev.off()
}

## Function to plot differences between two distributions
plotdiscrseq <- function(rdistr,pdistr,label=''){
    pmatrix <- rdistr-pdistr
    dims <- dim(pmatrix)
   ssqrt <- function(x){sign(x)*sqrt(abs(x))}
##    ssqrt <- function(x){sign(x)*(abs(x)^(1/3))}
    
    png(paste0(plotsdir,label,'_diff.png'))
    image2D(ssqrt(t(pmatrix)),x=1:(dims[2]),y=1:(dims[1]),xlab='slot',ylab='observation',zlim=c(-1,1), col=barpalette(51))
    dev.off()
}



## Define a function that plots:
## - sequences of overlap, relative entropies, means, stds
## - histograms for a selected set of trials
## and that outputs the final distributions, rel-entropies, overlaps, means, stds
comparison <- function(participant,maxtrials=200,trialstoshow=c(1:4, 99:102, 197:200),label='',params=rep(0.5,4),initial.distr=NULL,graphs=TRUE){
    obs <- observations(participant,maxtrials)

    ## sequence of frequency parameters of the changepoint-JD model. At
    ## every "reset" we set the first equal to the participant's initial
    ## distribution plus a value equal to 1/1000 of the minimum nonzero value
    ## ever assigned
    pdistr <- regularizedistr(participant,maxtrials)

    ## robot's reset distribution
    if(length(initial.distr)==0){## use this participant's
        initial.distr <- pdistr[1,]
    }
    else{## use another participant's or custom
        initial.distr <- regularizedistr(initial.distr,maxtrials)[1,]
    }
    
    rdistr <- robotdistribution(initial.distr,params,obs)

    ## calculate the overlap and relative entropy of the participant's distr.
    jsseq <- sapply(1:(maxtrials+1),function(i){jsd(rdistr[i,],pdistr[i,])})
    klseq <- sapply(1:(maxtrials+1),function(i){kld(rdistr[i,],pdistr[i,])})

    ## sequence of means and stds
    meanperson <- allmeans(pdistr)
    stdperson <- allstds(pdistr)
    meanrobot <- allmeans(rdistr)
    stdrobot <- allstds(rdistr)

## ## plot obzervations
## df <- data.frame(x=2:(maxtrials+1), y1=obs)
##     maxy <- N
##     miny <- 1
## g <- ggplot() + theme_classic() 
## g <- g + geom_point(data=df, aes(x,y1), colour='black') +
##     geom_line(data=df, aes(x,y1), colour='black') +
##     xlim(1,maxtrials+1) + ylim(miny,maxy) +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='observations',
##          title=paste0('participant ', participant))
## pdfname <- paste0(plotsdir,'observ_',participant,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()
## ##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)

## ## plot the overlap
## df <- data.frame(x=1:(maxtrials+1), y=overlap)
## g <- ggplot() + theme_classic() 
## g <- g + #geom_point(data=df, aes(x,y), colour=mygreen) +
##     geom_line(data=df, aes(x,y), colour=mygreen) +
##     xlim(1,maxtrials+1) + # ylim(0,1) +
##     theme_classic() +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='overlap',
##          title=paste0('participant ', participant,', stub ',stubbornness))
## pdfname <- paste0(plotsdir,'overlap_',participant,'-stub_',stubbornness,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()

    ## plots
if(graphs==TRUE){
    slabel <- substring(label,1,1)
    if(slabel!='' & slabel!='_'){label <- paste0('_',label)}
    pdfname <- paste0(plotsdir,'p',participant,label,'.pdf')
    pdf(pdfname, width = 148*mmtoin, height=148*0.6*mmtoin)
    
    ## plot sequence of means
    df <- data.frame(x=1:(maxtrials+1), y1=meanperson, y2=meanrobot, y3=c(NA,obs))
    maxy <- N #max(df$y1, df$y2)
    miny <- 1 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (maxtrials+1)/20/300*271
    personwidth <- (maxtrials+1)/20/300*223
    g <- ggplot() + theme_classic() 
    g <- g + annotation_raster(person, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
    g <- g + geom_line(data=df, aes(x,y3), colour='black', alpha=0.25)
    g <- g + #geom_point(data=df, aes(x,y1), colour=myred, alpha=0.33) +
        geom_line(data=df, aes(x,y1), colour=myred, alpha=0.75) +
                                        #geom_point(data=df, aes(x,y2), colour=myblue, alpha=0.33) +
        geom_line(data=df, aes(x,y2), colour=myblue, alpha=0.75) +
        xlim(1,maxtrials+1) + ylim(miny,maxy) +
        theme(aspect.ratio=0.5) +
        labs(x='observation',y='mean',
             title=paste0('participant #', participant,', means (black: plinko outcomes)'))
    g <- g + geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                           ymax=maxy,ymin = maxy-iconheight),
                       color=NA, fill=myred, alpha=0.75, stat='identity', position='identity') +
        geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                      ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                  color=NA, fill=myblue, alpha=0.75)
    print(g)

    ## plot sequence of stds
    df <- data.frame(x=1:(maxtrials+1), y1=stdperson, y2=stdrobot)
    maxy <- maxstd #max(df$y1, df$y2)
    miny <- 0 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (maxtrials+1)/20/300*271
    personwidth <- (maxtrials+1)/20/300*223
    g <- ggplot() + theme_classic() 
    g <- g + annotation_raster(person, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
    g <- g + #geom_point(data=df, aes(x,y1), colour=myred, alpha=0.33) +
        geom_line(data=df, aes(x,y1), colour=myred, alpha=0.75) +
                                        #geom_point(data=df, aes(x,y2), colour=myblue, alpha=0.33) +
        geom_line(data=df, aes(x,y2), colour=myblue, alpha=0.75) +
        xlim(1,maxtrials+1) + ylim(miny,maxy) +
        theme(aspect.ratio=0.5) +
        labs(x='observation',y='std',
         title=paste0('participant #', participant,', st. deviations'))
g <- g + geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.75, stat='identity', position='identity') +
        geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.75)
    print(g)

        ## plot h function
    nobs <- 0:maxtrials
    probd <- function(s,m,maxtrials,params){
	if(s>=m | m<0 | s<0){return(NA)}
	probchangepoint(s,m,maxtrials,params)
    }
    pmatrix <- sapply(nobs,function(i){sapply(nobs,function(j){probd(i,j, maxtrials,params)})})
    ##    png(paste0(plotsdir,label,'_optimh.png'))
    image2D(pmatrix,x=nobs,y=nobs,xlab='m',ylab='v',zlim=c(0,1),
            col=barpalettepos(25),NAcol=mygreen)

    title(paste0("part. #", participant, ", stub. = ", sprintf("%.3g",exp(params[1])),", changepoint:"))


    ## plot the relative entropy and Jansen-Shannon
    cols <- c('rel. entropy'=mygreen,'Jansen-Shannon'=myyellow)
    maxjs <- max(c(jsseq))
df <- data.frame(x=1:(maxtrials+1), y1=klseq, y2=jsseq)
g <- ggplot() + theme_classic() 
g <- g + #geom_point(data=df, aes(x,y), colour=myyellow) +
    geom_line(data=df, aes(x,y1, colour='rel. entropy'), alpha=0.75)
g <- g + #geom_point(data=df, aes(x,y), colour=myyellow) +
    geom_line(data=df, aes(x,y2, colour='Jansen-Shannon'), alpha=0.75)
g <- g + scale_color_manual(values=cols) +
    xlim(1,maxtrials+1) + ylim(0,maxjs) +
    theme(aspect.ratio=0.5, legend.title=element_blank(),
                   legend.background=element_blank(),
                   legend.justification=c(0,1),
                   legend.position=c(0,1),
) +
    labs(x='observation',y='discrepancy',
         title=paste0('participant #', participant), ', discrepancies')
    print(g)

        ## plot graphical discrepancy
    pmatrix <- rdistr-pdistr
    dims <- dim(pmatrix)
    ssqrt <- function(x){sign(x)*sqrt(abs(x))}
##    ssqrt <- function(x){sign(x)*(abs(x)^(1/3))}
    
    image2D(ssqrt(pmatrix),y=1:(dims[2]),x=1:(dims[1]),ylab='slot',xlab='observation',zlim=c(-1,1), col=barpalette(51))
    title(paste0("participant #", participant, ", sqrt(robot-participant)"))



          ## plot histograms for a range of trials
          if(!is.null(trialstoshow)){rangehist <- trialstoshow
              maxhist <- max(c(pdistr,rdistr))
              #pdf(pdfname,width = 148*mmtoin, height = 148*0.6*mmtoin)
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
             title=paste0('participant #', participant,', observation ',i,', rel-entr. = ',signif(klseq[i],2),', JS = ',signif(jsseq[i],2)))
g <- g + geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
print(g)
}}
dev.off()}

    return(list(distr=pdistr,robotdistr=rdistr,rentropy=klseq,js=jsseq,meanperson=meanperson,meanrobot=meanrobot,stdperson=stdperson,stdrobot=stdrobot))
}

comparisonall <- function(participant,maxtrials=200,trialstoshow=c(1:4, 99:102, 197:200),savedir='./',label='',params=rep(0.5,4),maxes,initial.distr=NULL){
    obs <- observations(participant,maxtrials)
    maxstd <- maxes[1] # max for std graph
    maxent <- maxes[2] # max for entopies graph
    maxdis <- maxes[3] # max for discrepancies
    maxhist <- maxes[4] # max for prob. histograms

    ## sequence of frequency parameters of the changepoint-JD model. At
    ## every "reset" we set the first equal to the participant's initial
    ## distribution plus a value equal to 1/1000 of the minimum nonzero value
    ## ever assigned
    pdistr <- regularizedistr(participant,maxtrials)

    ## robot's reset distribution
    if(length(initial.distr)==0){## use this participant's
        initial.distr <- pdistr[1,]
    }
    else{## use another participant's or custom
        initial.distr <- regularizedistr(initial.distr,maxtrials)[1,]
    }
    
    rdistr <- robotdistribution(initial.distr,params,obs)

    ## calculate the overlap and relative entropy of the participant's distr.
    jsseq <- sapply(1:(maxtrials+1),function(i){jsd(rdistr[i,],pdistr[i,])})
    klseq <- sapply(1:(maxtrials+1),function(i){kld(rdistr[i,],pdistr[i,])})

    ## sequence of means and stds
    meanperson <- allmeans(pdistr)
    stdperson <- allstds(pdistr)
    meanrobot <- allmeans(rdistr)
    stdrobot <- allstds(rdistr)

## ## plot obzervations
## df <- data.frame(x=2:(maxtrials+1), y1=obs)
##     maxy <- N
##     miny <- 1
## g <- ggplot() + theme_classic() 
## g <- g + geom_point(data=df, aes(x,y1), colour='black') +
##     geom_line(data=df, aes(x,y1), colour='black') +
##     xlim(1,maxtrials+1) + ylim(miny,maxy) +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='observations',
##          title=paste0('participant ', participant))
## pdfname <- paste0(plotsdir,'observ_',participant,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()
## ##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)

## ## plot the overlap
## df <- data.frame(x=1:(maxtrials+1), y=overlap)
## g <- ggplot() + theme_classic() 
## g <- g + #geom_point(data=df, aes(x,y), colour=mygreen) +
##     geom_line(data=df, aes(x,y), colour=mygreen) +
##     xlim(1,maxtrials+1) + # ylim(0,1) +
##     theme_classic() +
##     theme(aspect.ratio=0.5) +
##     labs(x='trial',y='overlap',
##          title=paste0('participant ', participant,', stub ',stubbornness))
## pdfname <- paste0(plotsdir,'overlap_',participant,'-stub_',stubbornness,'.pdf')
## save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
## dev.off()

    ## plots
    slabel <- substring(label,1,1)
    if(slabel!='' & slabel!='_'){label <- paste0('_',label)}
    pdfname <- paste0(savedir,'p',participant,label,'.pdf')
    pdf(pdfname, width = 148*mmtoin, height=148*0.6*mmtoin)
    
    ## plot sequence of means
    df <- data.frame(x=1:(maxtrials+1), y1=meanperson, y2=meanrobot, y3=c(NA,obs))
    maxy <- N #max(df$y1, df$y2)
    miny <- 1 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (maxtrials+1)/20/300*271
    personwidth <- (maxtrials+1)/20/300*223
    g <- ggplot() + theme_classic() 
    g <- g + annotation_raster(person, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
    g <- g + geom_line(data=df, aes(x,y3), colour='black', alpha=0.25)
    g <- g + #geom_point(data=df, aes(x,y1), colour=myred, alpha=0.33) +
        geom_line(data=df, aes(x,y1), colour=myred, alpha=0.75) +
                                        #geom_point(data=df, aes(x,y2), colour=myblue, alpha=0.33) +
        geom_line(data=df, aes(x,y2), colour=myblue, alpha=0.75) +
        xlim(1,maxtrials+1) + ylim(miny,maxy) +
        theme(aspect.ratio=0.5) +
        labs(x='observation',y='mean',
             title=paste0('participant #', participant,', means (black: plinko outcomes)'))
    g <- g + geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                           ymax=maxy,ymin = maxy-iconheight),
                       color=NA, fill=myred, alpha=0.75, stat='identity', position='identity') +
        geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                      ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                  color=NA, fill=myblue, alpha=0.75)
    print(g)

    ## plot sequence of stds
    df <- data.frame(x=1:(maxtrials+1), y1=stdperson, y2=stdrobot)
    maxy <- maxstd #max(df$y1, df$y2)
    miny <- 0 #min(df$y1, df$y2)
    iconheight <- (maxy-miny)/10
    robotwidth <- (maxtrials+1)/20/300*271
    personwidth <- (maxtrials+1)/20/300*223
    g <- ggplot() + theme_classic() 
    g <- g + annotation_raster(person, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = maxtrials+1, xmin = maxtrials+1-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
    g <- g + #geom_point(data=df, aes(x,y1), colour=myred, alpha=0.33) +
        geom_line(data=df, aes(x,y1), colour=myred, alpha=0.75) +
                                        #geom_point(data=df, aes(x,y2), colour=myblue, alpha=0.33) +
        geom_line(data=df, aes(x,y2), colour=myblue, alpha=0.75) +
        xlim(1,maxtrials+1) + ylim(miny,maxy) +
        theme(aspect.ratio=0.5) +
        labs(x='observation',y='std',
         title=paste0('participant #', participant,', st. deviations'))
g <- g + geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.75, stat='identity', position='identity') +
        geom_rect(aes(xmax=maxtrials+1-robotwidth,xmin=maxtrials+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.75)
    print(g)

        ## plot h function
    nobs <- 0:maxtrials
    probd <- function(s,m,maxtrials,params){
	if(s>=m | m<0 | s<0){return(NA)}
	probchangepoint(s,m,maxtrials,params)
    }
    pmatrix <- sapply(nobs,function(i){sapply(nobs,function(j){probd(i,j, maxtrials,params)})})
    ##    png(paste0(plotsdir,label,'_optimh.png'))
    image2D(pmatrix,x=nobs,y=nobs,xlab='m',ylab='v',zlim=c(0,1),
            col=barpalettepos(25),NAcol=myyellow)

    title(paste0("h for part. #", participant, ", stubbornness = ", sprintf("%.3g",exp(params[1]))))


    ## plot the relative entropy and Jansen-Shannon
    cols <- c('rel. entropy'=mygreen,'Jansen-Shannon'=myyellow)
df <- data.frame(x=1:(maxtrials+1), y1=klseq, y2=jsseq)
g <- ggplot() + theme_classic() 
g <- g + #geom_point(data=df, aes(x,y), colour=myyellow) +
    geom_line(data=df, aes(x,y1, colour='rel. entropy'), alpha=0.75)
g <- g + #geom_point(data=df, aes(x,y), colour=myyellow) +
    geom_line(data=df, aes(x,y2, colour='Jansen-Shannon'), alpha=0.75)
g <- g + scale_color_manual(values=cols) +
    xlim(1,maxtrials+1) + ylim(0,maxent) +
    theme(aspect.ratio=0.5, legend.title=element_blank(),
                   legend.background=element_blank(),
                   legend.justification=c(0,1),
                   legend.position=c(0,1),
) +
    labs(x='observation',y='discrepancy',
         title=paste0('participant #', participant), ', discrepancies')
    print(g)

        ## plot graphical discrepancy
    pmatrix <- rdistr-pdistr
    dims <- dim(pmatrix)
##    ssqrt <- function(x){sign(x)*sqrt(abs(x))}
    ssqrt <- function(x){x}
##    ssqrt <- function(x){sign(x)*(abs(x)^(1/3))}
    
    image2D(ssqrt(pmatrix),y=1:(dims[2]),x=1:(dims[1]),ylab='slot',xlab='observation',zlim=c(-maxdis,+maxdis), col=barpalette(51))
    title(paste0("participant #", participant, ", sqrt(robot-participant)"))



          ## plot histograms for a range of trials
          if(!is.null(trialstoshow)){rangehist <- trialstoshow
              ##maxhist <- max(c(pdistr,rdistr))
              #pdf(pdfname,width = 148*mmtoin, height = 148*0.6*mmtoin)
              for(j in 1:length(rangehist)){
                  i <- rangehist[j]
                  df <- data.frame(x=rep(1:N, 2),
                                   who=rep(c('person','robot'), each=N),
                                   y=c(pdistr[i,], rdistr[i,]))
                  maxy <- maxhist
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
             title=paste0('participant #', participant,', observation ',i,', rel-entr. = ',signif(klseq[i],2),', JS = ',signif(jsseq[i],2)))
g <- g + geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
print(g)
}
dev.off()}

    return(list(distr=pdistr,robotdistr=rdistr,rentropy=klseq,js=jsseq,meanperson=meanperson,meanrobot=meanrobot,stdperson=stdperson,stdrobot=stdrobot))
}



### function to calculate discrepancy on cubic grid
### time: ca 22000 s for grid of side 25
arraydiscrepancy <- function(participant,border=1e-6,gridpoints=11,label='',maxtrials=200,stubbornness=0.01){
    n <- maxtrials
    nparams <- 3 ##((nregions+1)^2+nregions+1)/2
    obs <- observations(participant,n)
    tdistr <- distribution(participant,n)
    ## we modify the participant's distributions adding plus a value equal
    ## to 1/100 of the minimum nonzero value ever assigned (to avoid zero
    ## probabilities in the robot)
    mini <- tdistr + min(c(tdistr[tdistr>0]))/1000
    pdistr <- t(sapply(1:(n+1),function(i){mini[i,]/sum(mini[i,])}))

    vpoints <- seq(border,1-border,length.out=gridpoints)
    tarray <- sapply(vpoints,function(k){sapply(vpoints,function(j){sapply(vpoints,function(i){
        discrepancyfast(logit(c(k, ## constant, 3rd array dim.
                            i, ## m coeff., columns
                            j)), ## s coeff., rows
                    pdistr,obs)})})})
    dim(tarray) <- rep(gridpoints,3)

    mint <- min(tarray[])
    maxt <- max(tarray[])

    pdf(paste0(plotsdir,label,'_discrepancy.pdf'))
    for(i in 1:gridpoints){
        image2D(tarray[,,i],zlim=c(mint,maxt),col=barpalettepos(50),
                NAcol=myred,
                x=vpoints,
                y=vpoints,
                xlab='invlogit(m coeff)',ylab='invlogit(s coeff)')
                title(paste0('intercept = ',logit(vpoints[i])))}
    dev.off()
    tarray}
## to know which dimension is which, check this:
## a <- sapply(1:2,function(k){sapply(1:3,function(j){sapply(1:4,function(i){i})})});dim(a) <- c(4,3,2);a
## , , 1
##      [,1] [,2] [,3]
## [1,]    1    1    1
## [2,]    2    2    2
## [3,]    3    3    3
## [4,]    4    4    4
##
## , , 2
##      [,1] [,2] [,3]
## [1,]    1    1    1
## [2,]    2    2    2
## [3,]    3    3    3
## [4,]    4    4    4

## testarray <- function(s,m,border=1e-6,gridpoints=11,label=''){
##     nparams <- 3 ##((nregions+1)^2+nregions+1)/2

##     vpoints <- seq(border,1-border,length.out=gridpoints)
##     tarray <- sapply(vpoints,function(k){sapply(vpoints,function(j){sapply(vpoints,function(i){
##         probchangepoint(s,m,200,1,logit(c(k, ## constant, 3rd array dim.
##                             i, ## m coeff., columns
##                             j)))})})})
##     dim(tarray) <- rep(gridpoints,3)

##     mint <- min(tarray[])
##     maxt <- max(tarray[])

##     pdf(paste0(plotsdir,label,'_discrepancy.pdf'))
##     for(i in 1:gridpoints){
##         image2D(tarray[,,i],zlim=c(0,1),col=barpalettepos(50),
##                 NAcol=myred,
##                 x=vpoints,
##                 y=vpoints,
##                 xlab='invlogit(m coeff)',ylab='invlogit(s coeff)')
##                 title(paste0('intercept = ',logit(vpoints[i])))}
##     dev.off()
##     tarray}
##  a<-testarray(100,100,1e-6,3,'test')

summaryparticipants <- function(participants=(1:40),maxtrials=200,label='',seed=999){
    slabel <- substring(label,1,1)
    if(slabel!='' & slabel!='_'){label <- paste0('_',label)}
    
    parlist <- list()
    reslist <- list()
    for(i in participants){
        message('participant ',i)
        optres <- reducediscrepancy(i,maxtrials,3,seed)
        warnlabel <- ''
        if(optres$convergence > 0){warnlabel='_NOTCONVERGED'}
        
        compres <- comparison(i,200,label=paste0(warnlabel,label),params=optres$par)
        
        saveRDS(list(optres=optres,compres=compres),paste0(plotsdir,'summary_p',i,warnlabel,label,'.rds'))
        message(' ')
    }
    message('finished')
}

generategraphsummary <- function(participants,dir,savedir,label){
    mat <- matrix(NA,length(participants),5)
    maxes <- matrix(-Inf,length(participants),4)
    j <- 0
    for(i in participants){j <- j+1
        filename <- paste0(dir,'summary_p',i,'_',label,'.rds')
        if(!file.exists(filename)){
        filename <- paste0(dir,'summary_p',i,'_NOTCONVERGED_',label,'.rds')}
        if(file.exists(filename)){
            mat[j,] <- c(i,readRDS(filename)$optres$par)

            summary <- c(i,readRDS(filename)$compres)
            maxes <- c(
                max(c(maxes[1], summary$stdperson,summary$stdrobot)),
                max(c(maxes[2], ##if(length(c(summary$rentropy[summary$rentropy==Inf]))>0){-Inf}else{summary$rentropy},
                      summary$js)),
                max(c(maxes[3], max(abs(summary$robotdistr-summary$distr)))),
                max(c(maxes[4], summary$distr,summary$robotdistr))
            )
        }}

    write.table(mat,paste0(plotsdir,'points_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
    write.table(maxes,paste0(plotsdir,'maxima_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')

    cat('maxima: ',maxes,'\n')
    
    for(i in participants){j <- j+1
        filename <- paste0(dir,'summary_p',i,'_',label,'.rds')
        warnlabel <- ''
        if(!file.exists(filename)){
        filename <- paste0(dir,'summary_p',i,'_NOTCONVERGED_',label,'.rds')}
        if(file.exists(filename)){

            optp <- c(i,readRDS(filename)$optres)
            if(optp$convergence > 0){warnlabel='_NOTCONVERGED'}
            compres <- comparisonall(i,200,savedir=savedir,label=paste0(warnlabel,label),params=optp$par,maxes=maxes)
        }}
}

saveparams <- function(participants,dir,label){
    setwd(dir)
    mat <- matrix(NA,length(participants),4)
    j <- 0
    for(i in participants){j <- j+1
        filename <- paste0('summary_p',i,'_',label,'.rds')
        if(!file.exists(filename)){
        filename <- paste0('summary_p',i,'_NOTCONVERGED_',label,'.rds')}
        if(file.exists(filename)){
        mat[j,] <- c(i,readRDS(filename)$optres$par)
        }}
    write.table(mat,paste0('points_',label,'.csv'),sep=',',row.names=F,col.names=F,na='Infinity')
}
