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
plotsdir <- './comparisons1/'

## load all data
dpath = "./data/"
fs = dir(path = dpath,pattern = "barData")
d <- read.csv(paste0(dpath,fs[1]))
for (f in fs[2:40]) {
dprime <- read.csv(paste0(dpath, f))
d <- rbind(d,dprime)
rm(dprime)
}
d <- d[,1:(ncol(d)-2)]
person <- readPNG('face3.png')
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

## sequence of robot distributions
robotdistribution <- function(priorfrequencies,stubbornness,obs){
    n <- length(obs)
    ## Johnson-Dirichlet model:
    ## stubbornness parameter and ensuing sequence
    L <- stubbornness + (0:(n+1))
    ## initialize
    ldistr <- matrix(NA, n+1, N)
    ldistr[1,] <- priorfrequencies
    ## Update the frequency parameters of the model from the observations
    ## as explained in plinkinetti180414
    for(i in 1:n){
        ldistr[1+i,] <- (L[i] * ldistr[i,] + replace(slots, obs[i], 1))/L[i+1]
    }
ldistr}

## sequences of means and stds
allmeans <- function(distrib){distrib %*% (1:N)}
allstds <- function(distrib){sqrt(distrib %*% (1:N)^2 - allmeans(distrib)^2)}


## Define a function that plots:
## - sequences of overlap, relative entropies, means, stds
## - histograms for a selected set of trials
## and that outputs the final distributions, rel-entropies, overlaps, means, stds
comparison <- function(participant,maxtrials=200,trialstoshow=c(1:10, 95:105, 190:200),stubbornness=0.1){
## number of trials to consider
    n <- maxtrials
    
    obs <- observations(participant,n)

    distr <- distribution(participant,n)

    ## sequence of frequency parameters of the JD model
    ## we set the first equal to the participant's initial distribution
    ## plus a value equal to 1/10 of the minimum nonzero value ever assigned
    mini <- min(c(distr[distr>0]))
    temp <- distr[1,] + mini/10
    temp <- temp/sum(temp)
    ldistr <- robotdistribution(temp,stubbornness,obs)

    ## calculate the overlap and relative entropy of the participant's distr.
    temp <-  distr * log(distr/ldistr)
    temp[is.nan(temp)] <- 0
    rentropy <- apply(temp,1,sum)
    overlap <- diag(ldistr %*% t(distr))

    ## sequence of means and stds
    meanperson <- allmeans(distr)
    stdperson <- allstds(distr)
    meanrobot <- allmeans(ldistr)
    stdrobot <- allstds(ldistr)

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

## plot the overlap
df <- data.frame(x=1:(n+1), y=overlap)
g <- ggplot() + theme_classic() 
g <- g + #geom_point(data=df, aes(x,y), colour=mygreen) +
    geom_line(data=df, aes(x,y), colour=mygreen) +
    xlim(1,n+1) + # ylim(0,1) +
    theme_classic() +
    theme(aspect.ratio=0.5) +
    labs(x='trial',y='overlap',
         title=paste0('participant ', participant,', stub ',stubbornness))
pdfname <- paste0(plotsdir,'overlap_',participant,'-stub_',stubbornness,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
dev.off()

## plot the relative entropy
df <- data.frame(x=1:(n+1), y=rentropy)
g <- ggplot() + theme_classic() 
g <- g + #geom_point(data=df, aes(x,y), colour=myyellow) +
    geom_line(data=df, aes(x,y), colour=myyellow) +
    xlim(1,n+1) +
    theme(aspect.ratio=0.5) +
    labs(x='trial',y='relative entropy',
         title=paste0('participant ', participant,', stub ',stubbornness))
pdfname <- paste0(plotsdir,'rentropy_',participant,'-stub_',stubbornness,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
dev.off()
##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)

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
         title=paste0('participant ', participant,', stub ',stubbornness))
g <- g + geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
pdfname <- paste0(plotsdir,'means_',participant,'-stub_',stubbornness,'.pdf')
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
         title=paste0('participant ', participant,', stub ',stubbornness))
g <- g + geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=n+1-robotwidth,xmin=n+1-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
pdfname <- paste0(plotsdir,'stds_',participant,'-stub_',stubbornness,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height = 148*0.6, units='mm', dpi = 300)
dev.off()
##ggsave(pdfname, width = 148, height = 148*0.6, units='mm', dpi = 300)


## plot histograms for a range of trials
if(!is.null(trialstoshow)){rangehist <- trialstoshow
maxhist <- max(c(distr,ldistr))
pdfname <- paste0(plotsdir,'histogram_',rangehist[1],'-',rangehist[length(rangehist)],'_',participant,'-stub_',stubbornness,'.pdf')
pdf(pdfname,width = 148*mmtoin, height = 148*0.6*mmtoin)
for(j in 1:length(rangehist)){
    i <- rangehist[j]
    df <- data.frame(x=rep(1:N, 2),
                     who=rep(c('person','robot'), each=N),
                     y=c(distr[i,], ldistr[i,]))
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
             title=paste0('participant ', participant,', stub ',stubbornness,', trial ',i,', overlap ',signif(overlap[i],2),', relentr. ',signif(rentropy[i],2)))
g <- g + geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
print(g)
}
dev.off()}

    return(list(distr=distr,robotdistr=ldistr,rentropy=rentropy,overlap=overlap,meanperson=meanperson,meanrobot=meanrobot,stdperson=stdperson,stdrobot=stdrobot))
}


comparerobots <- function(index,obs,changepoint,priorfrequencies=rep(1/N,N),stubbornness=0.1){
    n <- length(obs)
    distr1 <- robotdistribution(priorfrequencies,stubbornness,obs)
    distr2 <- robotdistribution(distr1[changepoint-1,],stubbornness,obs[changepoint:n])
    pdfname <- paste0(plotsdir,'compare_robots_',index,'-stub_',stubbornness,'.pdf')
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
    overlap <- diag(distr2 %*% t(distr1t))
return(c(rentropy[length(rentropy)],rentropy2[length(rentropy2)], overlap[length(overlap)]))
}
