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

## number of slots and sequence of zeros for later use
N <- 40
slots <- integer(N)

## number of trials to consider
n <- 200

## choose a participant
x <- 8


## sequence of observations for that participant
obs <- c(d[d$Participant==x & d$Slot.Number==1,'Ball.Position'])[1:n]

## sequence of distributions for that participant
## (consistency checked by for-loop)
## and find minimal non-zero value
distr <- matrix(d[d$Participant == x & d$Trial <= n+1, c('Participant.Slot.Estimate')], n+1, N, byrow=T)
mini <- min(c(distr[distr>0]))
##

## Johnson-Dirichlet model:
## parameter Lambda and ensuing sequence
L <- 0.1 + (0:(n+1))

## sequence of frequency parameters of the model
## we set the first equal to the participant's initial distribution
## plus a value equal to 1/10 of the minimum nonzero value ever assigned
ldistr <- matrix(NA, n+1, N)
temp <- distr[1,] + mini/10
ldistr[1,] <- temp/sum(temp)

##initialize relative entropy and overlap
rentropy <- rep(NA,n+1)
rentropy[1] <- 0
overlap <- rep(NA,n+1)
overlap[1] <- ldistr[1,] %*% distr[1,]

## Update the frequency parameters of the model from the observations
## as explained in plinkinetti180414
## and calculate the overlap and relative entropy of the participant's distr. 
for(i in 1:n){
    ldistr[1+i,] <- (L[i] * ldistr[i,] + replace(slots, obs[i], 1))/L[i+1]
    temp <-  distr[1+i,] * log(distr[i+1,]/ldistr[i+1,])
    temp[is.nan(temp)] <- 0
    rentropy[1+i] <- sum(temp)
    overlap[1+i] <- ldistr[1+i,] %*% distr[1+i,]
}

## plot the overlap
df <- data.frame(x=1:(n+1), y=overlap)
g <- ggplot() + theme_classic() 
g <- g + geom_point(data=df, aes(x,y), colour=mygreen) +
    geom_line(data=df, aes(x,y), colour=mygreen) +
    xlim(1,n+1) + # ylim(0,1) +
    theme_classic() +
    theme(aspect.ratio=0.5) +
    labs(x='trial',y='overlap',
         title=paste0('participant ', x))
pdfname <- paste0('overlap_',x,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height = 148*0.5, units='mm', dpi = 300)
dev.off()

## plot the relative entropy
df <- data.frame(x=1:(n+1), y=rentropy)
g <- ggplot() + theme_classic() 
g <- g + geom_point(data=df, aes(x,y), colour=myyellow) +
    geom_line(data=df, aes(x,y), colour=myyellow) +
    xlim(1,n+1) +
    theme(aspect.ratio=0.5) +
    labs(x='trial',y='relative entropy',
         title=paste0('participant ', x))
pdfname <- paste0('rentropy_',x,'.pdf')
save_plot(pdfname, g, base_width = 148, base_height = 148*0.5, units='mm', dpi = 300)
dev.off()
##ggsave(pdfname, width = 148, height = 148*0.5, units='mm', dpi = 300)


## plot histograms for a range of trials
rangehist <- 1:10
plots <- list()
maxhist <- max(c(distr,ldistr))
pdfname <- paste0('histogram_',rangehist[1],'-',rangehist[length(rangehist)],'_',x,'.pdf')
pdf(pdfname,width = 148*mmtoin, height = 148*0.5*mmtoin)
for(j in 1:length(rangehist)){
    i <- rangehist[j]
    df <- data.frame(x=rep(1:N, 2),
                     who=rep(c('person','robot'), each=N),
                     y=c(distr[i,], ldistr[i,]))
    maxy <- max(df$y)
    iconheight <- maxy/10
    robotwidth <- iconheight/maxy*N*0.5/300*271
    personwidth <- iconheight/maxy*N*0.5/300*223
    g <- ggplot() + theme_classic() 
g <- g + annotation_raster(person, xmax = N, xmin = N-robotwidth, 
                               ymax = maxy, ymin=maxy-iconheight, interpolate = T) +
        annotation_raster(robot, xmax = N, xmin = N-robotwidth, 
                          ymax = maxy*0.99-iconheight, ymin=maxy*0.99-2*iconheight, interpolate = T) 
    g <- g + geom_bar(data=df, aes(x=x,y=y,fill=who), alpha=0.5,stat='identity', position='identity') +
         scale_fill_manual(values=c(myred,myblue)) +
        scale_x_continuous(breaks=c(1,seq(5,N,length.out=8))) + 
        #xlim(1,N) + #ylim(0,maxhist) +
        theme(aspect.ratio=0.5,
              legend.title=element_blank(),
              legend.background=element_blank(),
              legend.justification=c(1,1),
              legend.position='none',
              ) +
        labs(x='slot',y='probability',
             title=paste0('participant ', x,', trial ',i,', overlap ',signif(overlap[i],3),', r. entr. ',signif(rentropy[i],3)))
g <- g + geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy,ymin = maxy-iconheight),
                          color=NA, fill=myred, alpha=0.5, stat='identity', position='identity') +
        geom_rect(aes(xmax=N-robotwidth,xmin=N-2*robotwidth,
                          ymax=maxy*0.99-iconheight,ymin = maxy*0.99-2*iconheight),
                          color=NA, fill=myblue, alpha=0.5)
print(g)
}
dev.off()


