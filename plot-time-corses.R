
library(parallel)

source('diskGrowthModel2.R')

if(exists('cl')) {
    stopCluster(cl)
    rm('cl')
}
cl <- makeCluster(detectCores(logical=FALSE))
invisible(clusterEvalQ(cl, source('diskGrowthModel2.R')))

n <- 2048
R <- 0.01
g <- 1
s <- 0
tend <- 300

ix.habitat <- 9
habitat.info <- (function() {
    load(file.path('output', paste0(ix.habitat, '.RData')), environment())
    list(h=h,
         h.beta=result[[1]]['h.beta'],
         h.pow=result[[1]]['h.pow'],
         spatial.ix=result[[1]]['spatial.ix'])
})()
image(habitat.info$h,useRaster=TRUE)
print(paste('(spatial factor)^(1/3):',
            0+(sum(habitat.info$h^2)*length(habitat.info$h))^(1/3)))

plot.dir <- 'plot'

png(file.path(plot.dir, 'landscape-for-fig2.png'),
    height=720, width=720,
    bg='transparent')
par(pty='s')
image(habitat.info$h,
      useRaster=TRUE,
      col=gray.colors(100),
      xaxt='n', yaxt='n', xlab='', ylab='',asp=1)
dev.off()


h.homogenous <- matrix(rep(1, 4), 2, 2)/4

for (type in c('hh', 'h0', '0h')) {
    if (substring(type, 1, 1) == 'h') {
        h.attraction <- habitat.info$h
    } else {
        h.attraction <- h.homogenous
    }
    if (substring(type, 2, 2) == 'h') {
        h.production <- habitat.info$h
    } else {
        h.production <- h.homogenous
    }
    ##
    system.time(
        res <- parLapply(cl, 1:100, function(ix, n, R, g, s, h.a, h.p, tend) {
            h.homogenous <- matrix(rep(1, 4), 2, 2)/4
            run(propagule.attraction=h.a,
                disperser.production=h.p,
                tm=1:tend,
                R=R,
                g=g,
                s=s,
                n=n,
                dt=1,
                dx=1,
                stop.when.complete=TRUE,
                show.progress=FALSE,
                output.populations=FALSE)$stats
        }, n=n, R=R, g=g, s=s,
        h.a=h.attraction, h.p=h.production,
        tend=tend)
    )
    save(res, file=paste('191118-time-course',
                         paste0('map', ix.habitat),
                         type,
                         'RData',
                         sep='.'))
}


for (type in c('00')) {
    if (substring(type, 1, 1) == 'h') {
        h.attraction <- habitat.info$h
    } else {
        h.attraction <- h.homogenous
    }
    if (substring(type, 2, 2) == 'h') {
        h.production <- habitat.info$h
    } else {
        h.production <- h.homogenous
    }
    ##
    system.time(
        res <- parLapply(cl, 1:100, function(ix, n, R, g, s, h.a, h.p, tend) {
            h.homogenous <- matrix(rep(1, 4), 2, 2)/4
            run(propagule.attraction=h.a,
                disperser.production=h.p,
                tm=1:tend,
                R=R,
                g=g,
                s=s,
                n=n,
                dt=1,
                dx=1,
                stop.when.complete=TRUE,
                show.progress=FALSE,
                output.populations=FALSE)$stats
        }, n=n, R=R, g=g, s=s,
        h.a=h.attraction, h.p=h.production,
        tend=tend)
    )
    save(res, file=paste('191118-time-course',
                         paste0('map', ix.habitat),
                         type,
                         'RData',
                         sep='.'))
}


z <- run(target.proportion=0.1,
                propagule.attraction=habitat.info$h,
                disperser.production=habitat.info$h,
                tm=1:10,
                R=R,
                g=g,
                s=s,
                n=n,
                dt=1,
                dx=1,
                stop.when.complete=TRUE,
                show.progress=FALSE,
    output.populations=TRUE)


model.settings <- c('h0', '0h', 'hh', '00')

## pdf('181105-time-course-draft.pdf',
##     height=7.5, width=11.25)
## layout(rbind(c(1,2,3),
##              c(4,5,6)))

## layout(cbind(1,2,3))
for (data.type in model.settings) {
    if (substring(data.type, 1, 1) == 'h') {
        h.attraction <- habitat.info$h
    } else {
        h.attraction <- h.homogenous
    }
    if (substring(data.type, 2, 2) == 'h') {
        h.production <- habitat.info$h
    } else {
        h.production <- h.homogenous
    }
    pops <- run(target.proportion=0.1,
                propagule.attraction=h.attraction,
                disperser.production=h.production,
                tm=1:tend,
                R=R,
                g=g,
                s=s,
                n=n,
                dt=1,
                dx=1,
                stop.when.complete=TRUE,
                show.progress=FALSE,
                output.populations=TRUE)
    ix.plot <- which.min((pops$stats[,'noccupied']/n^2 - 0.1)^2)
    n.plot <- n
    h <- h.production[rep(1:nrow(h.production), each=n/nrow(h.production)),
                      rep(1:ncol(h.production), each=n/ncol(h.production))]
    h <- h/sum(h)
    print(sum(h))
    print(sum(habitat.info$h))
    m <- drawOnMatrix(n, pops$populations[[ix.plot]])
    m[m == 0] <- NA
    m <- m*h[rep(1:nrow(h), each=n/nrow(h)),
             rep(1:ncol(h), each=n/ncol(h))]
    print(range(m, na.rm=TRUE))
    print(range(habitat.info$h, na.rm=TRUE))
    pdf(file.path(plot.dir, paste('191118-colony-distribution-example-10-percent-occupation',
              data.type,
              paste0('t', ix.plot),
              'pdf',
              sep='.')),
        height=5, width=5)
    par(pty='s')
    image(seq(0, n-1, length.out=ncol(habitat.info$h)),
          seq(0, n-1, length.out=ncol(habitat.info$h)),
          habitat.info$h,
          xlab='X',
          ylab='Y',
          col=gray.colors(100, start=1, end=0.75),
          asp=1,
          useRaster=TRUE)
    image(seq(0, n-1, length.out=n.plot),
          seq(0, n-1, length.out=n.plot),
          m,
          col=colorRampPalette(c(rgb(0.5, 0.5, 1), rgb(0.8, 0, 0)), space='Lab')(100),
          zlim=range(habitat.info$h)/4^2,
          useRaster=TRUE,
          add=TRUE)
    ## text(2000, 200, paste('t =', ix.plot),
    ##      adj=1,
    ##      cex=2)
    dev.off()
}
## layout(1)

line.cols <- c(type.hh=hsv(2/3,0.9,0.8),
               type.h0=gray(0.2),
               type.0h=hsv(0, 1, 0.8),
               type.00=gray(0.2))
for(data.type in model.settings) {
    col <- line.cols[paste('type', data.type, sep='.')]
    data.file <- paste(paste0('191118-time-course.map', ix.habitat),
                       data.type, 'RData',
                       sep='.')
    load(data.file)
    ys <- do.call(cbind, lapply(res, function(st) {
        v <- st[,'noccupied']/n^2
        c(v, rep(v[length(v)], tend - length(v)))
    }))
    zs <- do.call(cbind, lapply(res, function(st) {
        v <- st[,'noccupied']/n^2
        ## v <- c(v, rep(v[length(v)], tend - length(v)))
        approx(v, 1:length(v), xout=seq(0, 1, length.out=201))$y
    }))
    pdf(file.path(plot.dir, paste('191118-time-course',
                                  data.type,
                                  'pdf',
                                  sep='.')),
        height=5, width=5)
    plot(NA, xlim=c(1, 60), ylim=c(1/n^2, 1), log='',
         xlab='Time',
         ylab='Population size (covered proportion)')
    matlines(1:nrow(ys), ys, lty=1, col='gray')
    lines(c(rowMeans(zs), tend),
          c(seq(0, 1, length.out=201), 1),
          col=col, lwd=2)
    dev.off()
}
## layout(1)
## dev.off()
