
model.names <- c('null', 'source', 'destination', 'source-destination_p=1')

data.dir <- 'output'

result <- lapply(model.names, function(model.name) {
    props <- c(0.05, 0.5, 0.95)
    ##
    load(file.path(data.dir,
                   paste0(model.name, '.RData')))
    ##
    if ('spatial.ix12' %in% names(result[[1]][[1]])) {
        spatial.ix.label <- 'spatial.ix12'
    } else {
        spatial.ix.label <- 'spatial.ix'
    }
    res <- Filter(function(x) x[[1]][spatial.ix.label]^(1/3) <= 1.6, result)
    do.call(rbind, lapply(res, function(x) {
        n <- x[[1]]['n']
        spatial.ix <- x[[1]][spatial.ix.label]
        times <- do.call(rbind, lapply(x[[2]], function(y) {
            approx(y[,1]/n^2, 0:(nrow(y) - 1),
                   xout=props)$y
        }))
        first.second.period.ratio <-
            (times[,2] - times[,1])/(times[,3] - times[,2])
        v <- c(spatial.ix,
               median(times[,1]),
               quantile(times[,1], probs=c(0.25, 0.75)),
               median(times[,2]),
               quantile(times[,2], probs=c(0.25, 0.75)),
               median(first.second.period.ratio),
               quantile(first.second.period.ratio, probs=c(0.25, 0.75)))
        ##
        names(v) <- c('spatial.ix',
                      'median.0.05',
                      paste(paste0('q', c('25', '75')),
                            '0.05', sep='.'),
                      'median.0.5',
                      paste(paste0('q', c('25', '75')),
                            '0.5', sep='.'),
                      'period.ratio.median',
                      paste('period.ratio',
                            paste0('q', c('25', '75')),
                            props[2], sep='.'))
        return(v)
    }))
})
names(result) <- model.names

result.grate <- lapply(model.names, function(model.name) {
    ##
    load(file.path(data.dir,
                   paste0(model.name, '.RData')))
    ##
    if ('spatial.ix12' %in% names(result[[1]][[1]])) {
        spatial.ix.label <- 'spatial.ix12'
    } else {
        spatial.ix.label <- 'spatial.ix'
    }
    res <- Filter(function(x) x[[1]][spatial.ix.label]^(1/3) <= 1.6, result)
    m <- do.call(rbind, lapply(res, function(x) {
        n <- x[[1]]['n']
        ##
        spatial.ix <- x[[1]][spatial.ix.label]
        growth.rates <- unlist(lapply(x[[2]], function(y) {
            prop <- y[,'noccupied']/n^2
            ix.in.range <- (0.01 <= prop) & (prop <= 0.5)
            na.omit(diff(log(prop))[ix.in.range])
        }))
        p.dens <- density(growth.rates, bw='SJ')
        growth.rate.peak <- p.dens$x[which.max(p.dens$y)]
        ##
        c(spatial.ix=spatial.ix[[1]],
          growth.rate=growth.rate.peak)
    }))
    return(m)
})
names(result.grate) <- model.names

result.period <- lapply(model.names, function(model.name) {
    props <- c(0.05, 0.95)
    ##
    load(file.path(data.dir,
                   paste0(model.name, '.RData')))
    ##
    if ('spatial.ix12' %in% names(result[[1]][[1]])) {
        spatial.ix.label <- 'spatial.ix12'
    } else {
        spatial.ix.label <- 'spatial.ix'
    }
    res <- Filter(function(x) x[[1]][spatial.ix.label]^(1/3) <= 1.6, result)
    m <- do.call(rbind, lapply(res, function(x) {
        n <- x[[1]]['n']
        spatial.ix <- x[[1]][spatial.ix.label]
        times <- do.call(rbind, lapply(x[[2]], function(y) {
            approx(y[,1]/n^2, 0:(nrow(y) - 1),
                   xout=props)$y
        }))
        periods <- times[,2] - times[,1]
        v <- c(spatial.ix,
               median(periods),
               quantile(periods, probs=c(0.25, 0.75)))
        ##
        names(v) <- c('spatial.ix',
                      'period.median',
                      paste('period', paste0('q', c('25', '75')), sep='.'))
        return(v)
    }))
    return(m)
})
names(result.period) <- model.names

result.period.closing <- lapply(model.names, function(model.name) {
    props <- c(0.95, 1)
    ##
    load(file.path(data.dir,
                   paste0(model.name, '.RData')))
    ##
    if ('spatial.ix12' %in% names(result[[1]][[1]])) {
        spatial.ix.label <- 'spatial.ix12'
    } else {
        spatial.ix.label <- 'spatial.ix'
    }
    res <- Filter(function(x) x[[1]][spatial.ix.label]^(1/3) <= 1.6, result)
    m <- do.call(rbind, lapply(res, function(x) {
        n <- x[[1]]['n']
        spatial.ix <- x[[1]][spatial.ix.label]
        times <- do.call(rbind, lapply(x[[2]], function(y) {
            approx(y[,1]/n^2, 0:(nrow(y) - 1),
                   xout=props)$y
        }))
        periods <- times[,2] - times[,1]
        v <- c(spatial.ix,
               median(periods),
               quantile(periods, probs=c(0.25, 0.75)))
        ##
        names(v) <- c('spatial.ix',
                      'period.median',
                      paste('period', paste0('q', c('25', '75')), sep='.'))
        return(v)
    }))
    return(m)
})
names(result.period.closing) <- model.names

plot.target <- 'plot'

(function() {
    prop.percent <- 5
    prop <- format(prop.percent/100, digits=2)
    xlim <- c(1, 1.6)
    ylim <- c(0, 25)
    pdf(file.path(plot.target,
                  paste('establishment-phase-period.pdf',
                        sep='-')),
        height=5, width=5,
        useDingbats=FALSE)
    on.exit(dev.off())
    plot(NA, xlim=xlim, ylim=ylim,
         xlab=expression('Spatial index'^{1/3}),
         ylab=paste('Establishment phase period'),
         xaxt='n')
    axis(1, at=c(1, 1.2, 1.4, 1.6))
    val.main <- paste('median', prop, sep='.')
    ##
    for (p in list(list(type='null', col=gray(0.5), pch=1),
                   list(type='destination', col=gray(0.24), pch=17),
                   list(type='source', col=hsv(0,0.9,0.8), pch=15),
                   list(type='source-destination_p=1', col=hsv(2/3,0.54,0.72), pch=16)
                   )) {
        basecol.rgb <- col2rgb(p$col)
        basecol.hsv <- rgb2hsv(basecol.rgb[1],
                               basecol.rgb[2],
                               basecol.rgb[3])
        barcol <- hsv(basecol.hsv[1],
                      basecol.hsv[2]*0.5,
                      min(basecol.hsv[3]*1.5, 1))
        if (p$type == 'null') {
            n <- 1
            sp.ix <- rep(1, n)
        } else {
            n <- nrow(result[[p$type]])
            sp.ix <- result[[p$type]][,1]
        }
        matlines(matrix(rep(sp.ix^(1/3), each=2),
                        nrow=2, ncol=n),
                 rbind(result[[p$type]][1:n,paste('q25', prop, sep='.')],
                       result[[p$type]][1:n,paste('q75', prop, sep='.')]),
                 col=barcol,
                 lty=1,
                 lwd=2)
        matpoints(sp.ix^(1/3),
                  result[[p$type]][1:n,val.main],
                  pch=p$pch,
                  col=p$col)
    }
})()

(function() {
    xlim <- c(1, 1.6)
    ylim <- c(0.3, 0.7)
    pdf(file.path(plot.target, 'asymptotic-growth-rate.pdf'),
        height=5, width=5,
        useDingbats=FALSE)
    on.exit(dev.off())
    plot(NA, xlim=xlim, ylim=ylim,
         xlab=expression('Spatial index'^{1/3}),
         ylab='Asymptotic growth rate',
         xaxt='n')
    axis(1, at=c(1, 1.2, 1.4, 1.6))
    R <- 0.01
    g <- 1
    s <- 1/3
    lines(xlim, rep((2*pi*R*g^2)^(1/3), 2),
          col='darkgray',
          lwd=2,
          lty=1)
    ##
    x.theor <- seq(1, 1.6^3, length.out=101)
    y.theor <- (function() {
        lambda <- (2*pi*R*g^2*x.theor)^(1/3)
        z <- (sqrt(54-lambda^3*s^6) + sqrt(54))^(1/3)
        lambda/sqrt(6)*(z + lambda*s^2/z)
    })()
    lines(x=x.theor^(1/3), y=y.theor, lty=1, lwd=2,
          col=hsv(2/3,0.2,1))
    ##
    for (p in list(list(type='null', col=gray(0.5), pch=1),
                   list(type='destination', col=gray(0.24), pch=17),
                   list(type='source', col=hsv(0,0.9,0.8), pch=15),
                   list(type='source-destination_p=1', col=hsv(2/3,0.54,0.72), pch=16)
                   )) {
        if (p$type == 'null') {
            n <- 1
            sp.ix <- rep(1, n)
        } else {
            n <- nrow(result[[p$type]])
            sp.ix <- result[[p$type]][,1]
        }
        matpoints(sp.ix^(1/3),
                  result.grate[[p$type]][1:n,'growth.rate'],
                  pch=p$pch,
                  col=p$col)
    }
})()


(function() {
    xlim <- c(1, 1.6)
    ylim <- c(0, 100)
    pdf(file.path(plot.target, 'expansion-phase-period.pdf'),
        height=5, width=5,
        useDingbats=FALSE)
    plot(NA, xlim=xlim, ylim=ylim,
         xlab=expression('Spatial index'^{1/3}),
         ylab='Expansion phase period',
         xaxt='n')
    axis(1, at=c(1, 1.2, 1.4, 1.6))
    ##
    for (p in list(list(type='null', col=gray(0.5), pch=1),
                   list(type='destination', col=gray(0.24), pch=17),
                   list(type='source', col=hsv(0,0.9,0.8), pch=15),
                   list(type='source-destination_p=1', col=hsv(2/3,0.54,0.72), pch=16)
                   )) {
        if (p$type == 'null') {
            n <- 1
            sp.ix <- rep(1, n)
        } else {
            n <- nrow(result[[p$type]])
            sp.ix <- result[[p$type]][,1]
        }
        basecol.rgb <- col2rgb(p$col)
        basecol.hsv <- rgb2hsv(basecol.rgb[1],
                               basecol.rgb[2],
                               basecol.rgb[3])
        barcol <- hsv(basecol.hsv[1],
                      basecol.hsv[2]*0.5,
                      min(basecol.hsv[3]*1.5, 1))
        matlines(matrix(rep(sp.ix^(1/3), each=2),
                        nrow=2, ncol=n),
                 rbind(result.period[[p$type]][1:n,paste('period', 'q25', sep='.')],
                       result.period[[p$type]][1:n,paste('period', 'q75', sep='.')]),
                 col=barcol,
                 lty=1,
                 lwd=2)
        matpoints(sp.ix^(1/3),
                  result.period[[p$type]][1:n,'period.median'],
                  pch=p$pch,
                  col=p$col)
    }
    dev.off()
})()

(function() {
    xlim <- c(1, 1.6)
    ylim <- c(0, 210)
    pdf(file.path(plot.target, 'naturalization-phase-period.pdf'),
        height=5, width=5,
        useDingbats=FALSE)
    plot(NA, xlim=xlim, ylim=ylim,
         xlab=expression('Spatial index'^{1/3}),
         ylab='Naturalization phase period',
         xaxt='n')
    axis(1, at=c(1, 1.2, 1.4, 1.6))
    ##
    y <- result.period.closing
    for (p in list(list(type='null', col=gray(0.5), pch=1),
                   list(type='destination', col=gray(0.24), pch=17),
                   list(type='source', col=hsv(0,0.9,0.8), pch=15),
                   list(type='source-destination_p=1', col=hsv(2/3,0.54,0.72), pch=16)
                   )) {
        if (p$type == 'null') {
            n <- 1
            sp.ix <- rep(1, n)
        } else {
            n <- nrow(result[[p$type]])
            sp.ix <- result[[p$type]][,1]
        }
        basecol.rgb <- col2rgb(p$col)
        basecol.hsv <- rgb2hsv(basecol.rgb[1],
                               basecol.rgb[2],
                               basecol.rgb[3])
        barcol <- hsv(basecol.hsv[1],
                      basecol.hsv[2]*0.5,
                      min(basecol.hsv[3]*1.5, 1))
        matlines(matrix(rep(sp.ix^(1/3), each=2),
                        nrow=2, ncol=n),
                 rbind(y[[p$type]][1:n,paste('period', 'q25', sep='.')],
                       y[[p$type]][1:n,paste('period', 'q75', sep='.')]),
                 col=barcol,
                 lty=1,
                 lwd=2)
        matpoints(sp.ix^(1/3),
                  y[[p$type]][1:n,'period.median'],
                  pch=p$pch,
                  col=p$col)
    }
    dev.off()
})()
