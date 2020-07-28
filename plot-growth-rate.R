
model.names <- c('source-destination_p=0')

data.dir <- 'output'

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
            ## Adjusted from (0.01 <= prop) & (prop <= 0.5)
            ix.in.range <- (0.0 <= prop) & (prop <= 0.1)
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

n <- 2048
R <- 0.01
g <- 1
s <- 0
error.ix <- which(abs(result.grate[[1]][,"growth.rate"] - (2*pi*R*g^2*result.grate[[1]][,"spatial.ix"])^(1/3)) > 0.1)

source('diskGrowthModel2.R')

(function(){
    layout(rbind(c(1,2), c(3,4)))
    print(Map(function(ix) {
        load(file.path("190720-RData", "source-destination_p=0", paste(ix, "RData", sep=".")))
        seed <- result[[1]]['seed']
        set.seed(seed)
        h1.beta <- runif(1, -5, -2)
        h1.pow <- runif(1, 1, 5)
        h2.beta <- runif(1, -5, -2)
        h2.pow <- runif(1, 1, 5)
        h1 <- genHabitatMap(n, h1.beta, h1.pow)
        h2 <- genHabitatMap(n, h2.beta, h2.pow)
        ##
        ## image(h1, useRaster=TRUE)
        ## image(h2, UseRaster=TRUE)
        c(X=sum(h1*h2)*n^2, X1=sum(h1^2)*n^2, beta1=h1.beta, pow1=h1.pow, beta2=h2.beta, pow2=h2.pow,
          result[[1]])
    }, error.ix))
    layout(1)
})()

plot.target <- 'plot'

(function() {
    xlim <- c(0.5, 1.5)
    ylim <- c(0.1, 0.7)
    pdf(file.path(plot.target, 'asymptotic-growth-rate_p=0.pdf'),
        height=5, width=5,
        useDingbats=FALSE)
    on.exit(dev.off())
    plot(NA, xlim=xlim, ylim=ylim,
         xlab=expression('Spatial index'^{1/3}),
         ylab='Asymptotic growth rate',
         xaxt='n')
    axis(1, at=c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6))
    R <- 0.01
    g <- 1
    s <- 0
    lines(xlim, rep((2*pi*R*g^2)^(1/3), 2),
          col='darkgray',
          lwd=2,
          lty=1)
    ##
    x.theor <- seq(xlim[1]^3, xlim[2]^3, length.out=101)
    y.theor <- (function() {
        lambda <- (2*pi*R*g^2*x.theor)^(1/3)
        z <- (sqrt(54-lambda^3*s^6) + sqrt(54))^(1/3)
        lambda/sqrt(6)*(z + lambda*s^2/z)
    })()
    lines(x=x.theor^(1/3), y=y.theor, lty=1, lwd=2,
          col=hsv(2/3,0.2,1))
    ##
    for (p in list(
                   list(type='source-destination_p=0', col=hsv(2/3,0.54,0.72), pch=16)
                   )) {
        if (p$type == 'null') {
            n <- 1
            sp.ix <- rep(1, n)
        } else {
            n <- nrow(result.grate[[p$type]])
            sp.ix <- result.grate[[p$type]][,1]
        }
        matpoints(sp.ix^(1/3),
                  result.grate[[p$type]][1:n,'growth.rate'],
                  pch=p$pch,
                  col=p$col)
    }
})()
