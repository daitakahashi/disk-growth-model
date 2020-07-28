
library(Rcpp)
library(EBImage)

sourceCpp('diskGrowthModel2.cpp')

find.populationsOnRegion <- function(n, circles, image) {
    populationOnRegion(n, circles, image);
}

drawOnMatrix <- function(n, circles, target.size=n) {
    drawOnMatrixC(target.size, circles, target.size/n)
}

upscale.matrix <- function(m, new.dim) {
    current.dim <- nrow(m)
    outer(seq(0.5, new.dim - 0.5, length.out=new.dim)*current.dim/new.dim,
          seq(0.5, new.dim - 0.5, length.out=new.dim)*current.dim/new.dim,
          ## function(x,y) interp.im(m, y, x)),
          function(x,y) {
              ## No interporation...
              m[ceiling(cbind(x,y))]
          })
}

distributeInitialPopulation <- function(h, density=1e-6, N=NA, dx=1) {
    if (!missing(density)) {
        N <- max(1, rpois(1, length(h)*dx^2*density))
    }
    v <- sample.int(length(h), N, replace=TRUE, prob=h)
    x <- (floor((v - 1) %% nrow(h)) + runif(length(v), 0, 1))*dx + 0.5
    y <- (floor((v - 1) /  nrow(h)) + runif(length(v), 0, 1))*dx + 0.5
    cbind(X=x, Y=y, Radius=0, Age=0)
}

source('pinkNoise.R')

genHabitatMap <- function(hresolution, beta, pow, scale=1) {
    h <- pinkNoise(hresolution, beta=beta, scale=scale)
    h <- ((h - min(h))/(max(h) - min(h)))^pow
    h <- h/sum(h)
    return(h)
}

simpson.index <- function(m, xlength=n) {
  return(sum((m/sum(m))^2)/(nrow(m)/xlength))
}

asymptoticGrowthRate <- function(v, npop=NA) {
    if (max(v) < 0.01) {
        NA
    } else {
        if (!missing(npop)) {
            v <- v[npop >= 100 & npop <= 1000]
        } else {
            v <- v[-(1:5)]
            v <- v[v < 1e-2 & v > 1e-3]
        }
        w <- log(v)
        x <- 1:length(w)
        if (length(x) < 3) {
            return (NA)
        }
        g <- glm(w ~ x)
        coef(g)['x']
    }
}

run <- function(propagule.attraction, tm,
                R=0.01, g=1, s=1/3, n=2048, dt=1, dx=1,
                initial.density=1e-6,
                initial.population.number=NA,
                disperser.production=propagule.attraction,
                initial.populations=NA,
                stop.when.complete=FALSE,
                show.progress=TRUE,
                output.populations=TRUE,
                target.proportion=1,
                minimum.detectable.size=NA) {
    propagule.attraction <- propagule.attraction/sum(propagule.attraction)
    if (missing(initial.populations) || !('matrix' %in% class(initial.populations))) {
        if (missing(initial.population.number)) {
            initial.population.number <-
                max(1, rpois(1, sum(n^2*dx^2*initial.density)))
        }
        circles <- distributeInitialPopulation(propagule.attraction,
                                               N=initial.population.number,
                                               dx=dx*n/nrow(propagule.attraction))
    } else {
        circles <- initial.populations
    }
    disperser.pickup.freq <- ((nrow(disperser.production)/n)^2
                              *(R*dt*dx^2*n^2)
                              *upscale.matrix(disperser.production,
                                              n))
    gg <- g*dt/dx
    ##
    if (is.na(minimum.detectable.size)) {
        remove.detected.colonies <- FALSE
        minimum.detectable.size <- 0
    } else {
        remove.detected.colonies <- TRUE
    }
    if (length(tm) > 1) {
        niter <- as.integer(length(tm))
    } else {
        niter <- as.integer(tm)
    }
    if (output.populations == 'last') {
        output.populations <- FALSE
        keep.last <- TRUE
    } else {
        output.populations <- identical(output.populations, TRUE)
        keep.last <- FALSE
    }
    result <- runSimulation(niter,
                            as.numeric(gg),
                            s,
                            disperser.pickup.freq,
                            propagule.attraction,
                            circles,
                            stop.when.complete,
                            show.progress,
                            output.populations,
                            target.proportion,
                            remove.detected.colonies,
                            minimum.detectable.size)
    if (keep.last) {
        return(list(population.last=result[[2]],
                    stats=result[[1]]$stats))
    } else {
        return(result[[1]])
    }
}

