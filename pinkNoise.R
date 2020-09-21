
## reference:
## Jon Yearsley (2020). Generate spatial data (https://www.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data), MATLAB Central File Exchange. Retrieved September 21, 2020.
##   Lennon (2000). Red-shifts and red herrings in geographical ecology. Ecography.
##   Keitt (2000). Spectral representation of neutral landscapes. Landscape Ecology.
##
## The parameter beta is a scaling exponent of the power-spectral density,
## which is the same as beta in Lennon (2000) and satisfies beta = -(1 + 2*H) (Keitt, 2000).

noise2D.ift <- function(Sf) {
    size <- ncol(Sf)
    ## generate an antisymmetric phase-shift matrix phi
    if (size %% 2 == 0) {
        n <- size/2 - 1
        phi.h0 <- runif(n)
        phi.v0 <- runif(n)
        phi.h1 <- runif(n)
        phi.v1 <- runif(n)
        phi.in11 <- matrix(runif(n^2), nrow=n)
        phi.in21 <- matrix(runif(n^2), nrow=n)
        phi.in12 <- -t(apply(apply(phi.in21, 2, rev), 1, rev))
        phi.in22 <- -t(apply(apply(phi.in11, 2, rev), 1, rev))
        phi <- rbind(c(0, phi.h0, 0, -rev(phi.h0)),
                     cbind(phi.v0, phi.in11, phi.v1, phi.in12),
                     c(0, phi.h1, 0, -rev(phi.h1)),
                     cbind(-rev(phi.v0), phi.in21, -rev(phi.v1), phi.in22))
    } else {
        n <- floor(size/2)
        phi.h0 <- runif(n)
        phi.v0 <- runif(n)
        phi1 <- matrix(runif(2*n^2), nrow=n)
        phi <- rbind(c(0, phi.h0, -rev(phi.h0)),
                     cbind(phi.v0, phi1),
                     cbind(-rev(phi.v0), -t(apply(apply(phi1, 2, rev), 1, rev))))
    }
    return(Re(fft(sqrt(Sf)*complex(argument=2*pi*phi)/size^2,
                  inverse=TRUE)))
}

pinkNoise <- function(size, beta=-2) {
    freq.factors <- rep(1, size^2)
    freqs <- c(0:(floor(size/2)), -((ceiling(size/2) - 1):1))/size
    ## generate a frequnency-distribution matrix Sf
    Sf <- outer(freqs, freqs, function(fx, fy) {
        return ((fx^2 + fy^2)^(beta/2))
    })
    Sf[!is.finite(Sf)] <- 0
    Sf[outer(freqs, freqs, function(fx, fy) {
        return (sqrt(fx^2 + fy^2) < 1/size)
    })] <- 0
    ##
    Sf <- Sf*abs(freq.factors)
    return(noise2D.ift(Sf))
}
