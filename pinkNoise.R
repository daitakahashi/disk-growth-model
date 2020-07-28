
## reference:
##   "Generate spatial data (1/f noise)"
##       https://nl.mathworks.com/matlabcentral/fileexchange/5091-generate-spatial-data
##   Lennon (2000). Red-shifts and red herrings in geographical ecology. Ecography.
##   Keitt (2000). Spectral representation of neutral landscapes. Landscape Ecology.
## Also, we refer Chipperfield et al. (2011). for necessity of phase correction.
##
## A pink noise is a spatially-autocorrelated noise with a power-law--power-
## spectral-density distribution.
##
## The algorithm below first generates a radial-symmetric--power-law-shaped
## power-spectral--density-distribution Sf, then randomize phases by phy,
## and finally generate the spatial noise by IDFT.
##
## Note that the phase shifts are corrected, i.e., arranged to be "antisymmetric",
## to keep the IDFT results being real values (Re() at the end of the function is
## simply for type coercion).
## The correction is necessary to generate a spatial noise with a desired power spectrum.
##
## The parameter beta is a scaling exponent of the power-spectral density,
## which is the same as beta in Lennon (2000) and satisfies beta = -(1 + 2*H) (Keitt, 2000).
## Note that large values of the scaling exponent beta diverges intergration of Sf, leads to
## numerically undesirable results.
## 
## Argument scale (0 < scale <= 1) is an experimental high-pass cutoff frequency. A smaller
## value of the scale result in more frequent ups and downs by removing large-scale trends,
## as a result, give more peaks.
## Note that the noise is no longer self similar when the scale is not equal to 1.
##
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

pinkNoise <- function(size, beta=-2, scale=1, freq.factors=rep(1, size^2)) {
    freqs <- c(0:(floor(size/2)), -((ceiling(size/2) - 1):1))/size
    ## generate a frequnency-distribution matrix Sf
    Sf <- outer(freqs, freqs, function(fx, fy) {
        return ((fx^2 + fy^2)^(beta/2))
    })
    Sf[!is.finite(Sf)] <- 0
    Sf[outer(freqs, freqs, function(fx, fy) {
        return (sqrt(fx^2 + fy^2) < 1/scale/size)
    })] <- 0
    ##
    Sf <- Sf*abs(freq.factors)
    return(noise2D.ift(Sf))
}
