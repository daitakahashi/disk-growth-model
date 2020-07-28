
source('diskGrowthModel2.R')

library(parallel)

if (exists('cl')) {
    stopCluster(cl)
    rm(cl)
}
cl <- makeCluster(detectCores()/2, outfile='')

invisible(clusterEvalQ(cl, source('diskGrowthModel2.R')))


data.dir <- 'output'

n <- 2048
R <- 0.01
g <- 1
s <- 0

n.iterations <- 200

clusterExport(cl, 'data.dir')

## compute source-mediated, destination-mediated, and null models
set.seed(12345)
for (dispersal.type in c('source', 'null')) {
    random.seeds <- sample.int(.Machine$integer.max, n.iterations)
    iter.args <- Map(list, 1:n.iterations, random.seeds)
    ##
    print(paste0('start computing "', dispersal.type, '"'))
    system.time(
        res <- parLapply(cl, iter.args, function(arg, n, R, g, model.type) {
            ix <- arg[[1]]
            seed <- as.integer(arg[[2]])
            set.seed(seed)
            ##
            h.beta <- runif(1, -5, -2)
            h.pow <- runif(1, 1, 5)
            h <- genHabitatMap(n, h.beta, h.pow)
            if (model.type == 'source') {
                source.distribution <- h
                destination.distribution <- matrix(rep(1, 4), 2, 2)/4
            } else if (model.type == 'destination') {
                source.distribution <- matrix(rep(1, 4), 2, 2)/4
                destination.distribution <- h
            } else {
                source.distribution <- matrix(rep(1, 4), 2, 2)/4
                destination.distribution <- matrix(rep(1, 4), 2, 2)/4
            }
            ##
            res <- lapply(1:100, function(ix2) {
                run(propagule.attraction=destination.distribution,
                    disperser.production=source.distribution,
                    tm=1:500,
                    R=R,
                    g=g,
                    s=0,
                    n=n,
                    dt=1,
                    dx=1,
                    stop.when.complete=TRUE,
                    show.progress=FALSE,
                    output.populations=FALSE)$stats
            })
            result <- list(c(ix=ix, seed=seed,
                             h.beta=h.beta, h.pow=h.pow,
                             spatial.ix=sum(h^2)/sum(h)^2*length(h),
                             R=R, g=g, n=n),
                           res)
            dir.create(file.path(data.dir, model.type), showWarnings=FALSE)
            save(result, file=file.path(data.dir, model.type, paste0(ix, '.RData')))
            print(paste('index', ix, 'done'))
            return (result)
        }, n=n, R=R, g=g, model.type=dispersal.type)
    )
}

for (target.dir in Sys.glob(file.path(data.dir, '*'))) {
    fname <- file.path(data.dir, paste(basename(target.dir), 'RData', sep='.'))
    result <- lapply(Sys.glob(file.path(target.dir, '*.RData')),
                     function(x) {
                         load(x)
                         return (result)
                     })
    save(result, file=fname)
}
