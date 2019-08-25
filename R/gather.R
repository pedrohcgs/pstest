# Functions used to guarantee replicability of the bootstrap
# These are based on the harvestr package, but since harvestr is not on CRAN anymore, we just use these functions
# safe version of retrieving the .Random.seed -
get.seed.ps <- function(){
  if(exists(".Random.seed", envir=.GlobalEnv, mode="numeric")) {
    seed <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    class(seed) <- c("rng-seed", "integer")
    seed
  } else {
    NULL
  }
}

gather.ps <- function(x, seed=1234){
  set.seed(seed, kind="L'Ecuyer-CMRG", normal.kind="Inversion")
  r <- get.seed.ps()
  seeds <- vector('list', x)
  for(i in seq_len(x)) {
    r <-
      seeds[[i]] <-
      structure( parallel::nextRNGStream(r)
                 , RNGlevel='stream'
                 , class=c("rng-seed", "integer")
      )
  }
  structure(seeds, class=c('rng-seeds', 'list'))
}
