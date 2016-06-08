`read.sre` <- function(file) {
  data(srekey)                          # is this what we want?
  pipere <- "\\|\\s*$"                  # detect this is a pipe cmdline
  if (length(grep(pipere, file))) {
    cmd <- sub(pipere, "", file)
    x <- read.table(pipe(cmd))
  }
  else
    x <- read.table(file)               # for non-unix environments
  names(x) <- c("mcond", "adapt", "tcond", "gender",
                "model", "test", "chan", "dec", "score")
  if (class(x$dec)=="factor")
    x$dec <- x$dec=='t'                 # in case this was "t|f"
  ## now we must guess the trial set
  x$model <- as.character(x$model)
  mm <- paste("m", x$model, sep="")
  sre <- unique(srekey$sre[mm,])         # what SREs are the models from?
  sre <- sre[!is.na(sre)]
  if (length(sre) != 1) 
    warning(paste("Not a unique evaluation, it appears", sre))
  ## merge in the key
  x <- merge(x, srekey[[sre]])
  class(x) <- list("sre", class(x))
  invisible(x)
}

