## det.R.  DET plotting routines for NIST SRE-2008 and earlier
## (c) 2008 David A. van Leeuwen, TNO, ICSI

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, version 2 of the
## License.

## $Revision: 1.1 $

## converts propabilities as percentage to quantiles of the normal distribution
nd <- function (p) {
  qnorm(p/100)
}

## note, as per 1 jan 2007 we don't normalize the DCF anymore.
## Too ugly!
DCF <- function (fa, miss, prior=0.01, cfa=1, cmiss=10, norm=1) {
  norm*(cmiss*miss*prior + cfa*fa*(1-prior))
}

## This function read a plain NIST result file (9-column format since 2005),
## adds target (and language if available) information, and returns a
## data.frame that can be read by subsequent functions like det().
## This should now read multiple years (2005--2008 in pin-name .key format)
## You need an external program `add-target` that adds target information.
## Note, that since we actually `cat` and read from a pipe, you may add
## other commands on a pipeline, e.g.
## > x <- read("out/short-short-2006-0.gz | myfilter")
read <- function(file) {
  if (length(grep("\\.gz", file))) cat <- "zcat"
  else cat <- "cat"
  x <- read.table(pipe(paste(cat, file, "| tr 'A-Z' 'a-z' | add-target")))
  names(x) <- c("mcond", "adapt", "tcond", "gender",
                      "model", "test", "chan", "dec", "score",
                      "target")
  x$dec <- x$dec=='t'
  invisible(x)
}

## This function takes a NIST result data frame, and computes relevant
## statistics (Cllr, Cdet, eer) and prepares for plotting.
## You may specify a condition that will weight trials according to their
## inverse frequency, i.e., compensating for their proportion difference.
## Examples:
## Plot a det of all trials of x:
## > plot(det(x))
## Compute DET statistics for x, equalizing trial counts over gender:
## > xd <- det(x, list(gender))
## > xd$eer
## Just give me the equal error rate for female
## > det(subset(x,gender=="f"))$eer
## Give some basic performance statistics
## > data.frame(det(x)[1:9])
det <- function(x, cond) {
  if (!missing(cond)) {
    ct <- eval(substitute(cond.table(x, cond, target=TRUE)),x, parent.frame())
    nvar <- ncol(ct)-1
    mf <- tapply(ct$Freq, ct$target, mean)
    w <- as.vector(mf[ct$target]/ct$Freq)
    ct <- transform(ct, weight=w)
    x <- merge(x, ct, by=names(ct)[1:nvar]) # add weight = 1/Freq
    x <- x[order(x$score),]             #sort by score
  } else {
    x <- transform(x, weight=1)
    x <- x[order(x$score),]             #sort by score
    ct <- NULL                          #cond.table
  }
  t <- as.numeric(x$target)
  w <- x$weight
  nt <- sum(t*w)                        # number of target trials
  nn <- sum((1-t)*w)                    # number of non-target trials
  miss <- rangecheck(cumsum(t*w)/nt) # the cummlative (weighted) probs
  fa <- rangecheck(1 - cumsum((1-t)*w)/nn)
  nmiss <- sum((x$target & !x$dec)*w) # weighted counts...
  nfa <- sum((!x$target & x$dec)*w)
  ## adcf confidence interval
  if (!is.null(ct)) {                   # estimated through normal
    afa <- nfa/nn                       # approximation and weighted counts
    ci <- qnorm(0.975)*sqrt(afa*(1-afa)/nn)
    afa.lci <- afa-ci
    afa.uci <- afa+ci
    amiss <- nmiss/nt
    ci <- qnorm(0.975)*sqrt(amiss*(1-amiss)/nt)
    amiss.lci <- amiss-ci
    amiss.uci <- amiss+ci
  } else {                              # "true" binomial confidence intervals
    fa.ci <- data.frame(ci.binomial(nfa, nn))
    afa <- fa.ci$p.x.n                  # actual false alarms
    afa.lci <- fa.ci$p.lci
    afa.uci <- fa.ci$p.uci
    miss.ci <- data.frame(ci.binomial(nmiss, nt))
    amiss <- miss.ci$p.x.n              # actual misses
    amiss.lci <- miss.ci$p.lci
    amiss.uci <- miss.ci$p.uci
  }
  ## actual, minimum DCF
  adcf <- DCF(afa, amiss)
  dcf <- DCF(fa, miss)
  mi <- which.min(dcf)
  min.score <- x$score[mi]
  ## Cllr
  cllr <- Cllr(x)
  x <- opt.llr(x, laplace=F)
  cllr.min <- Cllr(x, opt=T)
  ## EER
  eeri <- which.min(abs(fa-miss))
  ## means of target and non-target scores
  mt <- mean(x$score[x$target])
  mn <- mean(x$score[!x$target])
  res <- list(Cllr=cllr, Cllr.min=cllr.min, eer=50*(fa[eeri]+miss[eeri]),
              Cdet=adcf, Cdet.min=dcf[mi],
              mt=mt, mn=mn, 
              nt=nt, nn=nn, n=nt+nn,
              afa=afa, amiss=amiss, afa.lci=afa.lci, afa.uci=afa.uci,
              amiss.lci=amiss.lci, amiss.uci=amiss.uci,
              mfa=fa[mi], mmiss=miss[mi], atscore=min.score,
              fa=fa, miss=miss, data=x, cond.table=ct)
  class(res) <- "det"
  invisible(res)
}

## This function tabulates exisiting conditons in a trial set.
## Empty conditions are removed.  Optionally, target information is
## implicitly conditioned in. 
## `cond' is a list of factors, e.g.
## > cond.table(xm, list(gender, mc))
cond.table <- function(x, cond, target=F) {
  if (!missing(cond)) {
    l <- eval(substitute(cond),x, parent.frame())
    nl <- as.list(1:ncol(x))
    names(nl) <- names(x)
    if (target) {                       # I am sure this can be done better
      t <- as.data.frame(table(data.frame(l,x$target)))
      names(t)[ncol(t)-1] <- "target"
      nf <- ncol(t)-2
    } else {
      t <- as.data.frame(table(l))
      nf <- ncol(t)-1
    }
    names(t)[1:nf] <- c(names(nl)[unlist(eval(substitute(cond), nl))])
  }
  t <- t[t$Freq>0,]
  t
}


summary.det <- function(x) data.frame(x[1:10])

## This function plot a DET curve.
## When called with x=NULL, no data is plotted (useful for plotting a frame)
## When called with nr=1 and lty=1 (defaults), a new frame is created
## Argument `x' is an object of class "det"
## Argument 'nr' (or 'col') set the color according to the current palette()
## Argument `lty' sets the line type. 
plot.det <- function(x, nr=1, lty=1, col=nr, optimize=T,
                     xmin=0.1, xmax=50, ymin=0.1, ymax=50,
                     xlab="false alarm probability (%)",
                     ylab="miss probability (%)",
                     ...) {
  xlim <- c(nd(xmin), nd(xmax))         # 0.5 % seems accurately enough
  ylim <- c(nd(ymin), nd(ymax))
  par(pty="s", cex.axis=1)
  if (is.null(x))                      # only produce frame...
    xdata <- ydata <- numeric(0)        # empty data set
  else {
    attach(x)
    size <- length(fa)
    if (optimize) {
      changes <- diff(diff(fa)!=0)<0 | diff(diff(miss)!=0)<0
      sample <- c(T, changes, T)
    } else sample <- 1:size
    xdata <- qnorm(fa[sample])
    ydata <- qnorm(miss[sample])
    detach(x)
  }
  if (nr==1 && lty==1) {                # first time, plot everything..
    plot(xdata, ydata, type="l", xaxt='n',yaxt='n', xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim, lwd=2, col=col, lty=lty, ...)
    ## draw the grid and axes
    l <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 60, 80, 90, 95, 98, 99, 99.5, 99.8, 99.9)
    qnl=qnorm(l/100)
    for (d in 1:2) axis(d,qnl,l)
    abline(h=qnl, lty=3)
    abline(v=qnl, lty=3)  
    ## y=x for equal error rate
    abline(coef=c(0,1), lty=3)
  }
  else
    lines(xdata, ydata, type="l", lty=lty, lwd=2, col=col)
  ## actual DCF
  if (is.null(x)) return(NULL)
  rect(qnorm(x$afa.lci), qnorm(x$amiss.lci),
       qnorm(x$afa.uci), qnorm(x$amiss.uci), border=col, lwd=2)
  points(qnorm(x$mfa), qnorm(x$mmiss), pch=1, cex=2, col=col, lwd=2)
  summary(x)
}

## this function compensates for numerical rounding problems while using
## weighted trial counts. 
rangecheck <- function(x) {
  x[x<0]=0
  x[x>1]=1
  x
}

## this function might be handy for putting a legend in the plot.
## Example:
## plot(det(x))
## plot(det(y),2)
## legend.det(c("our system", "their system"))
legend.det <- function(legend, where="ur", order=1:length(legend), col,
                       lty, ...) {
  if (where == "ur") {
    x <- nd(50)
    y <- nd(50)
    xj <- 1
    yj <- 1
  } else {
    x <- nd(0.1)
    y <- nd(0.1)
    xj <- 0
    yj <- 0
  }
  n <- length(legend)
  if (missing(col)) 
    col=1:n
  if (missing(lty))
    lty=rep(1,n)
  legend(x, y, xjust=xj, yjust=yj, legend=legend[order], col=col[order], lwd=2, bg="white", lty=lty[order])
}

## plot several det curves according to specified conditions, starting with
## an overall, weighted per condition, plot.
## Example:
## x <- merge.meta(x,meta.2008)
## plot.cond(x, list(mtype,ttype,tmic))
plot.cond <- function(x, cond, nr=1, equalize=T, ...) {
  ## make an empty canvas
  r <- plot.det(NULL, ...)
  if (!missing(cond)) {
    ct <- eval(substitute(cond.table(x, cond, target=FALSE)),x, parent.frame())
    nf <- ncol(ct)-1                   # number of factors
    for (i in 1:nrow(ct)) {
      s <- merge(x,ct[i,], by=names(ct)[1:nf]) # poor man's subset
      nr <- nr+1
      p <- plot(det(s),nr)
      if (nf>1)
        row.names(p) <- apply(ct[i,1:nf], 1, paste, collapse=" ")
      else
        row.names(p) <- ct[i,1]
      r <- rbind(r, p)
    }
    ## plot overall DET last, in order to paint it "on top"
    nr <- nr+1
    if (equalize)
      p <- eval(substitute(plot(det(x, cond), nr, col=1, ...)),x, parent.frame())
    else
      p <- plot(det(x), nr, col=1, ...)
  } else p <- plot(det(x), 2, col=1, ...)
  row.names(p) <- "all"
  r <- rbind(p,r)                       # but put all on top...
  legend.det(legend=row.names(r), order=order(r$eer, decreasing=T))
  r
}

## From here some routines for calculating Cllr, minCllr and plotting APE. 

## Niko's LLR cost function
Cllr <- function(x, opt=F) {
  if (opt) {                            # do we want optimum?
    if (is.null(x$opt.llr)) x <- opt.llr(x, F) # did we have optimum?
    x$score <- x$opt.llr
  }
  if (is.null(x$weight)) x <- transform(x, weight=1) # did we have weights?
  attach(x)
  c.miss <- 1/log(2) * mean(log(1+exp(-score[target]))*weight[target])
  c.fa <- 1/log(2) * mean(log(1+exp(score[!target]))*weight[!target])
  detach(x)
  (c.miss+c.fa)/2
}

## after Niko's opt_logr.m
opt.llr <- function(x, laplace=T) {
  o <- order(x$score)
  p.ideal <- as.numeric(x$target[o]) ## ideal posterior
  if (is.null(x$weight)) x <- transform(x, weight=1) # did we have weights?  
  w.ideal <- x$weight[o]
  nt <- sum(p.ideal)
  nn <- sum(1-p.ideal)
  if (laplace) {
    p.ideal <- c(1,0,p.ideal,1,0) ## lapace's rule of succession
    w.ideal <- c(1,1,w.ideal,1,1)
  }
  p.opt <- monoreg(p.ideal,w=w.ideal)$yf
  if (laplace) 
    p.opt <- p.opt[3:(length(p.opt)-2)]
  post.log.odds <- log(p.opt)-log(1-p.opt)
  log.prior.odds <- log(nt/nn)
  llrs <- post.log.odds - log.prior.odds
  llrs[o] <- llrs
  transform(x, opt.llr=llrs)
}

read.meta <- function() {
  m <- read.table("ndx/short2-short3-2008.ndx")
  names(m) <- c("model", "s", "sphchan", "mlang",
                "mtype", "mmic", "tlang", "ttype", "tmic")
  rownames(m) <- paste(m$model, m$sphchan, sep="-")
  m
}

merge.meta <- function(x,m=meta.2008) {
  row.names(x) <- paste(x$model, x$test, x$chan, sep="-")
  merge(x,m,by=0)
}

## routines for PAV plots...  Not yet in weighted version...

# p = sigmoid(log odds)
sigmoid <- function(lo) 1/(1+exp(-lo))

# again shamelessly copied
bayes.error.rate <- function(tar, non, lpo) {
  n <- length(lpo)
  pmiss <- pfa <- vector(mode="numeric", length=n)
  for (i in 1:n) {
    p <- sigmoid(tar+lpo[i])
    pmiss[i] <- mean(1-sign(p-0.5))/2
    p <- sigmoid(non+lpo[i])
    pfa[i] <- mean(1-sign(0.5-p))/2
  }
  pmiss*sigmoid(lpo) + pfa*sigmoid(-lpo)
}

## ape.plot() makes an Applied Probability of Error plot for a number
## of "det" data structures.  Note: it has not been checked that
## weighted trials have been implemented correctly!
ape.plot <- function(data, legend="Ape plot", bar=T, bw=F) {
  def.par <- par(no.readonly=T)          # save parameters...
  if (class(data)=="det") data <- list(data)
  nr <- 1+as.numeric(bar)
  nsys <- length(data)
  layout(matrix(c(1,1:nsys,rep(nsys+1,nsys+1)),nr,nsys+1, byrow=T),
         c(1,rep(5,nsys)), c(3,2))
  legend <- rep(legend, nsys)
  clog <- clog.min <- vector(mode="numeric", length=nsys)
  lpo <- seq(-7, 7, by=0.1)
  pe <- pe.min <- pe.ref <- matrix(nrow=length(lpo), ncol=nsys)
  for (i in 1:nsys) {
    x <- data[[i]]$data
    x <- opt.llr(x, laplace=F) ## add optimal llr columns
    clog[i] <- Cllr(x)
    clog.min[i] <- Cllr(x, opt=T)
    pe[,i] <- bayes.error.rate(x$score[x$target], x$score[!x$target], lpo)
    pe.min[,i] <- bayes.error.rate(x$opt.llr[x$target],
                                   x$opt.llr[!x$target], lpo)
    pe.ref[,i] <- bayes.error.rate(0, 0, lpo);
  }
  max.e <- max(pe)
  par(cex=1)
  if (bw) {
    col <- rep(1,4);
    lty <- c(1,2,3,2)
  } else {
    col <- c("red", "green", "black", "purple")
    lty <- c(1,1,2,2)
  }
  for (i in 1:nsys) {
    if (i==1)
      par(mar=c(5,4,2,1))
    else
      par(mar=c(5,2,2,1))
    plot(lpo, pe[,i], col=col[1], type="l", lwd=2, main=legend[[i]],
         ylim=c(0,max.e), panel.first=grid(2, NULL, lty=lty[1]),
         xlab="log prior odds", ylab="probability of error")
    lines(lpo, pe.min[,i], col=col[2], lty=lty[2], lwd=2)
    lines(lpo, pe.ref[,i], col=col[3], lty=lty[3], lwd=2)
    abline(v=log(1/9.9), col=col[4], lty=lty[4])
  }
  if (bar) {
    bars <- rbind(clog.min, clog-clog.min)
    par(mar=c(3,4,0,0))
    if (bw) {
      col <- grey((2:1)/3)
    } else {
      col <- c("green", "red")
    }
    b <- barplot(bars, col=col, beside=F, ylab="Cllr (bits)")
    m <- max(clog)
    mx <- mean(b)
    legend(mx, 0.1*m, xjust=0.5, yjust=0,
           legend=c("calibration loss", "discrimination loss"),
           fill=rev(col), bg="white")
  }
  par(def.par)                          # reset parameters
}

## routine to do copy a graph to .eps
peps <- function(file) 
  dev.print(postscript, file=file, horizontal=F, paper="special", width=7, height=7)


## Finally, read in some constants and libraries
meta.2008 <- read.meta()
source("contrib/binomial-confidence.R")
library(fdrtool)                        # for monoreg

