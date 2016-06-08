README for the SRE-tools
========================

This package contains a set of routines that can compute various
performance measures for NIST style speaker recognition evaluation
system output.  The main routines are written in R, an advanced
interpreted statistical programming language.  R, like this package,
is free software, and both are released under the GNU General Public
License, version 2.  This package also comes with some helper files
and contributed software.  

This README is written in the understanding that you know what has
been going on in SRE-2008.

Features
--------

The routines in this package allow
 - flexible plotting of DET curves
 - calculation of EER, minCdet, Cdet, minCllr and Cllr
 - adding legends to the plot in the order of the EER
 - plotting APE-plots, with proper legend
 - merging of meta-information to the scores, and manipulation of
   these tables
 - conditioning trials on meta-information 
 - weighting trials based on a chosen conditioning, before calculation
   of DET parameters and performance metrics.  This is pretty new,
   although George Doddington's DET_tools.pm has similar capability.

Installation
------------

Make sure you have a running version of R on your system.  It is
available from http://www.r-project.org , but it may already have been
instyalled on you system.  On Debian/GNU Linux installing could be
accomplished by issuing

$ sudo apt-get install r-base

or something similar.  On other Linux systems, you may have to compile
from sources or click your way though a whole lot of menus to find the
R package.  On other Unix system, your mileage may vary.  On non-Unix
systems, you're currently left on your own.  Feedback is appreciated.

Unpack the tarball in a directory where you have access to the
scores.  If you do _not_ have the package `fdrtool' installed (see
later for detection of this), you might have to install it using

$ export R_LIBS=$HOME/lib/R
$ mkdir -p $R_LIBS
$ R CMD INSTALL contrib/fdrtool_1.2.4.tar.gz

or similar.  On the Mac you can install the package from the menu
"Packages & Data", choose the "CRAN sources", and search for
"fdrtool".

The R tools rely on an external program `bin/add-target', so make
sure that "./bin" is in your path, e.g., 

$ PATH=./bin:$PATH

If this has been set up correctly, it should be possible to run the
small R session: (">" indicates the R prompt, "#" is a comment)

$ R
> source("det.R")			# load the tools
> x <- read("prim/tno-1") 		# read example SRE output
> plot(det(x), main="a lot of trials")	# make a DET plot over all results
       Cllr  Cllr.min      eer       Cdet   Cdet.min       mt        mn    nt
1 0.2501407 0.2237935 5.623803 0.03375549 0.03221064 6.341288 -5.142138 20449
     nn     n
1 78327 98776
> peps("sillydet.eps")			# print this graph to .eps
> q()					# quit, or ctrl-D
Save workspace image? [y/n/c]: y	# answer y
$

If you get an error at the "source" statement that a library cannot be
loaded, you don't have 'fdrtool' installed.  If you can't get the
installation of this package going, you can consider replacing the
line in det.R

  p.opt <- monoreg(p.ideal,w=w.ideal)$yf

to 

  p.opt <- isoreg(p.ideal)$yf

so then hopefully you will get most of the code running.  In this case,
the value of Cllr.min in the condition-weighted analysis will be
incorrect.  

Using the package, introduction to R
------------------------------------

The package contains several R functions.  

1) The first function is 

read(file)

which reads a NIST submission file `file', adds target information
(though the external program `add-target'), and returns an R `data
frame'---which is basically the table with the NIST submission file
information in it.  

Examples: (the ">" at the beginning is the R-prompt)

> source("det.R")
> x <- read("prim/tno-1")

reads the TNO submission into data frame x.  ("<-" is the assignment
operator.)

> sum(x$target) 
[1] 20449

There are 20449 target trials in the test.  (The [1] indicates that
this is the first value of an array of numbers.)

> table(x$target, x$gender)
       
            f     m
  FALSE 47184 31143
  TRUE  12159  8290

This operation gives a contingency table or target and non-target
trials, for male and female trials. 

> mean(x$score)
[1] -2.764793

The mean score over all trials.  Not very useful, but

> tapply(x$score, x$target, mean)
    FALSE      TRUE 
-5.142138  6.341288

gives an impression where target and non-target scores hang out.  We
can also condition this on speaker's sex

> tapply(x$score, list(x$gender,x$target), mean)
      FALSE     TRUE
f -4.500746 5.949503
m -6.113895 6.915922

But usually, we don't look too deeply into the plain data.  What I
often seem to do, is to wonder what information there actually is in
the data frame.  I found that a indexing a data frame with [1,] like in 

> x[1,]
   mcond adapt  tcond gender model  test chan   dec     score target
1 short2     n short3      f 10017 fbnep    a FALSE -5.569686  FALSE

gives me the first row in the table---but I am usually more interested
in the column names---the identifiers I can use after the "x$". 


2) The second function in det.R is 

det(x, cond)

which takes a NIST submission data frame `x', and computes basic
detection statistics, and returns them in a data structure of class
"det".  Optionally, conditioning factors for equalizing the weights of
several sub-conditions of the data can be given in `cond'.  (This
subject will be treated later). 

Examples:

> xd <- det(x)
> summary(xd)
       Cllr  Cllr.min      eer       Cdet   Cdet.min       mt        mn    nt
1 0.2501407 0.2237935 5.623803 0.03375549 0.03221064 6.341288 -5.142138 20449
     nn     n
1 78327 98776

3) summary.det() is a function in det.R and shows some basic
performance measures that are recorded in the class "det".  "mt" and
"mn" mean mean of target and non-target scores, and "nt" and "nn" are
the number or target and non-target trials.  Note, that this is not a
very sensible example, since we've been plainly pooling over all
trials, but here we can show some handy R tricks again:

> x <- merge.meta(x,meta.2008)

or shorter,

> x <- merge.meta(x)

`meta.2008' is a data frame containing meta-information given before
the evaluation, which was read from "det.R" at the beginning of the
session.  The function merge.meta() from det.R adds additional
meta-information for each trial to the data frame `x'.

> x.1 <- subset(x, mtype=="interview" & ttype=="interview")
> xd.1 <- det(x.1)
> summary(xd.1)
       Cllr  Cllr.min      eer       Cdet   Cdet.min       mt        mn    nt
1 0.2375118 0.2124601 5.625438 0.03009871 0.02370526 4.468742 -4.134552 11540
     nn     n
1 22641 34181

The `subset' function from R can be used to select trials from `x'.
This example showed the NIST common condition 1.  We tend to use the
letters `m' and `t' for model (or train) and test segment.

4) Well, now it is time to plot this more sensible sub-set of the trials:

> plot(xd.1)
       Cllr  Cllr.min      eer       Cdet   Cdet.min       mt        mn    nt
1 0.2375118 0.2124601 5.625438 0.03009871 0.02370526 4.468742 -4.134552 11540
     nn     n
1 22641 34181

(These results confirm the official results in the SRE-2008 results
release file "results/1/tno_1.1.det.png"---which is always
encouraging).  We get a summary for free here---and a nice DET plot,
with 95% confidence boxes around the actual decision operating point,
a circle at Cdet.min, and plenty of room for a legend...  

The function plot.det (which is called implicitly for plot() on an
object of class "det") has many options, but the most useful ones are

plot(det, nr=1, lty=1, ...)

`det' is the data structure of class "det", `nr' is used for plotting
colour, and `lty' for the line type.  If nr==1 and lty==1, a new
plotting canvas is created and the DET grid is made, otherwise, the
DET curve is added to an existing plot.  For example, do

> plot(xd,2)

to add the previous DET in the same graph. 

With the meta data you can filter existing score data frames, but it
is also possible to use an external filter program while reading the
data.  The trick is that the function read() actually `(z)cat's the file
and reads from a pipeline, so you can say:

> x.6 <- read("prim/tno-1 | nc-2008-filter --cond 6")
> plot(det(x.6), main="NIST condition 6")
       Cllr  Cllr.min      eer       Cdet   Cdet.min      mt       mn   nt
1 0.2518594 0.2181511 5.637018 0.03285049 0.02861965 9.74785 -5.26091 2678
     nn     n
1 33218 35896

Here we used the example filter program `nc-2008-filter' that filters
trials of a particular NIST common condition according to the trial
key.

5) Trials weighting.  In SRE-2008 the different acoustic conditions in
the requires core test ``short2-short3' have all different number of
trials.  The proportions vary a lot.  This can be appreciated from the
det.R helper function cond.table()

cond.table(x, cond, target)

This function produces a frequency table of trials in conditions
specified by the list of factors in `cond', possibly cross-classified
with target/non-target condition as well, if `target' is TRUE.  Example:

> cond.table(x, list(mtype, mmic, ttype, tmic), target=T) 
       mtype mmic     ttype tmic target  Freq
1  interview  mic interview  mic  FALSE 22641
4  phonecall  phn interview  mic  FALSE  4850
8  phonecall  phn phonecall  mic  FALSE  6982
13 interview  mic phonecall  phn  FALSE 10636
16 phonecall  phn phonecall  phn  FALSE 33218
17 interview  mic interview  mic   TRUE 11540
20 phonecall  phn interview  mic   TRUE  2500
24 phonecall  phn phonecall  mic   TRUE  1472
29 interview  mic phonecall  phn   TRUE  1105
32 phonecall  phn phonecall  phn   TRUE  3832

Only conditions with actual trials in that condition are returned.  Out
of the 32 possible combinations, only 10 occurred in short2-short3.
The problem now is that naive pooling of trials like in
"plot(det(x))" shows performance that is very unequally weighted, e.g., more
than half the target trials come from the int-int condition.

The function det() can equalize these unequal trial counts.  By saying

> xd <- det(x, list(mtype, mmic, ttype, tmic)) ## may take some time...

the contributions of each of the 5 different acoustical conditions for
which there are trials available are weighted equally.  This is
accomplished by weighting each trial according to the inverse of the
relative frequency of the trial's condition.  Now, we can obtain more
sensible over-all-trials performance metrics:

> summary(xd)
       Cllr  Cllr.min      eer       Cdet   Cdet.min       mt        mn    nt
1 0.2328931 0.1988914 4.998876 0.02826864 0.02599788 6.341288 -5.142138 20449
     nn     n
1 78327 98776

6) Finally in order to stimulate researchers to carry out
over-all-trials evaluation, we have a utility function

plot.cond(x, cond, ...)

that will split the trial data in `x' according to the conditions
obtained from the list of factors in `cond', and will make separate
conditioned DET plots and analysis.  It will also show an
over-all-conditions, equally weighted, analysis.  For example:

> plot.cond(x.1, gender)
         Cllr  Cllr.min      eer       Cdet   Cdet.min       mt        mn    nt
all 0.2303026 0.2045689 5.413335 0.02960904 0.02311025 4.468742 -4.134552 11540
f   0.2768975 0.2522662 6.719676 0.03283787 0.02663579 3.952384 -3.545275  6639
m   0.1837077 0.1486914 4.040198 0.02638020 0.01855576 5.168213 -4.949086  4901
       nn     n
all 22641 34181
f   13137 19776
m    9504 14405

will take the trial data from `x.1' (interview-interview trials), and
split according to gender.  But more interestingly, the "all" pooled
condition is actually an equal-weight mix of the "f" and "m"
conditions.  Because there were more female trials and female
performance is worse than male, the original trial-based pooling led
to a slightly pessimistic view on the performance, compared to what we
believe is better: equalizing the condition weights. 

As a final example, we can make the contributions of the 5 acoustical
conditions in SRE-2008 visible with the single command

> plot.cond(x, list(mtype, ttype, tmic), main="Condition weights equalized")
                             Cllr  Cllr.min      eer       Cdet   Cdet.min
all                     0.2328931 0.1988914 4.998876 0.02826864 0.02599788
interview interview mic 0.2375118 0.2124601 5.625438 0.03009871 0.02370526
phonecall interview mic 0.2411139 0.1575142 4.395876 0.02972124 0.02110144
phonecall phonecall mic 0.2378819 0.1448310 4.009232 0.02358522 0.01577032
interview phonecall phn 0.2260216 0.1899311 5.349262 0.02787835 0.02392979
phonecall phonecall phn 0.2219362 0.1906779 4.900495 0.03005969 0.02479888
                               mt        mn    nt    nn     n
all                      6.341288 -5.142138 20449 78327 98776
interview interview mic  4.468742 -4.134552 11540 22641 34181
phonecall interview mic  5.597978 -6.482080  2500  4850  7350
phonecall phonecall mic  8.958881 -6.747623  1472  6982  8454
interview phonecall phn  5.606874 -5.251124  1105 10636 11741
phonecall phonecall phn 11.671634 -5.260909  3832 33218 37050
> dev.print(pdf, "tno-1-all.pdf")

Note, that all performance metrics, Cdet, Cdet.min, EER, Cllr and even
Cllr.min are calculated considering the weighting of conditions.  

Further reading
---------------

More arguments and a more mathematical description of the conditioned
weighting is included as doc/cond-weight.pdf . 

For novices, it can be useful to read the NIST SRE-2008 evaluation plan, 
http://www.nist.gov/speech/tests/sre/2008/sre08_evalplan_release4.pdf

David van Leeuwen, david.vanleeuwen@gmail.com
Berkeley, CA, 9 June 2008. 
