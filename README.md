# survival-sacrifice-model

Primal-dual interior point method for estimation in the
survival-sacrifice model.

We give an Rcpp script for computing the real MLE for the two
distribution functions in the so-called survival-sacrifice model.

The method is described in the manuscript:  Antonio E. Gomes, Piet Groeneboom
and Jon A. Wellner. Nonparametric estimation of the lifetime and
disease onset distributions for a survival-sacrifice model:
https://www.stat.washington.edu/jaw/RESEARCH/PAPERS/survsacr.pdf, 2001.

This is presently the most efficient method known for computing the real
maximum likelihood estimator. The obvious alternative, the EM algorithm,
as already described in Bruce W. Turnbull and Toby J. Mitchell,
Nonparametric estimation of the distribution of time to onset for
specific diseases in survival/sacrifice experiments. Biometrics,
40(1):41, 1984. URL https://doi.org/10.2307/2530742 is very unwieldy for
the present model because of the very large number of "candidate" points
of mass it has to check in a 2-dimensional search, apart from being slow
by itself.

In contrast, the primal-dual interior point method has
super-linear convergence and only deals with 2n parameters (n being
sample size).

The present method is demonstrated by the R file surv_sacr.R for the data 
set of 109 mice in Table 1 of the above paper of Turnbull and Mitchell,
which is the data file RFM109.txt in the directory.

The program expects a data file with three columns of the same type as RFM109.txt.
The first column contains the observations on the times (ties allowed).
The second and third column contain indicators, where the first indicator equals 1
if time of onset of the disease studied is less than or equal to C (a censoring time,
corresponding to death from an unrelated cause) and zero otherwise, and the second
indicator equals 1 if the time of death from the disease studied is less than or equal
to C, and zero otherwise. The program gives back the maximum likelihood estimates of the
two distribution function. The are output$MLE1 and output$MLE2. The demonstration program
also gives a drawing of these functions.

See also: vander Laan, M.J. and Jewell, N.P (2003). Current status and right-censored
data structures when observing a marker at the censoring time. Annals of
Statistics 31, 512–535, in particular pages 514-515.
