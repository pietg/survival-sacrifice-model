# survival-sacrifice-model

We give an R script for computing the real MLE for the two
distribution functions in the so-called survival-sacrifice model by the primal-dual interior point method and the EM algorithm. It can be run for the data file MLE_RFM109.txt from the R script MLE_RFM109.R and uses the R package Rcpp.

The primal-dual interior point method is described in the manuscript:  Antonio E. Gomes, Piet Groeneboom
and Jon A. Wellner. Nonparametric estimation of the lifetime and
disease onset distributions for a survival-sacrifice model:
http://dutiosb.twi.tudelft.nl/%7Epietg/survsac-v2.pdf (2018)

The primal-dual interior point method is presently the most efficient method known for computing the real
maximum likelihood estimator. The obvious alternative, the EM algorithm,
as already described in Bruce W. Turnbull and Toby J. Mitchell,
Nonparametric estimation of the distribution of time to onset for
specific diseases in survival/sacrifice experiments. Biometrics,
40(1):41, 1984. URL https://doi.org/10.2307/2530742 is very unwieldy for
the present model because of the very large number of "candidate" points
of mass it has to check in a 2-dimensional search, apart from being slow
by itself.

In contrast, the primal-dual interior point method has
super-linear convergence and only deals with 2n parameters (n=sample size).

The present method is demonstrated by the R file surv_sacr.R for the data 
set of 109 mice in Table 1 of the above paper of Turnbull and Mitchell,
which is the data file RFM109.txt in the directory.

The program expects a data file with three columns of the same type as RFM109.txt.
The first column contains the observations on the times (ties allowed).
The second and third column contain indicators, where the first indicator equals 1
if time of onset of the disease studied is less than or equal to Y (an observation time)
and zero otherwise, and the second indicator equals 1 if the time of death from the
disease studied is less than or equal to Y, and zero otherwise. The program gives back
the maximum likelihood estimates of the two distribution functions. They are output1$MLE1
and output1$MLE2 for the primal-dual algorithm and output2$MLE1
and output2$MLE2 for the EM algorithm. The demonstration program also gives a drawing of
these functions.
The codes for running the computation are in primal_dual.cpp and EM.cpp.

See also: vander Laan, M.J. and Jewell, N.P (2003). Current status and right-censored
data structures when observing a marker at the censoring time. Annals of
Statistics 31, 512–535, in particular pages 514-515.
The picture of the survival functions can be compared with Figure 2 on p. 46 of the
1984 paper of Turnbull and Mitchell (see above).

For comparison, the EM algorihm and its result on the same data is shown. The results
of the two methods are exactly the same for the points where the values are uniquely
determined by the likelihood function, but slightly different between these points,
as can be seen from the plot. A counter of the 10,000 iterations of EM is shown for "keeping faith".
The log likelihoods of the results of the two methods are seen by typing
output1$loglikelihood and output2$loglikelihood. Note that they are the same in the
decimals shown and also equal to the log likelihood reported in Turnbull and Mitchell (1984).

The R script simulation.R gives a simulation, using the primal-dual interior point method, for 1000 samples of size 1000 for the model of Example 1 of the manuscript http://dutiosb.twi.tudelft.nl/%7Epietg/survsac-v2.pdf (2018) of Gomes, Groeneboom and Wellner, which is also studied by van der Laan, Jewell and Peterson (1979), Efficient estimation of the lifetime and disease onset distribution, Biometrika, 84, 539-554. The simulation script gives an example of the computations of Table 3, p. 26, of Gomes, Groeneboom and Wellner (2018).
