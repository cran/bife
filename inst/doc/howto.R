## ---- eval = FALSE-------------------------------------------------------
#  data(psid, package = "bife")
#  head(psid)

## ---- eval = FALSE-------------------------------------------------------
#  ##    ID LFP KID1 KID2 KID3     INCH AGE TIME
#  ## 1:  1   1    1    1    1 58807.81  26    1
#  ## 2:  1   1    1    0    2 41741.87  27    2
#  ## 3:  1   1    0    1    2 51320.73  28    3
#  ## 4:  1   1    0    1    2 48958.58  29    4
#  ## 5:  1   1    0    1    2 53634.62  30    5
#  ## 6:  1   1    0    0    3 50983.13  31    6

## ---- eval = FALSE-------------------------------------------------------
#  library(bife)
#  stat <- bife(LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID, psid, "probit")
#  summary(stat)

## ---- eval = FALSE-------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.7144874  0.0562414 -12.704   < 2e-16 ***
#  ## KID2      -0.4114824  0.0515524  -7.982  1.44e-15 ***
#  ## KID3      -0.1298719  0.0415477  -3.126   0.00177 **
#  ## log(INCH) -0.2417756  0.0541720  -4.463  8.08e-06 ***
#  ## AGE        0.2319779  0.0375351   6.180  6.40e-10 ***
#  ## I(AGE^2)  -0.0028846  0.0004989  -5.781  7.41e-09 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6058.88,
#  ## null deviance= 8152.05,
#  ## nT= 5976, N= 664
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -1.121

## ---- eval = FALSE-------------------------------------------------------
#  apes_stat <- get_APEs(stat)
#  summary(apes_stat)

## ---- eval = FALSE-------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.0927845  0.0161225  -5.755  8.67e-09 ***
#  ## KID2      -0.0534358  0.0096773  -5.522  3.36e-08 ***
#  ## KID3      -0.0168654  0.0052808  -3.194   0.00140 **
#  ## log(INCH) -0.0313974  0.0109470  -2.868   0.00413 **
#  ## AGE        0.0301250  0.0212329   1.419   0.15596
#  ## I(AGE^2)  -0.0003746  0.0099167  -0.038   0.96987
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE-------------------------------------------------------
#  stat_bc <- bias_corr(stat)
#  summary(stat_bc)

## ---- eval = FALSE-------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.6309003  0.0555075 -11.366   < 2e-16 ***
#  ## KID2      -0.3635502  0.0511328  -7.110  1.16e-12 ***
#  ## KID3      -0.1149806  0.0413488  -2.781   0.00542 **
#  ## log(INCH) -0.2139642  0.0536616  -3.987  6.68e-05 ***
#  ## AGE        0.2052755  0.0373054   5.503  3.74e-08 ***
#  ## I(AGE^2)  -0.0025520  0.0004962  -5.144  2.70e-07 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 6058.88,
#  ## null deviance= 8152.05,
#  ## nT= 5976, N= 664
#  ##
#  ## ( 7173 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -0.969

## ---- eval = FALSE-------------------------------------------------------
#  apes_stat_bc <- get_APEs(stat_bc)
#  summary(apes_stat_bc)

## ---- eval = FALSE-------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## KID1      -0.1015555  0.0146046  -6.954  3.56e-12 ***
#  ## KID2      -0.0585204  0.0090041  -6.499  8.07e-11 ***
#  ## KID3      -0.0185083  0.0052532  -3.523  0.000426 ***
#  ## log(INCH) -0.0344416  0.0102247  -3.368  0.000756 ***
#  ## AGE        0.0330430  0.0186548   1.771  0.076513 .
#  ## I(AGE^2)  -0.0004108  0.0086844  -0.047  0.962272
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## ---- eval = FALSE-------------------------------------------------------
#  library(data.table)
#  setDT(psid)
#  psid[, LLFP := shift(LFP), by = ID]

## ---- eval = FALSE-------------------------------------------------------
#  dyn <- bife(LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) | ID, psid, "probit")
#  dyn_bc <- bias_corr(dyn, L = 1L)
#  summary(dyn_bc)

## ---- eval = FALSE-------------------------------------------------------
#  ## binomial - probit link
#  ##
#  ## LFP ~ LLFP + KID1 + KID2 + KID3 + log(INCH) + AGE + I(AGE^2) |
#  ##     ID
#  ##
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## LLFP       1.0025730  0.0473072  21.193   < 2e-16 ***
#  ## KID1      -0.4741452  0.0679084  -6.982  2.91e-12 ***
#  ## KID2      -0.1958494  0.0625931  -3.129  0.001754 **
#  ## KID3      -0.0754092  0.0505115  -1.493  0.135461
#  ## log(INCH) -0.1947087  0.0621154  -3.135  0.001721 **
#  ## AGE        0.2009782  0.0477735   4.207  2.59e-05 ***
#  ## I(AGE^2)  -0.0024145  0.0006293  -3.837  0.000125 ***
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  ##
#  ## residual deviance= 4774.57,
#  ## null deviance= 6549.14,
#  ## nT= 4792, N= 599
#  ##
#  ## ( 1461 observation(s) deleted due to missingness )
#  ## ( 6896 observation(s) deleted due to perfect classification )
#  ##
#  ## Number of Fisher Scoring Iterations: 6
#  ##
#  ## Average individual fixed effect= -1.939

## ---- eval = FALSE-------------------------------------------------------
#  apes_dyn_bc <- get_APEs(dyn_bc)
#  summary(apes_dyn_bc)

## ---- eval = FALSE-------------------------------------------------------
#  ## Estimates:
#  ##             Estimate Std. error z value Pr(> |z|)
#  ## LLFP       0.1826348  0.0209169   8.731   < 2e-16 ***
#  ## KID1      -0.0752484  0.0114039  -6.598  4.15e-11 ***
#  ## KID2      -0.0310820  0.0074600  -4.166  3.09e-05 ***
#  ## KID3      -0.0119677  0.0053992  -2.217  0.026653 *
#  ## log(INCH) -0.0309009  0.0091310  -3.384  0.000714 ***
#  ## AGE        0.0318959  0.0177688   1.795  0.072646 .
#  ## I(AGE^2)  -0.0003832  0.0079898  -0.048  0.961749
#  ## ---
#  ## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

