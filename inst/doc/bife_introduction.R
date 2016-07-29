## ---- eval=FALSE---------------------------------------------------------
#  glm(y ~ X + factor(i), family = binomial())

## ---- eval=FALSE---------------------------------------------------------
#  time.bife   <- system.time(bife(  y ~ x + d | id, model = "logit", bias.corr = "ana"))[3]
#  time.clogit <- system.time(clogit(y ~ x + d + strata(id))                          )[3]
#  time.glm    <- system.time(glm(   y ~ x + d + 0 + factor(id), family = binomial()) )[3]

## ---- echo=FALSE, results='asis'-----------------------------------------
# Load 'bife'
library("bife")

# Load results --- store.N and store.T
time.n <- time.n
time.t <- time.t

# N and T vector
N.vector <- rep(100, 10)
T.vector <- rep(10, 10)

# Bind results
results <- cbind("N" = time.n[, 1], "T" = T.vector, time.n[, 2:4], "N" = N.vector, time.t)

# Print results
knitr::kable(results)

## ---- echo=FALSE, warning=FALSE, fig.show='hold'-------------------------
# Load package
library("ggplot2")

# Colour palette for colour-blind
cb.Palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Transform to data.frame
time.n <- data.frame(time.n)
time.t <- data.frame(time.t)

# Regression N
plot.N <- data.frame(N = time.n[["N"]])
plot.N$bife.corr <- fitted(lm(bife.corr ~ N, data = time.n))
plot.N$clogit <- fitted(lm(clogit ~ N, data = time.n))
plot.N$glm<- fitted(lm(glm ~ N + I(N^2) + I(N^3), data = time.n))

# Regression T
plot.T <- data.frame(T = time.t[["T"]])
plot.T$bife.corr <- fitted(lm(bife.corr ~ T, data = time.t))
plot.T$clogit <- fitted(lm(clogit ~ T + I(T^2) + I(T^3), data = time.t))
plot.T$glm<- fitted(lm(glm ~ T, data = time.t))

# Plot N
p <- ggplot(plot.N) +
            ylab(NULL) +
            xlim(100, 1000) +
            ylim(0, 1) +
            theme_bw() +
            theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
            geom_line(aes(N, bife.corr, colour = "bife.corr"), size = 0.5) +
            geom_line(aes(N, clogit, colour = "clogit"), size = 0.5) +
            geom_line(aes(N, glm, colour = "glm"), size = 0.5) +
            scale_color_manual("", values = cb.Palette)

# Plot T
q <- ggplot(plot.T) +
            ylab(NULL) +
            xlim(10, 100) +
            ylim(0, 1) +
            theme_bw() +
            theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
            geom_line(aes(T, bife.corr, colour = "bife.corr"), size = 0.5) +
            geom_line(aes(T, clogit, colour = "clogit"), size = 0.5) +
            geom_line(aes(T, glm, colour = "glm"), size = 0.5) +
            scale_color_manual("", values = cb.Palette)

# Print
p
q

## ---- echo=FALSE, results='asis'-----------------------------------------
# Load results
results.psid <- results.psid

# Change the order of the cols and rename cols
results.psid <- cbind(results.psid[, 1], results.psid[, 3], results.psid[, 2], results.psid[, 4])
colnames(results.psid) <- c("bife",  "glm", "bife.corr", "clogit")

# Print results
knitr::kable(results.psid)

## ---- echo=FALSE, warning=FALSE------------------------------------------
# Load data
psid <- psid

## ------------------------------------------------------------------------
mod.logit <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 | ID, data = psid, bias.corr = "ana")
summary(mod.logit)

## ------------------------------------------------------------------------
apeff.bife(mod.logit, discrete = c("KID1", "KID2", "KID3"), bias.corr = "semi")

## ------------------------------------------------------------------------
mod.probit <- bife(LFP ~ AGE + I(INCH / 1000) + KID1 + KID2 + KID3 | ID, 
                   data = psid, bias.corr = "ana", model = "probit")
summary(mod.probit)

## ------------------------------------------------------------------------
apeff.bife(mod.probit, discrete = c("KID1", "KID2", "KID3"), bias.corr = "semi")

## ---- echo=FALSE, results='asis'-----------------------------------------
# Load results
results.acs <- results.acs

# Change the order of the cols and rename cols
results.acs <- cbind(results.acs[, 1], results.acs[, 3], results.acs[, 2], results.acs[, 4])
colnames(results.acs) <- c("bife",  "glm", "bife.corr", "clogit")

# Print results
knitr::kable(results.acs)

## ---- echo=FALSE, warning=FALSE------------------------------------------
# Load survival
library("survival")

# Load data
acs <- acs

## ------------------------------------------------------------------------
print(try(clogit(LFP ~ AGEP + I(PINCP / 1000) + FER + strata(ST), data = acs)))

## ------------------------------------------------------------------------
mod.logit <- bife(LFP ~ AGEP + I(PINCP / 1000) + FER | ST, data = acs, bias.corr = "no")
summary(mod.logit)

## ------------------------------------------------------------------------
apeff.bife(mod.logit, discrete = "FER")

