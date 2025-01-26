library(MASS)
library(BayesS5)
library(dlbayes)
library(devtools)
library(lars)
library(ROCR)
library(glmnet)
library(bvartools)
library(CVTuningCov)
library(rstanarm)
install_url('https://cran.r-project.org/src/contrib/Archive/BayesPen/BayesPen_1.0.tar.gz')
require(BayesPen)


# Read the data
data(Boston); attach(Boston)
fit_ols <- lm(medv ~ . - 1, data = Boston)
summary(fit_ols)



# Low Dimension Uncorrelated Variables
data(Boston); attach(Boston)
X = cbind(crim, zn, indus, chas, nox, rm, age, dis, rad, tax, ptratio, black, lstat)
y = medv
y = y - mean(y)
n = nrow(X)
set.seed(123)
## Add independent 20 noise variables generated from Gaussian distribution
X = cbind(X, matrix(rnorm(20*n), n, 20))
## Scale X
X = scale(X)
p = ncol(X)
## MCMC
fit_mcmc <- stan_glm(y ~ X - 1, family = gaussian(link = "identity"), prior = normal(0, 5), iter = 5000, warmup = 500, thin = 10, cores = 4)
fit_mcmc$


## Frequentist Lasso
fit_lasso <- cv.glmnet(X, y, family = "gaussian", alpha = 1, nlambda = 40, intercept = FALSE)
### ROC
total_step <- length(fit_lasso$lambda)
tp_total <- c(); fp_total <- c()
for (i in 1:total_step) {
  mat <- as.matrix(coef(fit_lasso, fit_lasso$lambda[i]))
  tp_total[i] <- sum(mat[2:14,] != 0)
  fp_total[i] <- sum(mat[15:total_step] != 0)
}
fpr_lasso <- sort(fp_total/20)
tpr_lasso <- sort(tp_total/13)
max_step <- 20
sum(tp_total[1:max_step])/max_step ### MS-O
sum(fp_total[1:max_step])/max_step ### MS-N


## Joint
fit_joint <- BayesPen.lm(y, X, prior = list(a1 = 0.1, b1 = 0.1, a2 = 0.1, b2 = 0.1), nIter = 5000, burnIn = 500, joint = TRUE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_joint)
### MS-O
original <- fit_joint$joint.path[1:max_step, 1:13]
tp <- rowSums(original)
(mso <- sum(tp)/max_step)
### MS-N
noise <- fit_joint$joint.path[1:max_step, 14:p]
fp <- rowSums(noise)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_joint$joint.path[, 1:13])
fp_total <- rowSums(fit_joint$joint.path[, 14:p])
fpr_joint <- fp_total/20
tpr_joint <- tp_total/13


## Marginal
fit_marginal <- BayesPen.lm(y, X, prior = list(a1 = 0.1, b1 = 0.1, a2 = 0.1, b2 = 0.1), nIter = 5000, burnIn = 500, joint = FALSE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_marginal)
### MS-O
original <- fit_marginal$coefs[1:max_step, 1:13]
tp <- rowSums(original != 0)
sum(tp)/max_step
### MS-N
noise <- fit_marginal$coefs[1:max_step, 14:p]
fp <- rowSums(noise != 0)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_marginal$coefs[, 1:13] != 0)
fp_total <- rowSums(fit_marginal$coefs[, 14:p] != 0)
fpr_marginal <- fp_total/20
tpr_marginal <- tp_total/13

plot(fpr_joint, tpr_joint, type = "l", lwd = 2, col = "red", main = "ROC Curves in Low Dimension using Independent Variables", xlab = "1 - Specificity", ylab = "Sensiticity")
lines(fpr_marginal, tpr_marginal, type = "l", lty = 2, lwd = 2, col = "blue")
lines(fpr_lasso, tpr_lasso, type = "l", lty = 3, lwd = 2, col = "darkgreen")
legend(0.8, 0.3, legend = c("Joint", "Marginal", "Lasso"), col = c("red", "blue", "darkgreen"), lty = 1:3, cex = 0.8)



# Low Dimension Correlated setting
data(Boston); attach(Boston)
X = cbind(crim, zn, indus, chas, nox, rm, age, dis, rad, tax, ptratio, black, lstat)
y = medv
y = y - mean(y)
n = nrow(X)
set.seed(123)
## Add independent 20 noise variables generated from Gaussian distribution
cov <- as.matrix(AR1(20, rho = 0.5))
X = cbind(X, t(cov%*%t(matrix(rnorm(20*n), n, 20))))
## Scale X
X = scale(X)
p = ncol(X)


## Frequentist Lasso
fit_lasso <- cv.glmnet(X, y, family = "gaussian", alpha = 1, nlambda = 34, intercept = FALSE)
### ROC
total_step <- length(fit_lasso$lambda)
tp_total <- c(); fp_total <- c()
for (i in 1:total_step) {
  mat <- as.matrix(coef(fit_lasso, fit_lasso$lambda[i]))
  tp_total[i] <- sum(mat[2:14,] != 0)
  fp_total[i] <- sum(mat[15:total_step] != 0)
}
fpr_lasso <- sort(fp_total/20)
tpr_lasso <- sort(tp_total/13)
sum(tp_total[1:max_step])/max_step ### MS-O
sum(fp_total[1:max_step])/max_step ### MS-N


## Joint
fit_joint <- BayesPen.lm(y, X, prior = list(a1 = 0.1, b1 = 0.1, a2 = 0.1, b2 = 0.2), nIter = 5000, burnIn = 500, joint = TRUE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_joint)
### MS-O
original <- fit_joint$joint.path[1:max_step, 1:13]
tp <- rowSums(original)
(mso <- sum(tp)/max_step)
### MS-N
noise <- fit_joint$joint.path[1:max_step, 14:p]
fp <- rowSums(noise)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_joint$joint.path[, 1:13])
fp_total <- rowSums(fit_joint$joint.path[, 14:p])
fpr_joint <- fp_total/20
tpr_joint <- tp_total/13


## Marginal
fit_marginal <- BayesPen.lm(y, X, prior = list(c(0.1, 0.1, 0.1, 0.1)), nIter = 5000, burnIn = 500, joint = FALSE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_marginal)
### MS-O
original <- fit_marginal$coefs[1:max_step, 1:13]
tp <- rowSums(original != 0)
sum(tp)/max_step
### MS-N
noise <- fit_marginal$coefs[1:max_step, 14:p]
fp <- rowSums(noise != 0)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_marginal$coefs[, 1:13] != 0)
fp_total <- rowSums(fit_marginal$coefs[, 14:p] != 0)
fpr_marginal <- fp_total/20
tpr_marginal <- tp_total/13

plot(fpr_joint, tpr_joint, type = "l", lwd = 2, col = "red", main = "ROC Curves in Low Dimension using Correlated Variables", xlab = "1 - Specificity", ylab = "Sensiticity")
lines(fpr_marginal, tpr_marginal, type = "l", lty = 2, lwd = 2, col = "blue")
lines(fpr_lasso, tpr_lasso, type = "l", lty = 3, lwd = 2, col = "darkgreen")
legend(0.8, 0.3, legend = c("Joint", "Marginal", "Lasso"), col = c("red", "blue", "darkgreen"), lty = 1:3, cex = 0.8)



# High Dimension
data(Boston); attach(Boston)
X = cbind(crim, zn, indus, chas, nox, rm, age, dis, rad, tax, ptratio, black, lstat)
y = medv
y = y - mean(y)
n = nrow(X)
set.seed(123)
## Add independent 400 noise variables generated from Gaussian distribution
X = cbind(X, matrix(rnorm(400*n), n, 400))
## Scale X
X = scale(X)
p = ncol(X)
## piMoM
fit_piMoM <- S5(X, y)
fit_piMoM$OBJ


## Frequentist Lasso
fit_lasso <- cv.glmnet(X, y, family = "gaussian", alpha = 1, nlambda = 418, intercept = FALSE)
### ROC
total_step <- length(fit_lasso$lambda)
tp_total <- c(); fp_total <- c()
for (i in 1:total_step) {
  mat <- as.matrix(coef(fit_lasso, fit_lasso$lambda[i]))
  tp_total[i] <- sum(mat[2:14,] != 0)
  fp_total[i] <- sum(mat[15:total_step] != 0)
}
fpr_lasso <- sort(fp_total/400)
tpr_lasso <- sort(tp_total/13)
max_step <- 20
sum(tp_total[1:max_step])/max_step ### MS-N
sum(fp_total[1:max_step])/max_step ### MS-N


## Joint
fit_joint <- BayesPen.lm(y, X, prior = list(c(0.1, 0.1, 0.1, 0.1)), nIter = 5000, burnIn = 500, joint = TRUE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_joint)
### MS-O
original <- fit_joint$joint.path[1:max_step, 1:13]
tp <- rowSums(original)
(mso <- sum(tp)/max_step)
### MS-N
noise <- fit_joint$joint.path[1:max_step, 14:p]
fp <- rowSums(noise)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_joint$joint.path[, 1:13])
fp_total <- rowSums(fit_joint$joint.path[, 14:p])
fpr_joint <- fp_total/400
tpr_joint <- tp_total/13


## Marginal
fit_marginal <- BayesPen.lm(y, X, prior = list(c(0.1, 0.1, 0.1, 0.1)), nIter = 5000, burnIn = 500, joint = FALSE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_marginal)
### MS-O
original <- fit_marginal$coefs[1:max_step, 1:13]
tp <- rowSums(original != 0)
sum(tp)/max_step
### MS-N
noise <- fit_marginal$coefs[1:max_step, 14:p]
fp <- rowSums(noise != 0)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_marginal$coefs[, 1:13] != 0)
fp_total <- rowSums(fit_marginal$coefs[, 14:p] != 0)
fpr_marginal <- fp_total/400
tpr_marginal <- tp_total/13

plot(fpr_joint, tpr_joint, type = "l", lwd = 2, col = "red", main = "ROC Curves in High Dimension using Independent Variables", xlab = "1 - Specificity", ylab = "Sensiticity")
lines(fpr_marginal, tpr_marginal, type = "l", lty = 2, lwd = 2, col = "blue")
lines(fpr_lasso, tpr_lasso, type = "l", lty = 3, lwd = 2, col = "darkgreen")
legend(0.8, 0.3, legend = c("Joint", "Marginal", "Lasso"), col = c("red", "blue", "darkgreen"), lty = 1:3, cex = 0.8)



# High Dimension Correlated Variables
data(Boston); attach(Boston)
X = cbind(crim, zn, indus, chas, nox, rm, age, dis, rad, tax, ptratio, black, lstat)
y = medv
y = y - mean(y)
n = nrow(X)
set.seed(123)
## Add independent 400 noise variables generated from Gaussian distribution
cov <- as.matrix(AR1(400, rho = 0.5))
X = cbind(X, t(cov%*%t(matrix(rnorm(400*n), n, 400))))
## Scale X
X = scale(X)
p = ncol(X)


## Frequentist Lasso
fit_lasso <- cv.glmnet(X, y, family = "gaussian", alpha = 1, nlambda = 400, intercept = FALSE)
### ROC
total_step <- length(fit_lasso$lambda)
tp_total <- c(); fp_total <- c()
for (i in 1:total_step) {
  mat <- as.matrix(coef(fit_lasso, fit_lasso$lambda[i]))
  tp_total[i] <- sum(mat[2:14,] != 0)
  fp_total[i] <- sum(mat[15:total_step] != 0)
}
fpr_lasso <- sort(fp_total/400)
tpr_lasso <- sort(tp_total/13)
max_step <- 20
sum(tp_total[1:max_step])/max_step ### MS-O
sum(fp_total[1:max_step])/max_step ### MS-N


## Joint
fit_joint <- BayesPen.lm(y, X, prior = list(c(0.1, 0.1, 0.1, 0.1)), nIter = 5000, burnIn = 500, joint = TRUE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_joint)
### MS-O
original <- fit_joint$joint.path[1:max_step, 1:13]
tp <- rowSums(original)
(mso <- sum(tp)/max_step)
### MS-N
noise <- fit_joint$joint.path[1:max_step, 14:p]
fp <- rowSums(noise)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_joint$joint.path[, 1:13])
fp_total <- rowSums(fit_joint$joint.path[, 14:p])
fpr_joint <- fp_total/400
tpr_joint <- tp_total/13


## Marginal
fit_marginal <- BayesPen.lm(y, X, prior = list(c(0.1, 0.1, 0.1, 0.1)), nIter = 5000, burnIn = 500, joint = FALSE, force = NULL, max.refit = 5000)
BayesPen.plot(fit_marginal)
### MS-O
original <- fit_marginal$coefs[1:max_step, 1:13]
tp <- rowSums(original != 0)
sum(tp)/max_step
### MS-N
noise <- fit_marginal$coefs[1:max_step, 14:p]
fp <- rowSums(noise != 0)
sum(fp)/max_step
### ROC
tp_total <- rowSums(fit_marginal$coefs[, 1:13] != 0)
fp_total <- rowSums(fit_marginal$coefs[, 14:p] != 0)
fpr_marginal <- fp_total/400
tpr_marginal <- tp_total/13

plot(fpr_joint, tpr_joint, type = "l", lwd = 2, col = "red", main = "ROC Curves in High Dimension using Correlated Variables", xlab = "1 - Specificity", ylab = "Sensiticity")
lines(fpr_marginal, tpr_marginal, type = "l", lty = 2, lwd = 2, col = "blue")
lines(fpr_lasso, tpr_lasso, type = "l", lty = 3, lwd = 2, col = "darkgreen")
legend(0.8, 0.3, legend = c("Joint", "Marginal", "Lasso"), col = c("red", "blue", "darkgreen"), lty = 1:3, cex = 0.8)