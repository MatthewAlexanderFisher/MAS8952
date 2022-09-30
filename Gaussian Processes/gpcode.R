# this file is a collection of all the official code of MAS8952:

dmvn <- function(x, mean, covariance) {
    # computes log density of MVN
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    n <- ncol(x)

    cholesky <- chol(covariance)
    solved_cholesky <- backsolve(cholesky, t(x) - mean, transpose = TRUE)
    quad_form <- colSums(solved_cholesky**2)
    output <- -sum(log(diag(cholesky))) - 0.5 * n * log(2 * pi) - 0.5 * quad_form
    names(output) <- rownames(x)
    return(output)
}


rmvn <- function(N, mean, covariance) {
    # samples from the MVN

    n <- length(mean) # dimension of random vector

    standard_normal_sample <- rnorm(n * N)
    standard_mvn_sample <- matrix(standard_normal_sample, nrow = N)

    cholesky <- chol(covariance)

    output <- t(cholesky) %*% t(standard_mvn_sample) + mean

    return(t(output))
}

condmvn <- function(a, mu1, mu2, K11, K12, K22) {
    # computes conditional of MVN

    inv_K22 <- solve(K22)

    tilde_mu <- mu1 + K12 %*% inv_K22 %*% (a - mu2)
    tilde_K <- K11 - K12 %*% inv_K22 %*% t(K12)

    output <- list("cond_mean" = as.vector(tilde_mu), "cond_covariance" = tilde_K)

    return(output)
}