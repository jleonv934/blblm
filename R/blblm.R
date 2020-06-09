#' @import purrr
#' @import stats
#' @import furrr
#' @import future
#' @importFrom magrittr %>%
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' @title blblm
#' @description Package allows users to perform bootstrap linear regression in a timely manner. This is acheived through the use of C functions to perform linear regression. Parallel Processing is also avaliable for use and the number of cores that are used can be specified as a parameter of the blblm function.
#' @param formula
#' This function performs bootstraped lineear regression
#'
#' @param data data is a data frame to perform lm
#' @param m    number of subsets made from data
#' @param B   number of bootstraps
#' @param n_cores specify number of cores if parallelization is used.
#'
#'
#' @export
blblm <- function(formula, data, m = 10, B = 5000, n_cores = 1) {
  data_list <- split_data(data, m)
  # if else to initiate parallel process with specified cores
  if (n_cores != 1) {
    plan(multiprocess, workers = n_cores)
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
    )
  }
  else {
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B)
    )
  }

  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#'
#' @param data data is entered and subset
#' @param m number of subsets
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#'
#' @param formula formula to fit the bootstraped model
#' @param data subset data
#' @param n number of obs in data
#' @param B number of bootstraps
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#'
#' @param formula formula to fit the bootstraped model
#' @param data subset data
#' @param n number of obs in data
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula formula to fit the bootstraped model
#' @param data subset data
#' @param freqs frequency
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  # Creating model matrix
  X <- model.matrix(reformulate(attr(terms(formula), "term.labels")), data)
  # extracting response variable from data
  y <- as.matrix(data[, all.vars(formula)[1]])
  #T ransforming data to use weights from freq
  w <- freqs
  rw <- sqrt(w)
  rw <- as.vector(rw)
  X_tilde <- rw * X
  y_tilde <- rw * y
  # fitting fast lm model from Rcpp Armadillo
  fit <- fast_lm(X_tilde, y_tilde)
  list(coef = blbcoef(fit, formula), sigma = blbsigma(fit))
}



#' blbfcoef
#'
#' @param fit individual lm model fit with cpp armadillo function
#' @param formula formula to fit the bootstraped model
#'
#' @return coef
blbcoef <- function(fit, formula) {
  cf <- coef(fit) # extracting coefficents
  S <- attr(terms(formula), "term.labels") #extracting variable names from formula
  rownames(cf) <- c("(Intercept)", S) #naming coefficients
  cf
}


#' compute sigma from fit
#'
#' @param fit individual lm model
blbsigma <- function(fit) {
  sqrt(sum(fit$res^2) / fit$df.residual)
}

#' @title printblblm
#' @param x lm bootstrap model taken as a parameter
#'
#' @param ... more parameters
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @title sigmablblm
#' @param object lm bootstrap model taken as a parameter
#'
#' @param confidence specify weather CI should be reurned if so confidence = TRUE
#' @param level insert level to compute CI
#' @param ... more parameters
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - level
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}


#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  results <- as.matrix(t(map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())))
  colnames(results) <- c("(Intercept)", attr(terms(object$formula), "term.labels"))
  results
}




#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  ps <- seq(parm) + 1

  mybiglist <- list()
  for (p in ps) {
    name <- parm[p - 1]
    tmp <- map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2))) #
    mybiglist[[name]] <- tmp
  }
  mybiglist
}


#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
