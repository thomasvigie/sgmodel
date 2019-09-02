# This file solves a stochastic growth  model by vectorization for a given (user-supplied) vector of parameters ---------------------------------



# Package description ----------------------------
#' sgmodel: A package for computating the solutions to a generic stochastic growth model.
#'
#' The sgmodel package provides three important functions:
#' \code{sgmod}, \code{util} and \code{Markovmoments}.
#'
#' @section The \code{sgmodel} function:
#' The \code{sgmodel} function solves a standard stochastic growth model using value function iteration. The stochastic component follows an autoregressive process of order one, and is discretized by a finite state Markov process.
#' @section The \code{util} function:
#'  It computes values for various uility functions encountered in economic theory.
#' @section The \code{Markovmoments} function:
#' It computes the four moments of a finite state Markov chain: expectation, variance, autocovariance and autocorrelation.
#' @import ggplot2
#' @docType package
#' @name package_sgmodel
NULL


# Function that solves the stochastic growth model -----------------------------------------

#' Sgmodel
#'
#' The function \code{sgmodel} computes the solutions to a generic stochastic growth model after discretizing the distribution of the stochastic element.
#' @param grid A numerical value, the number of capital grid points to consider for k (t). Default value set to 1000.
#' @param utiltype The type of preference for the \code{util} function. Can be "log", "CRRA", "CARA", "Cobb-Douglas", "CES". See description of \code{util} for details. Default type set to "log".
#' @param utilparam Numerical value, preference parameter for the \code{util} function. See description of \code{util} for details. Default set to 1.
#' @param A Numerical value, preference parameter for the \code{util} function. See description of \code{util} for details. Default set to 1.
#' @param depre Numerical value for the depreciation parameter. Must be between 0 and 1. Default value set to 1.
#' @param discount Numerical value for the discount factor. Must be (strictly) between 0 and 1. Default value set to 0.95.
#' @param prod Numerical value for the Cobb-Douglas production function. Must be (strictly) between 0 and 1. Default value set to 0.3.
#' @param states umerical value for the number of states of the Markov process approximating the TFP process. Default value set to 2.
#' @param m Numerical value for the \code{Rtauchen} function. See description of \code{Rtauchen} for details. Default value set to 3.
#' @param rho Autocorrelation of the TFP AR(1) process, used to approximate the process with a Markov process.
#' @param sigma Standard deviation of the white noise in the TFP process, used to approximate the process with a Markov process.
#' @param ... Additional arguments.
#' @return The function returns a list containing:
#' \item{Capital grid }{Vector of values for capital.}
#' \item{Savings }{ Vector of size (\code{grid} x \code{States}) indicating which coordinates of the capital grid are the optimal savings decision.}
#' \item{Consumption }{Vector of size (\code{grid} x \code{States}) indicating the optimal consumption decisions using the optimal savings decision, and given the capital level of the corresponding coordinate of \code{Capital grid}.}
#' \item{Z }{States of the TFP process.}
#' \item{PTM }{The probability transition matrix of the process.}
#' \item{Production parameter }{The exponent on capital in the Cobb-Douglas production function.}
#' \item{Utility type }{The type of utility function. See the details of "util" for the available types}
#' \item{Discount factor }{The discount factor used in the model.}
#' \item{Depreciation }{The depreciation rate of capital used in the model.}
#' \item{Rho }{Autocorrelation of the TFP AR(1) process.}
#' \item{Sigma }{Standard deviation of the white noise in the TFP process.}
#' @examples
#' model <- sgmodel(grid= 100, rho = 0.2, sigma = 0.02)
#'
#' grid <- 200
#' utiltype <- "CRRA"
#' utilparam <- 4
#' A <- 1
#' depre <- 0.03
#' discount <- 0.95
#' prod <- 0.3
#' states <- 5
#' m <- 10
#' rho <- 0.2
#' sigma <- 0.02
#' model <- sgmodel(grid, utiltype, utilparam, A, depre, discount, prod, states, m, rho, sigma)
#' @references  Tauchen G (1986), Finite state markov-chain approximations to univariate and vector autoregressions.
#' \emph{Economics letters}, \bold{20}(2), 177--181.
#'
#' Merton R. C (1971), Optimum consumption and portfolio rules in a continuous-time model.
#' \emph{Journal of Economic Theory}, \bold{3}(4), 373--413.
#' URL \url{http://www.sciencedirect.com/science/article/pii/002205317190038X}
#' @export
sgmodel <- function(grid, utiltype, utilparam, A, depre, discount, prod, states, m, rho, sigma, ... ) {

  options(warn = - 1)

  if (missing(A)) {
    A <- 1
    }
  if (missing(grid)) {
    grid <- 1000
    }
  if (missing(utiltype)) {
    utiltype <- "log"
    }
  if (missing(depre)) {
    depre <- 1
    }
  if (missing(discount)) {
    discount <- 0.95
    }
  if (missing(prod)) {
    prod <- 0.3
    }
  if (missing(states)) {
    states <- 2
    }
  if (missing(m)) {
    m <- 3
    }

  if (depre < 0 || depre > 1) {
    stop("depre should be strictly between 0 and 1 !")
    }
  if (discount <= 0 || discount >= 1) {
    stop("discount should be strictly between 0 and 1 !")
    }
  if (prod <= 0 || prod >= 1) {
    stop("prod should be strictly between 0 and 1 !")
    }

  z <- exp( rbind( as.matrix( Rtauchen::Tgrid (states, sigma, rho, m)  ) )  )
  gamma <- Rtauchen::Rtauchen (states, sigma, rho, m)
  TFP <- (1 - discount * (1 - depre)) / (discount * prod)   # makes the non stochastic steady state level of capital equal to 1
  ksslow <- ( ( (discount * prod * TFP * z[1])) / (1 - discount * (1 - depre))) ^ (1 / (1 - prod))
  ksshigh <- ( ( (discount * prod * TFP * z[nrow(z)])) / (1 - discount * (1 - depre))) ^ (1 / (1 - prod))
  a <- log(ksslow) / log(10)
  b <- log(ksshigh) / log(10)
  kgrid <- ramify::logspace(a, b, n = grid) # the grid is spaced in logs. The package 'ramify' contains that function

  rep_col <- function(x, n){
    matrix(rep(x, each = n), ncol = n, byrow = TRUE)
  }

  kcol <- matrix(rep(kgrid, states), grid * states, 1)
  K <- rep_col(kcol, grid)

  rep_row <- function(x, n) {
    matrix(rep(x, each = n), nrow = n)
  }   # functions to repeat a vector row-wise

  Kprime <- data.matrix(rep_row(t(kgrid), states * grid))
  y <- TFP * kgrid ^ (prod) %*% t(z)  #dim = [ngrid . nstates]
  ycol <- matrix(y, grid * states, 1)   #dim = [ ngrid*nstates . 1 ]
  Y <- rep_col(ycol, grid) #dim = [ ngrid*nstates . ngrid ]
  C <- Y + (1 - depre) * K - Kprime
  U <- util(C, A, prefparam = utilparam, type = utiltype, ngoods = 1)   # we allow for different types of utility functions, reported in the util function
  U <- ifelse(C < 0, -Inf, U)   # put a utility of -infinity if consumption is negative
  rm(C)
  V <- matrix(0, nrow = grid, ncol = states)
  v <- matrix(V, nrow = grid * states, ncol = 1 )
  dev <- 1
  iter <- 0
  d <- matrix (0, nrow = grid * states, ncol = 1)

  # Main loop --------------------------------

  e <- matrix(1, nrow = grid, ncol = 1)

  while (dev > 0.001) {
    V <- matrix(v, nrow = grid, ncol = states)
    Vprime <- (gamma %*% t(V)) %x% e
    B <- U + discount * Vprime
    v1 <- apply(B, 1, max)
    dnew <- which(B == rep_col(v1, grid), arr.ind = TRUE)
    dnew <- dnew[order(dnew[, 1] ), 2 ]   # Because the 'which' function reports row in a messy order, I had to reorder them so that the columns picked correspond to each row properly
    dnew <- matrix(dnew, nrow = grid * states, 1)
    dev <- max(abs(d - dnew))
    d <- dnew
    iter <- iter + 1
    v <- v1
  }

  cons <- ycol + ( 1 - depre ) * kcol - data.matrix(kgrid[dnew])
  interest <- prod * kgrid ^ (prod - 1) %*% t(z) - depre
  interest <- matrix(interest, nrow = states * grid, ncol = 1)

  res <- list(kgrid, dnew, cons, z, gamma, prod, utiltype, discount, depre, rho, sigma)
  names(res) <- c("Capital grid", "Savings", "Consumption", "Z", "PTM", "Production parameter", "Utility type", "Discount factor", "Depreciation", "Rho", "Sigma")
  return (res)
}

# Auxiliary functions -----------------------------------------------

#' Markovmoments
#'
#' The function \code{Markovmoments} computes the expectation, variance, autocovariance and autocorrelation of a Markov process.
#' @param states A numerical vector with the states of the Markov process.
#' @param ptm The probability transition matrix, a square matrix of dimension length(states) whose columns sum to one.
#' @param ... Additional arguments.
#' @return It returns a list containing:
#' \item{Expectation }{The mean of the process.}
#' \item{Variance }{The variance of the process.}
#' \item{Autocovariance }{The autocovariance of the process.}
#' \item{Autocorrelation }{The autocorrelation of the process.}
#' \item{Stationary distribution }{The stationary distribution of the process, used for the computation of the moments.}
#' @examples
#' a <- c(-1, 1)
#' A <- matrix(c(0.5, 0.6,
#'               0.5, 0.4), 2, 2)
#' Markovmoments(a, A)
#' @export
Markovmoments <- function(states, ptm, ...) {

  if ( sum(colSums(ptm)) != length(states)) {
    stop("The columns of the probability transition matrix do not sum to one !")
  }

  pi0 <- matrix( 0, nrow = length(states), ncol = 1)
  pi1 <- matrix( 1 / length(states), nrow = length(states), ncol = 1)
  dev <- 0
  while (dev > 0.0001) {
    pi1 <- ptm %*% pi0
    dev <- norm(pi1 - pi0)
    iter <- iter + 1
  }
  statiodist <- pi1
  mu <- crossprod(statiodist, states)
  var <- crossprod(statiodist, states ^ 2) - mu ^ 2
  auto <- states * ptm %*% states
  autocovar <- crossprod(statiodist, auto) - mu ^ 2
  autocorr <- autocovar / var

  res <- list(mu, var, autocovar, autocorr, statiodist)
  names (res) <- c("Expectation", "Variance", "Autocovariance", "Autocorrelation", "Stationary distribution")
  return(res)
}


# Utility function that allows for different types and parameters ---------------------------------------

#' Util
#'
#' The function \code{util} computes values for different types of utility functions and different parameters. See \code{sgmodel_vignette} for detailed functional forms.
#' @param x A numeric vector of length \emph{ngoods} with values to compute utility for.
#' @param A A numerical value that will premultiply the utility function. Default value set to 1.
#' @param prefparam A numerical value, the preference parameter applied to the utility function depending on \emph{type}.
#' @param type A character for the Type of utility function. Can be "log", "CRRA", "CARA", "Cobb-Douglas", "CES". Default type set to "log".
#' @param ngoods Numerical value for the number of goods to consider. Default value set to 1.
#' @param ... Additional arguments.
#' @return A numerical value, the utility function evaluated at the arguments.
#' @examples
#' x <- c(exp(1), exp(1))
#' A <- 2
#' type <- "log"
#' ngoods <- 2
#' util(x = x, A = A, type = type, ngoods = ngoods)
#' @references   Merton R. C (1971), Optimum consumption and portfolio rules in a continuous-time model.
#' \emph{Journal of Economic Theory}, \bold{3}(4), 373--413.
#' URL \url{http://www.sciencedirect.com/science/article/pii/002205317190038X}.
#' @export
util <- function (x, A, prefparam, type = c("log", "CRRA", "CARA", "Cobb-Douglas", "CES" ), ngoods, ... )  {

  if (missing(prefparam)) {
  type <- "log"
  prefparam <- 1
    }
  if (missing (type)) {
    type <- "log"
    }
  if (missing(ngoods)) {
    ngoods <- 1
    }
  if (type == "CES" && missing(ngoods)) {
    ngoods <- 2
    }
  if (type == "CES" && missing(prefparam)) {
    prefparam <- 2
    }
  if (missing(A)) {
    A <- 1
    }
  if (ngoods == 1) {
    if (type == "CRRA")       {
      utility <- A * (x ^ (1 - prefparam)) / (1 - prefparam)
      }
    else if (type == "CARA")   {
      utility <- A * exp( - prefparam * x)
      }
    else if (type == "log")    {
      utility <- A * log(x)
      }
    else if (type == "Cobb-Douglas")    {
      utility <- A * x ^ prefparam
    }
    else if (type == "CES")    {
      utility <- A * x
    }
  }
  if (ngoods > 1) {
    if (type == "CRRA")       {
      utility <- A * sum ( x ^ (as.vector(1 - prefparam) ) / (1 - prefparam) )
      }
    else if (type == "CARA")   {
      utility <- A * sum(exp( - prefparam * x))
      }
    else if (type == "log")    {
      utility <- A * sum(log(x))
      }
    else if (type == "Cobb-Douglas")    {
      utility <- A * prod(x ^ prefparam)
      }
    else if (type == "CES") {
      utility <- A * (sum(x ^ (1 / prefparam))) ^ (prefparam)
      }
  }
  return(utility)
}


#' print_sgmod
#'
#' The function \code{print_sgmod} prints results of the \code{sgmodel} function.
#' @param x A \code{sgmodel} object.
#' @param ... Additional arguments.
#' @return The function prints the call of the function, the \emph{Savings}, \emph{Consumption} and \emph{Capital grid} vectors  from \code{sgmodel}.
#' @examples
#' grid <- 200
#' utiltype <- "CRRA"
#' utilparam <- 4
#' A <- 1
#' depre <- 0.03
#' discount <- 0.95
#' prod <- 0.3
#' states <- 3
#' m <- 5
#' rho <- 0.2
#' sigma <- 0.02
#' model <- sgmodel(grid, utiltype, utilparam, A, depre, discount, prod, states, m, rho, sigma)
#' print_sgmod(model)
#' @export
print_sgmod <- function(x, ...)  {

  cat("Call:\n")
  print(x$call)
  cat("\nStationary stochastic growth model:\n")
  cat("\nSavings decision rule:\n")
  print( (x$"Capital grid"[x$Savings]) )
  cat("\nConsumption decision rule:\n")
  print(x$Consumption)
  cat("\nCapital level:\n")
  print(x$"Capital grid")
}



#' summary_sgmod
#'
#' The function \code{summary_sgmod} prints a summary for results of the \code{sgmodel} function.
#' @param object A \code{sgmodel} object.
#' @param ... Additional arguments.
#' @return It returns a list with the model parameters. It includes:
#' \item{Utility function }{The type of utility function. See the details of \code{util} for the available types}
#' \item{Capital share }{The exponent on capital in the Cobb-Douglas production function.}
#' \item{Discount factor }{The discount factor used in the model.}
#' \item{Depreciation }{The depreciation rate of capital used in the model.}
#' \item{Rho }{Autocorrelation of the TFP AR(1) process.}
#' \item{Sigma }{Standard deviation of the white noise in the TFP process.}
#' \item{Number of TFP states }{Number of states of the TFP process.}
#' @examples
#' grid <- 200
#' utiltype <- "CRRA"
#' utilparam <- 4
#' A <- 1
#' depre <- 0.03
#' discount <- 0.95
#' prod <- 0.3
#' states <- 3
#' m <- 3
#' rho <- 0.2
#' sigma <- 0.02
#' model <- sgmodel(grid, utiltype, utilparam, A, depre, discount, prod, states, m, rho, sigma)
#' summary_sgmod(model)
#' @export
summary_sgmod <- function(object, ...)  {

  descri <- list( object$"Utility type",
               object$"Production parameter",
               object$"Discount factor",
               object$Depreciation,
               object$Rho,
               object$Sigma,
               length(object$Z) )
  names(descri) <- c("Utility function", "Capital share", "Discount factor", "Depreciation", "Rho", "Sigma", "Number of TFP states")
  res <- list(description = descri, call = object$call)
  class(res) <- "summary.sgmod"
  res
}

#' print.summary_sgmod
#'
#' The function \code{print.summary_sgmod} prints a summary for a \code{sgmodel} object.
#' @param x An object of class \code{sgmod}.
#' @param ... Additional arguments.
#' @return It returns a list with the model parameters. It includes:
#' \item{Utility function }{The type of utility function. See the details of \code{util} for the available types}
#' \item{Capital share }{The exponent on capital in the Cobb-Douglas production function.}
#' \item{Discount factor }{The discount factor used in the model.}
#' \item{Depreciation }{The depreciation rate of capital used in the model.}
#' \item{Rho }{Autocorrelation of the TFP AR(1) process.}
#' \item{Sigma }{Standard deviation of the white noise in the TFP process.}
#' \item{Number of TFP states }{Number of states of the TFP process.}
#' @examples
#' grid <- 200
#' utiltype <- "CRRA"
#' utilparam <- 4
#' A <- 1
#' depre <- 0.03
#' discount <- 0.95
#' prod <- 0.3
#' states <- 3
#' m <- 4
#' rho <- 0.2
#' sigma <- 0.02
#' model <- sgmodel(grid, utiltype, utilparam, A, depre, discount, prod, states, m, rho, sigma)
#' summary_sgmod(model)
#' @export
print.summary_sgmod <- function(x, ...)  {

  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$description)
}


#' plot_sgmod
#'
#' The function \code{plot_sgmod} returns a plot of the \code{Savings} value of a \code{sgmodel} object on the \code{Capital grid} value.
#' @param x A \code{sgmod} object.
#' @param ... Additional arguments.
#' @return It returns a plot using \code{ggplot} that graphs the \code{Savings} decisions from the \code{sgmodel} object on the \code{Capital grid}. The plot shows as many facets as \code{length(Z)} where \code{Z} is the vector of states of the TFP process.
#' @examples
#' model <- sgmodel( grid = 100, rho = 0.2, sigma = 0.02)
#' plot_sgmod(model)
#' grid <- 200
#' utiltype <- "CRRA"
#' utilparam <- 4
#' A <- 1
#' depre <- 0.03
#' discount <- 0.95
#' prod <- 0.3
#' states <- 5
#' m <- 2
#' rho <- 0.2
#' sigma <- 0.02
#' model <- sgmodel(grid, utiltype, utilparam, A, depre, discount, prod, states, m, rho, sigma)
#' plot_sgmod(model)
#' @references    Wickham H (2009), ggplot2: Elegant Graphics for Data Analysis.
#' URL \url{http://ggplot2.org}
#' @export
plot_sgmod <- function(x, ...)  {

  nstates <- length(as.vector(x$"Savings")) / length(as.vector(x$"Capital grid"))
  Savings <- x$"Capital grid"[x$Savings]
  dat <- matrix( Savings, nrow = length(x$"Capital grid"), ncol = nstates )
  data <- data.frame(x$"Capital grid", dat)
  for (i in 1:nstates)  {
    names(data)[i + 1] <- paste("Kt1", i, sep = "_")
  }
  names(data)[1] <- "Kt"
  D <- matrix(0, nrow = length(as.vector(x$"Savings")), ncol = 3)
  ind <- 1
  for (i in 1:nstates)  {
    D[ind:(i * length(as.vector(x$"Capital grid"))), ] <- cbind(data[, 1], data[, (i + 1)], rep(i, length(as.vector(x$"Capital grid"))))
    ind <- ind + length(as.vector(x$"Capital grid"))
  }
  data <- data.frame(D)
  names(data) <- c("Kt", "Kt1", "Zt")
  ggplot2::ggplot(data = data) +
    geom_point(mapping = aes_string(x = "Kt", y = "Kt1")) +
    facet_wrap(~ Zt) +
    ylab("Savings decision") + xlab("Capital at time t")
}



