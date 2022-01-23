#' lm_diag_plots() is a general function for plotting residual diagnostics for an lm() object
#' By: Prof. Erik B. Erhardt, UNM Statistics, erike@stat.unm.edu, StatAcumen.com
#' 02/02/2020
#'
#' @param fit               linear model object returned by lm()
#' @param rc_mfrow          number of rows and columns for the graphic plot, default is c(1, 3); use "NA" for a single plot with 3 columns
#' @param which_plot        default plot numbers for lm()
#' @param n_outliers        number to identify in plots from lm() and qqPlot()
#' @param sw_qqplot         T/F for whether to show the QQ-plot
#' @param sw_boxcox         T/F for whether to show Box-Cox transformation
#' @param sw_constant_var   T/F for whether to assess constant variance
#' @param sw_collinearity   T/F for whether to assess multicollinearity between predictor variables
#' @param sw_order_of_data  T/F for whether to show residuals by order of data
#' @param sw_addedvar       T/F for whether to show added-variables plot
#' @param sw_plot_set       NULL to accept other plot options, or "simple" to exclude boxcox, constant var, collinearity, order of data, and added-variable plots. "simpleAV" to add back in the added-variable plots.  "all" includes all possible plots in this function.
#'
#' @return NULL, invisibly
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ cyl + disp + hp + gear, data = datasets::mtcars)
#' lm_diag_plots(fit)
#' mod <- formula(mpg ~ cyl + disp + hp + gear)
#' fit <- lm(mod, data = datasets::mtcars)
#' lm_diag_plots(fit)
lm_diag_plots <-
  function(
    fit              = NULL
  , rc_mfrow         = c(1, 3)
  , which_plot       = c(4, 6, 1)
  , n_outliers       = 3
  , sw_qqplot        = TRUE
  , sw_boxcox        = TRUE
  , sw_constant_var  = TRUE
  , sw_collinearity  = TRUE
  , sw_order_of_data = TRUE
  , sw_addedvar      = TRUE
  , sw_plot_set      = c(NA, "simple", "simpleAV", "all")[1]
  ) {
  # main version is in ./notes, copy to ./worksheet and ./homework

  ### Function arguments
  ## fit               linear model object returned by lm()
  ## rc_mfrow          number of rows and columns for the graphic plot, default is c(1, 3); use "NA" for a single plot with 3 columns
  ## which_plot        default plot numbers for lm()
  ## n_outliers        number to identify in plots from lm() and qqPlot()
  ## sw_qqplot         T/F for whether to show the QQ-plot
  ## sw_boxcox         T/F for whether to show Box-Cox transformation
  ## sw_constant_var   T/F for whether to assess constant variance
  ## sw_collinearity   T/F for whether to assess multicollinearity between predictor varaibles
  ## sw_order_of_data  T/F for whether to show residuals by order of data
  ## sw_addedvar       T/F for whether to show added-variables plot
  ## sw_plot_set       NULL to accept other plot options, or "simple" to exclude boxcox, constant var, collinearity, order of data, and added-variable plots. "simpleAV" to add back in the added-variable plots.  "all" includes all possible plots in this function.

  ## Changes
  # 2/2/2020 created function
  # 2/14/2020 ask=FALSE to remove prompt "Hit <Return> to see next plot:"

  # sw_plot_set version.  NA accepts argument switch settings.
  if(!is.na(sw_plot_set)) {
    if(sw_plot_set == "all") {
      which_plot = c(4, 6, 1, 2, 3, 5)
    }
    if(sw_plot_set == "simple") {
      sw_boxcox        <- FALSE
      sw_constant_var  <- FALSE
      sw_collinearity  <- FALSE
      sw_order_of_data <- FALSE
      sw_addedvar      <- FALSE
    }
    if(sw_plot_set == "simpleAV") {
      sw_boxcox        <- FALSE
      sw_constant_var  <- FALSE
      sw_collinearity  <- FALSE
      sw_order_of_data <- FALSE
    }
  }

  # variable names
  var_names <- names(fit$model)[-1]
  # display settings
  if (is.na(rc_mfrow[1])) {
    n_plots <-
      sw_qqplot           +
      length(which_plot)  +
      length(var_names)   +
      sw_boxcox           +
      sw_constant_var     +
      sw_collinearity     +
      sw_order_of_data    +
      sw_addedvar
    rc_mfrow <- c(ceiling(n_plots / 3), 3)
    #rc_mfrow <- c(1, 3)
  }
  op <- par(no.readonly = TRUE) # the whole list of settable par
  par(mfrow = rc_mfrow)


  # Normal quantile plot (QQ-plot)
  #library(car)
  if(sw_qqplot) {
    car::qqPlot(as.numeric(fit$residuals), las = 1, id = list(n = n_outliers), main = "QQ Plot", ylab = "Residuals")
  }

  #library(nortest)
  nortest::ad.test(fit$residuals)


  # default: Fitted, Cook's distance (with cutoff), and Leverage (with cutoffs)
  for(i_plot in which_plot) {
    par(ask = FALSE)  # do not ask for next plot
    plot(fit, which = i_plot, id.n = n_outliers)
    if (i_plot == 4) {
      Di_large <- 4 / (dim(fit$model)[1] - dim(fit$model)[2] - 1)
      abline(h = Di_large, col = "blue", lty = 3)  # horizontal line
    }
    if (i_plot == 6) {
      lev_large <- c(2, 3) * dim(fit$model)[2] / dim(fit$model)[1]
      abline(h = Di_large    , col = "blue", lty = 3)  # horizontal line
      abline(v = lev_large[1], col = "blue", lty = 3)  # vertical line
      abline(v = lev_large[2], col = "blue", lty = 2)  # vertical line
    }
  }

  # residuals plotted vs each main effect
  if (length(var_names)) {
    for(i_plot in 1:length(var_names)) {
      m_lab <- paste("Residuals vs.", var_names[i_plot])
      plot(x = fit$model[,var_names[i_plot]], y = fit$residuals, main = m_lab, ylab = "Residuals", xlab = var_names[i_plot])
      abline(h = 0, col = "gray75", lty = 3)  # horizontal line at zero

      if((class(fit$model[,var_names[i_plot]]) %in% c("numeric", "integer"))) {
        # use simple smoother if less than 4 observations, otherwise use splines
        if (length(unique(fit$model[,var_names[i_plot]])) < 4) {
          # Loess
          #lines(predict(loess(fit$residuals ~ fit$model[,var_names[i_plot]], enp.target=1)), col="red", lwd=1)
          # Friedman's SuperSmoother
          lines(
            supsmu(
              x = fit$model[,var_names[i_plot]]
            , y = fit$residuals
            )
          , col = "red"
          , lwd = 1
          )
        } else {
          # if the IQR is 0, we need to increase the tol to be strictly positive
          ss_tol <- 1e-6 * IQR(fit$model[,var_names[i_plot]]) # default
          if (ss_tol == 0) { ss_tol <- 1e-6 * diff(range(fit$model[,var_names[i_plot]])) }
          if (ss_tol == 0) { ss_tol <- 0.1 }
          lines(
            smooth.spline(
              fit$residuals ~ fit$model[,var_names[i_plot]]
            , spar = 0.8
            , tol = ss_tol
            )
          , col = "red"
          , lwd = 1
          , lty = 1
          )
        }

      }
    }
  }

  # residuals vs order of data
  if(sw_order_of_data) {
    # order of data (not always interesting)
    plot(fit$residuals, main = "Residuals vs Order of data", ylab = "Residuals")
    abline(h = 0, col = "gray75", lty = 3)  # horizontal line at zero
  }

  # Box-Cox transformation suggestion
  # only if all values are positive
  if(sw_boxcox) {
    if(min(fit$model[,1] > 0)){
      #library(car)  # car::boxCox relies on family="bcPower", but "bcPower" is a function in the car package
      bcPower <- car::bcPower   # load the required function into the environment
      try(  # this may not work if the model function isn't in the environment, or if other objects are not available, too
        car::boxCox(fit, lambda = seq(-3, 3, length = 101), main = "Box-Cox power transformation")
      )
      rm(bcPower)               # remove it
      abline(v = 1  , col = "orange", lty = 3, lwd = 2)  # horizontal line at zero
    }
  }

  # Evaluate homoscedasticity
  if (length(var_names)) {
    if(sw_constant_var) {
      #library(car)
      # non-constant error variance test
      print(car::ncvTest(fit))
      # plot studentized residuals vs. fitted values
      try(
        car::spreadLevelPlot(fit, sub = "(Homoscedasticity)")
      )
    }
  }

  # Evaluate Collinearity
  if (length(var_names)) {
    if(sw_collinearity) {
      if (length(var_names) >= 2) {
        #library(car)
        vif_val <- car::vif(fit) # variance inflation factors
        dotchart(vif_val, main = "Collinearity", xlab = "Variance Inflation Factor (VIF)", sub = "Not as useful with interactions")
        abline(v = 0  , col = "gray75", lty = 3)  # horizontal line at zero
        abline(v = 2^2, col = "blue"  , lty = 2)  # vertical line
        abline(v = 10 , col = "blue"  , lty = 3)  # vertical line
      } else {
        warning("Collinearity plot only available with at least two predictor (x) variables.")
      }
    }
  }

  # Evaluate Partial regression residual plot (added-variables plot)
  if (length(var_names)) {
    if(sw_addedvar) {
      #library(car)
      car::avPlots(fit, id = list(n = n_outliers), layout = rc_mfrow, ask = FALSE)
    }
  }


  par(op) # reset plotting options

  invisible(NULL)

  ## Useful list of diags: http://www.statmethods.net/stats/rdiagnostics.html
} # end of lm_diag_plots()




#### Visual comparison of whether sampling distribution is close to Normal via Bootstrap
# a function to compare the bootstrap sampling distribution with
#   a normal distribution with mean and SEM estimated from the data
bs_one_samp_dist <- function(dat, N = 1e4) {
  n <- length(dat);
  # resample from data
  sam <- matrix(sample(dat, size = N * n, replace = TRUE), ncol=N);
  # draw a histogram of the means
  sam_mean <- colMeans(sam);
  # save par() settings
  old_par <- par(no.readonly = TRUE)
  # make smaller margins
  par(mfrow=c(2,1), mar=c(3,2,2,1), oma=c(1,1,1,1))
  # Histogram overlaid with kernel density curve
  hist(dat, freq = FALSE, breaks = 6
      , main = "Plot of data with smoothed density curve")
  points(density(dat), type = "l")
  rug(dat)

  hist(sam_mean, freq = FALSE, breaks = 25
      , main = "Bootstrap sampling distribution of the mean"
      , xlab = paste("Data: n =", n
                   , ", mean =", signif(mean(dat), digits = 5)
                   , ", se =", signif(sd(dat)/sqrt(n)), digits = 5))
  # overlay a density curve for the sample means
  points(density(sam_mean), type = "l")
  # overlay a normal distribution, bold and red
  x <- seq(min(sam_mean), max(sam_mean), length = 1000)
  points(x, dnorm(x, mean = mean(dat), sd = sd(dat)/sqrt(n))
       , type = "l", lwd = 2, col = "red")
  # place a rug of points under the plot
  rug(sam_mean)
  # restore par() settings
  par(old_par)
}


#### Visual comparison of whether sampling distribution is close to Normal via Bootstrap
# a function to compare the bootstrap sampling distribution
#   of the difference of means from two samples with
#   a normal distribution with mean and SEM estimated from the data
bs_two_samp_diff_dist <- function(dat1, dat2, N = 1e4) {
  n1 <- length(dat1);
  n2 <- length(dat2);
  # resample from data
  sam1 <- matrix(sample(dat1, size = N * n1, replace = TRUE), ncol=N);
  sam2 <- matrix(sample(dat2, size = N * n2, replace = TRUE), ncol=N);
  # calculate the means and take difference between populations
  sam1_mean <- colMeans(sam1);
  sam2_mean <- colMeans(sam2);
  diff_mean <- sam1_mean - sam2_mean;
  # save par() settings
  old_par <- par(no.readonly = TRUE)
  # make smaller margins
  par(mfrow=c(3,1), mar=c(3,2,2,1), oma=c(1,1,1,1))
  # Histogram overlaid with kernel density curve
  hist(dat1, freq = FALSE, breaks = 6
      , main = paste("Sample 1", "\n"
                    , "n =", n1
                    , ", mean =", signif(mean(dat1), digits = 5)
                    , ", sd =", signif(sd(dat1), digits = 5))
      , xlim = range(c(dat1, dat2)))
  points(density(dat1), type = "l")
  rug(dat1)

  hist(dat2, freq = FALSE, breaks = 6
      , main = paste("Sample 2", "\n"
                    , "n =", n2
                    , ", mean =", signif(mean(dat2), digits = 5)
                    , ", sd =", signif(sd(dat2), digits = 5))
      , xlim = range(c(dat1, dat2)))
  points(density(dat2), type = "l")
  rug(dat2)

  hist(diff_mean, freq = FALSE, breaks = 25
      , main = paste("Bootstrap sampling distribution of the difference in means", "\n"
                   , "mean =", signif(mean(diff_mean), digits = 5)
                   , ", se =", signif(sd(diff_mean), digits = 5)))
  # overlay a density curve for the sample means
  points(density(diff_mean), type = "l")
  # overlay a normal distribution, bold and red
  x <- seq(min(diff_mean), max(diff_mean), length = 1000)
  points(x, dnorm(x, mean = mean(diff_mean), sd = sd(diff_mean))
       , type = "l", lwd = 2, col = "red")
  # place a rug of points under the plot
  rug(diff_mean)
  # restore par() settings
  par(old_par)
}



# Function to plot t-distribution with shaded p-value
t_dist_pval <- function(t_summary) {

  ### Usage
  ## t_summary <- t.test(datasets::mtcars$mpg, mu = 20, data = datasets::mtcars)
  ## t_summary <- t.test(mpg ~ am, mu = 0, data = datasets::mtcars)
  ## t_dist_pval(t_summary)

  par(mfrow=c(1,1))
  lim_extreme <- max(4, abs(t_summary$statistic) + 0.5)
  lim_lower <- -lim_extreme;
  lim_upper <-  lim_extreme;
  x_curve <- seq(lim_lower, lim_upper, length=200)
  y_curve <- dt(x_curve, df = t_summary$parameter)
  plot(x_curve, y_curve, type = "n"
      , ylab = paste("t-dist( df =", signif(t_summary$parameter, 3), ")")
      , xlab = paste("t-stat =", signif(t_summary$statistic, 5)
                   , ", Shaded area is p-value =", signif(t_summary$p.value, 5)))
  if ((t_summary$alternative == "less")
      | (t_summary$alternative == "two.sided")) {
    x_pval_l <- seq(lim_lower, -abs(t_summary$statistic), length=200)
    y_pval_l <- dt(x_pval_l, df = t_summary$parameter)
    polygon(c(lim_lower, x_pval_l, -abs(t_summary$statistic))
          , c(0, y_pval_l, 0), col="gray")
  }
  if ((t_summary$alternative == "greater")
      | (t_summary$alternative == "two.sided")) {
    x_pval_u <- seq(abs(t_summary$statistic), lim_upper, length=200)
    y_pval_u <- dt(x_pval_u, df = t_summary$parameter)
    polygon(c(abs(t_summary$statistic), x_pval_u, lim_upper)
          , c(0, y_pval_u, 0), col="gray")
  }
  points(x_curve, y_curve, type = "l", lwd = 2, col = "black")
  invisible(NULL)
}




# best subset, returns results sorted by BIC
f_bestsubset <-
  function(
    form       # model formula
  , dat        # dataset
  , nbest = 5  # number of models to return for each model size, default is 5
  ) {

### Usage
  ##   # best subset selection on our model
  ##   i_best <-
  ##     f_bestsubset(
  ##       form = formula(mpg ~ cyl + disp + hp + gear)
  ##     , dat  = datasets::mtcars
  ##     )
  ##
  ##     op <- options(); # saving old options
  ##     options(width=100) # setting command window output text width wider
  ##   i_best %>% print(n = Inf, width = Inf)
  ##     options(op); # reset (all) initial options


  # all output in the bs "best subset" object
  bs <-
    leaps::regsubsets(
      form
    , data    = dat
    , nvmax   = 30
    , nbest   = nbest
    , method  = "exhaustive"
    )

  # selected output in named columns
  bs2 <-
    cbind(
      summary(bs)$which   # columns indicating which model terms are included
    , SIZE  = (rowSums(summary(bs)$which) - 1)  # number of terms in model
    , rss   = summary(bs)$rss                   # residual sum of squares
    , r2    = summary(bs)$rsq                   # R^2
    , adjr2 = summary(bs)$adjr2                 # Adj-R^2
    , cp    = summary(bs)$cp                    # Cp
    , bic   = summary(bs)$bic                   # BIC
    ) %>%
    as_tibble() %>%
    # sort models ascending by BIC (best model at top)
    arrange(
      bic
    )

  # return sorted table
  return(bs2)
}



# Graphical Assessment of Multivariate Normality
f_mnv_norm_qqplot <- function(x, name = "") {
  # creates a QQ-plot for assessing multivariate normality

  ### Usage
  ## f_mnv_norm_qqplot(shells[shells$sex == "F", 2:4], "Female")
  ## f_mnv_norm_qqplot(shells[shells$sex == "M", 2:4], "Male")

  x <- as.matrix(x)         # n x p numeric matrix
  center <- colMeans(x)     # centroid
  n <- nrow(x)
  p <- ncol(x)
  cov <- cov(x)
  d <- mahalanobis(x, center, cov) # distances
  qqplot(
      qchisq(ppoints(n), df = p)
    , d
    , main = paste("QQ Plot MV Normality:", name)
    , ylab = "Mahalanobis D2 distance"
    , xlab = "Chi-squared quantiles"
  )
  abline(a = 0, b = 1, col = "red")
}



# Correlation plot with ellipses
f_plot_corr_ellipse <- function(
    corr
  , outline = FALSE
  , col = "grey"
  , upper.panel = c("ellipse", "number", "none")[2]
  , lower.panel = c("ellipse", "number", "none")[1]
  , diag = c("none", "ellipse", "number")[3]
  , digits = 2
  , bty = "n"
  , axes = FALSE
  , xlab = ""
  , ylab = ""
  , asp = 1
  , cex.lab = par("cex.lab")
  , cex = 0.75 * par("cex")
  , mar = c(0, 0, 1, 0)  # 0.1 + c(2, 2, 4, 2)
  , ...
  ) {

  # this is a modified version of the plotcorr function from the ellipse package
  # this prints numbers and ellipses on the same plot but upper.panel and lower.panel changes what is displayed
  # diag now specifies what to put in the diagonal (numbers, ellipses, nothing)
  # digits specifies the number of digits after the . to round to
  # unlike the original, this function will always print x_i by x_i correlation rather than being able to drop it
  # modified by Esteban Buz
  # see http://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/

  if (!require('ellipse', quietly = TRUE, character = TRUE)) {
    stop('Need the ellipse library.  Run: install.packages("ellipse")')
  }
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ #diag behavior
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ #lower half of plot
      if (lower.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { #upper half of plot
      if (upper.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}

# EOF
