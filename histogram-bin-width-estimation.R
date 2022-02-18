# Find the precision of a floating point number x's
# last significant digit after rounding x to m significant digits.
# This function is useful in combination with the R base round function.
# Ex.) Input of x = 0.01, m = 2 returns 3
# Ex.) Input of x = 0.00001, m = 2 returns 6
# Ex.) Input of x = 11, m = 2 returns 0
# Ex.) Input of x = 10, m = 2 returns 0
# Ex.) Input of x = 100, m = 2 returns -1
get_precision = function(x, m) {
  x = abs(signif(x, m))
  int_part = as.integer(x)
  
  if (int_part == 0) {
    num_zeros_after_decimal_before_sig_figs = as.integer(abs(log10(x))) - 1
    return(m + num_zeros_after_decimal_before_sig_figs)
  }
  else {
    int_part_num_digits = as.integer(log10(int_part) + 1)
    return(m - int_part_num_digits)
  }
}

# Use the precision of the first s significant digits 
# from the range of x as the rounding value.
# output: list. Contains range_data and precision.
get_rounding_info = function(x, s) {
  range_data = max(x) - min(x)
  stopifnot("x has no variability."= !isTRUE(all.equal(range_data, 0)))
  return(list(range_data = range_data, precision = get_precision(range_data, s)))
}

# input:
# x: numeric. Data.
get_bin_widths = function(x) {
  rounding_info = get_rounding_info(x, 1)
  range_data = rounding_info[["range_data"]]
  precision = rounding_info[["precision"]]
  mag = 10^(-precision - 2)
  h = seq(from = mag, to = range_data + mag, by = mag)
  return(h)
}

# input:
# h: numeric. Bin Width.
# output:
# numeric. Vector of 11 possible origins.
get_origin_seq = function(x, h) {
  rounding_value = get_rounding_info(x, 2)[["precision"]]
  # Get possible origins
  min_x = min(x)
  starting_val = round(min_x, rounding_value)
  if (starting_val > min_x) {
    starting_val = starting_val - 0.5 * 10^-rounding_value
  }
  return(starting_val - h * seq(1, 0, -0.1))
}

# For each possible bin width, compute some
# viable origins for a histogram.
# This an important step because histograms can
# look very different depending on the choice of origin!
# input: 
# x: numeric. The numeric data.
# output:
# list. First element is vector of bin widths.
# 2nd element is matrix of possible origins with
# each column corresponding to a different bin width.
get_origins = function(x) {
  h = get_bin_widths(x)

  origins = vapply(X = h, FUN = get_origin_seq, FUN.VALUE = numeric(11L), x = x, USE.NAMES = FALSE)
  
  origins_for_each_bin_width = list(h, origins)
  names(origins_for_each_bin_width) = c("h", "origins")

  return(origins_for_each_bin_width)
}

get_breaks = function(x, origin, h) {
  return(seq(from = origin, to = max(x) + h, by = h))
}

# This is a helper function that is not vectorized.
get_UCV1 = function(x, origin, h){
  # sample size
  n = length(x) 
  
  breaks = get_breaks(x, origin, h)
  
  counts = hist(x = x, breaks = breaks, right = FALSE, plot = FALSE)$counts
  
  UCV = 2/((n - 1)*h) - (n + 1)/(n^2*(n-1)*h) * sum(counts^2)
  
  return(UCV)
}

# Compute the Unbiased Cross-Validation Criterion (UCV).
# Minimize this criterion to find an optimal bin width
# for a histogram. 
# Input:
# x: numeric dataset
# origin: numeric vector. Reference point(s) for first bin.
# origin can be of length 1 or longer. If origin has a length
# greater than 1, then for each element in origin, compute
# the UCV.  origin can also be a matrix where each column
# corresponds to different origins for a fixed bin width.
# h: bin width. numeric vector. 
# Output:
# fix
# numeric vector of same length as the number of bins
# that the chosen bin width would produce, i.e. of
# the number of breaks minus 1 or rather
# ceiling((max(x) + h - origin) / h) - 1

get_UCV = function(x){
  origins_for_each_bin_width = get_origins(x)
  
  UCV = matrix(data = NA_real_, 
               nrow = NROW(origins_for_each_bin_width$origins),
               ncol = NCOL(origins_for_each_bin_width$origins))
  
  for (j in seq_len(NCOL(origins_for_each_bin_width$origins))) {
    for (i in seq_len(NROW(origins_for_each_bin_width$origins))) {
      UCV[i, j] = get_UCV1(x, 
               origins_for_each_bin_width$origins[i, j], 
               origins_for_each_bin_width$h[j])
    }
  }
  return(list(h = origins_for_each_bin_width$h, origins = origins_for_each_bin_width$origins, UCV = UCV))
}

# Given the numeric data x, first find good bin widths to test.
# Then, feed these bin widths into get_UCV to evaluate them.
minimize_UCV = function(x) {
  UCV = get_UCV(x)
  
  inds = arrayInd(which.min(UCV[["UCV"]]), dim(UCV[["UCV"]]))
  
  return(list(origin = UCV[["origins"]][inds], h = UCV[["h"]][inds[, 2]], min_UCV = UCV[["UCV"]][inds]))
}

get_num_bins = function(x, origin, h) {
  return(ceiling((max(x) + h - origin) / h) - 1)
}

# Plot a histogram.
# input:
# x: numeric. Data.
# base_hist: logical. Should the result of the base
# hist call be returned instead?
optimal_hist = function(x, base_hist = FALSE, ...) {
  hist_parameters = c(list(...), list(x = x))

  if (!exists("plot", where = hist_parameters)) {
    hist_parameters$plot = TRUE
  }
  if (hist_parameters$plot == TRUE) {
    if (!exists("xlab", where = hist_parameters)) {
      hist_parameters$xlab = deparse(substitute(x))
    }
    if (!exists("main", where = hist_parameters)) {
      hist_parameters$main = c("Histogram of ", deparse(substitute(x)))
    }
  }

  # Create an optimal histogram.
  if (!base_hist & !exists("nclass", where = hist_parameters) & !exists("breaks", where = hist_parameters)) {
    
    optimal_parameters = minimize_UCV(x)
    
    num_bins = get_num_bins(x, optimal_parameters$origin, optimal_parameters$h)
    
    hist_parameters$breaks = seq(from = optimal_parameters$origin, by = optimal_parameters$h, length.out = num_bins + 1)
    
    opt_hist = do.call(hist, hist_parameters)
    
    opt_hist$xname = deparse(substitute(x))
    
    return(c(UCV = optimal_parameters[["min_UCV"]], opt_hist))
  } 
  # Create a custom, non-optimal histogram.
  else {
    
    custom_hist = do.call(hist, hist_parameters)
    custom_hist$xname = deparse(substitute(x))
    
    # If the bin widths are all the same . . . 
    if (custom_hist$equidist) {
      
      h = custom_hist$breaks[2] - custom_hist$breaks[1]

      return(c(UCV = get_UCV1(x = x, 
                              origin = custom_hist$breaks[1], 
                              h = h), 
               custom_hist))
    } 
    # Else the bin widths are different and we can't
    # get the UCV . . .
    else {
      return(custom_hist)
    }
  }
}
