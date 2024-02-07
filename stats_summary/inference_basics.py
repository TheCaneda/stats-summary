inferential_stats_summary = """
Basic Inferential Statistics Summary
------------------------------------

|| Expected Value (Mean) ||
Definition: The expected value (EV) or mean of a random variable is a fundamental
measure that provides an average outcome of a process in the long run.
How to Calculate:
- For a discrete random variable, EV = Σ[xi * P(xi)], where xi is each possible
  value the variable can assume, and P(xi) is the probability of xi.
- For a continuous random variable, EV is calculated as the integral of x * f(x)
  dx over all possible values of x, where f(x) is the probability density
  function of x.

|| Variance ||
Definition: Variance measures the spread of a set of data points or a random
variable around its mean, indicating how much the values differ from the mean.
How to Calculate:
- For a sample: s² = Σ[(xi - x̄)²] / (n - 1), where xi is each sample value,
  x̄ is the sample mean, and n is the sample size.
- For a population: σ² = Σ[(xi - μ)²] / N, where xi is each population value,
  μ is the population mean, and N is the population size.

|| Standard Deviation ||
Definition: The standard deviation is the square root of the variance, providing
a measure of the dispersion of data points in the same units as the data.
How to Calculate:
- For a sample: s = sqrt(s²), where s² is the sample variance.
- For a population: σ = sqrt(σ²), where σ² is the population variance.

|| Confidence Interval ||
Definition: A confidence interval (CI) provides a range of values within which
the true population parameter (mean, proportion) is expected to lie with a
certain level of confidence.
How to Calculate:
- For the mean with known population variance: CI = x̄ ± Z*(σ/√n), where x̄ is
  the sample mean, Z is the Z-score corresponding to the desired confidence
  level, σ is the population standard deviation, and n is the sample size.
- For the mean with unknown population variance (using t-distribution): CI =
  x̄ ± t*(s/√n), where t is the t-score from the t-distribution for the desired
  confidence level and degrees of freedom (n-1), and s is the sample standard
  deviation.

|| p-Value ||
Definition: The p-value measures the probability of observing the collected data,
or something more extreme, assuming the null hypothesis is true.
How to Calculate: Calculated based on the test statistic (Z-score, t-score, etc.)
and the corresponding distribution. The exact calculation depends on the
statistical test being performed.

These concepts form the backbone of inferential statistics, allowing for the
analysis and interpretation of data, and the making of inferences about a
population based on sample data.
"""
