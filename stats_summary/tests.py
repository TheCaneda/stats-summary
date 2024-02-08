from thorough_examples import thorough_examples

summaries = {
    "mann-whitney-test": """
    The Mann-Whitney U test is a non-parametric statistical test used to compare TWO INDEPENDENT samples to 
    determine whether there is a difference in their central tendency, specifically their MEDIANS. 
    It is particularly useful in the following contexts:

    * Non-Normal Distribution: When the data do not follow a normal distribution, which violates the assumption
    required for parametric tests like the t-test.

    * Ordinal Data: The test can be applied to ordinal data, where the data can be ranked but the intervals between 
    ranks are not necessarily equal.

    * Unequal Variances: It is also appropriate when the two groups have variances that are not assumed to be equal, 
    another assumption underlying many parametric tests.

    * Small Sample Sizes: The Mann-Whitney U test can be used for small sample sizes, where the central limit theorem might not apply, 
    and thus, normality cannot be assumed.

    The test works by ranking all the observations from both groups together and then comparing the sums of ranks between the groups. 
    The null hypothesis of the Mann-Whitney U test is that the distributions of both groups are equal. 
    If the test finds a significant difference in the sum of ranks, it suggests that one group tends to have higher or lower 
    values than the other, indicating a difference in their distributions.
    """,
    "z-test": """
    The z-test is a statistical test used to determine whether two population means are different when the variances 
    are known and the sample size is large. It is applicable in situations such as:

    * Comparing the mean of a sample to a known population mean (one-sample z-test).
    * Comparing the means of two independent samples (two-sample z-test).
    * Comparing observed proportions to expected proportions in large samples.

    The z-test assumes the distribution of the sample mean is normally distributed. 
    This assumption is valid either by the Central Limit Theorem for large sample sizes or if the original population is 
    normally distributed for any sample size. The test provides a means to test hypotheses about population means or 
    proportions, offering insights into whether observed differences are statistically significant.
    """,
    "t-test": """
    The t-test is a type of inferential statistic used to determine if there is a significant difference between 
    the means of two groups, which may be related in certain features. It is used in the following contexts:

    * Comparing the means of two independent groups (independent samples t-test).
    * Comparing the means of two related groups (paired t-test).
    
    The t-test is designed for situations where the population variance is unknown, making it especially useful for 
    small sample sizes. It assumes that the data follows a normal distribution and that the samples have similar variances. 
    The t-test helps in making decisions about the means of different populations, under the assumption of normally 
    distributed data.
    """,
    "f-test": """
    The F-test is used primarily to compare statistical models that have been fitted to a data set, 
    in order to identify the model that best fits the population from which the data were sampled. Specific uses include:

    * Testing the equality of variances across groups (variance ratio test).
    * Assessing the overall significance of a regression model (ANOVA F-test).

    The F-test evaluates whether the group variances are significantly different from each other, indicating differences in distribution 
    shapes or dispersion. It is critical for ensuring that any observed differences are indeed significant and not due to random chance.
    """,
    "chi-square-test": """
    The chi-square test is a non-parametric statistical test that measures the association between categorical variables. It has two key uses:

    * Testing independence between two variables in a contingency table.
    * Goodness-of-fit testing to see how well an observed distribution matches an expected distribution.

    The chi-square test does not assume a normal distribution for the data, making it a valuable tool for categorical data analysis. 
    It provides a method to gauge the discrepancies between expected and observed frequencies, offering insights into whether any observed 
    differences are statistically significant.""",
    "wilcoxon-test": """
    The Wilcoxon test, also known as the Wilcoxon Signed-Rank Test, is a non-parametric statistical test used to compare two related samples, 
    matched samples, or repeated measurements on a single sample. It assesses whether their population mean ranks differ. It's particularly useful for:

    Data that does not meet the normal distribution assumption.

    Situations where the sample size is too small to rely on the Central Limit Theorem.

    Analyzing ordinal data or non-ordinal data that do not meet the assumptions of parametric tests.

    The test ranks the differences between pairs, ignoring the sign, and then compares the sum of the ranks for positive differences against 
    the sum of the ranks for negative differences. It's an alternative to the paired sample t-test when data do not satisfy the normality assumption.
    """,
    "mann-whitney-test": """
    The Mann-Whitney U test is a non-parametric statistical test used to compare two independent samples to determine whether there is a difference
     in their central tendency, specifically their medians. It is particularly useful in the following contexts:

    * Non-Normal Distribution: When the data do not follow a normal distribution, which violates the assumption required for parametric tests like the t-test.

    * Ordinal Data: The test can be applied to ordinal data, where the data can be ranked but the intervals between ranks are not necessarily equal.

    * Unequal Variances: It is also appropriate when the two groups have variances that are not assumed to be equal, another assumption underlying many parametric tests.

    * Small Sample Sizes: The Mann-Whitney U test can be used for small sample sizes, where the central limit theorem might not apply, and thus, normality cannot be assumed.

    The test works by ranking all the observations from both groups together and then comparing the sums of ranks between the groups. The null hypothesis of
     the Mann-Whitney U test is that the distributions of both groups are equal. If the test finds a significant difference in the sum of ranks, 
     it suggests that one group tends to have higher or lower values than the other, indicating a difference in their distributions.
    """,
    "anova": """
    ANOVA, or Analysis of Variance, is a collection of statistical models used to analyze the differences among group means and their associated procedures, 
    helping to identify if the means of three or more independent samples differ significantly. It's widely used for:

    * Comparing means across multiple groups to determine if at least one group mean differs from the others.
    * Situations where multiple independent variables might influence a single quantitative dependent variable.

    ANOVA assumes that the data within and across groups are normally distributed and that the variances across groups are equal. 
    The test calculates an F-statistic, a ratio of variance between the groups to variance within the groups, to determine whether the observed 
    differences are significant.
    """,
    "kruskal-wallis-test": """
    The Kruskal-Wallis Test is a non-parametric method for testing whether samples originate from the same distribution. 
    It is used for comparing medians across three or more groups, serving as a non-parametric alternative to one-way ANOVA when the data do not meet ANOVA's 
    normality or homogeneity of variance assumptions. Applications include:

    * Situations with ordinal data or when the scale of measurement does not justify the assumption of interval scaling.
    * Data that do not meet the requirements of normal distribution, which is crucial for ANOVA.

    The test ranks all group data together and then evaluates the differences between the ranks of each group. 
    A significant Kruskal-Wallis test suggests that at least one sample median differs from the others.
    """,
    "friedman-test": """
    The Friedman test is a non-parametric alternative to repeated measures ANOVA, used to detect differences in treatments across multiple test attempts. 
    It's applicable for:

    * Comparing three or more related or matched groups.
    * Situations where the data do not meet the assumptions necessary for repeated measures ANOVA, including the normal distribution of residuals.

    The test involves ranking each block of related measurements together and then analyzing these ranks. It's particularly useful for analyzing data from 
    within-subjects designs when parametric test assumptions cannot be satisfied.
    """,
    "correlation": """
    Correlation measures the strength and direction of the linear relationship between two quantitative variables. It's used for:

    * Determining how two variables vary together and understanding the strength of their relationship.
    * Preliminary analyses that inform the development of more detailed statistical models.

    The Pearson correlation coefficient measures linear correlation, while Spearman's rank correlation coefficient and Kendall's tau are used for ordinal 
    data or when the assumption of linearity and normality is not met. These measures provide insights into the predictive relationships between variables.
    """,
    "regression": """
    Regression analysis is a powerful statistical method that allows for the prediction of a dependent variable based on the values of one or more
    independent variables. It's used for:

    * Modeling the relationships between variables and forecasting.
    * Understanding the impact of multiple independent variables on a dependent variable.

    Linear regression deals with linear relationships, while logistic regression is used for binary outcomes, and multiple regression can handle multiple predictors. 
    The analysis provides coefficients that represent the relationship between each independent variable and the dependent variable, 
    facilitating prediction and interpretation of interactions.
    """,
    "time-series": """
    Time series analysis involves statistical techniques for analyzing time order data points to understand the underlying structure and function 
    for forecasting. It's used for:

    * Analyzing trends, seasonal patterns, and cyclical fluctuations in data collected over time.
    * Forecasting future values based on past observations.

    Time series models, like ARIMA and Exponential Smoothing, account for the time-dependent nature of the data, allowing for the analysis of trends, 
    seasonality, and autocorrelation. These models are essential in fields like finance, economics, and weather forecasting, where understanding 
    temporal dynamics is crucial.
    """,
}

code_examples = {
    "mann-whitney-test": """

    # two different teaching methods on two different student groups (traditional x interactive)

    group_a = np.array([ 75, 85, 80, 70, 90, 95, 65, 80 ])
    group_b = np.array([ 85, 90, 75, 88, 80, 84, 82, 95 ])

    full_array = [n for n in group_a] + [n for n in group_b]
    full_array

    ranked_array = scipy.stats.rankdata(full_array)

    ranked_a = ranked_array[:len(group_a)]
    ranked_b = ranked_array[len(group_a):]

    U, p = scipy.stats.mannwhitneyu(group_a, group_b, alternative='two-sided') # use the samples, not the RANKED samples here
    U, p

    """,
    "wilcoxon": """

    # before and after a training program on a group of employees

    before = [ 100, 105, 98, 87, 110, 103, 91, 95, 102, 106 ]
    after = [ 108, 110, 99, 89, 115, 105, 93, 97, 105, 108 ]

    diffs = [ n-m for (n,m) in zip(before, after) ]
    diffs

    scipy.stats.wilcoxon(diffs)

    """,
    "kruskall-wallis": """
    
    # independent samples

    diet_a = [ 3, 2, 4, 5, 2 ]
    diet_b = [ 4, 6, 5, 7, 6 ]
    diet_c = [ 5, 4, 6, 7, 8 ]

    scipy.stats.kruskal(diet_a, diet_b, diet_c)

    """,
    "kolmogorov": """

    reaction_times = [ 5.2, 4.8, 6.1, 5.7, 5.4, 5.9, 4.9, 5.3, 6.2, 5.8 ]

    # scipy.stats.kstest(reaction_times, scipy.stats.norm.cdf, N=len(reaction_times)) #wrong call
    scipy.stats.kstest(reaction_times, 'norm', args=(np.mean(reaction_times), np.std(reaction_times, ddof=1))) #correct call

    """,
    "friedman": """
    
    # samples have relations

    goals_first_tri = [ 3, 2, 4, 5, 2 ]
    goals_second_tri = [ 4, 6, 5, 7, 6 ]
    goals_third_tri = [ 5, 4, 6, 7, 8 ]

    scipy.stats.friedmanchisquare(
        goals_first_tri,
        goals_second_tri,
        goals_third_tri
    )

    """,
    "z-test": """

    |||||||||||||||||||||||||||| PROPORTION |||||||||||||||||||||||||||||||||\

    import numpy as np
    from scipy.stats import norm

    # Sample sizes
    n1, n2 = 1200, 1500

    # Number of successes (purchases)
    x1, x2 = 150, 180

    # Proportions of successes
    p1, p2 = x1 / n1, x2 / n2

    # Pooled proportion
    pooled_p = (x1 + x2) / (n1 + n2)

    # Standard error of the difference in proportion
    std_error = np.sqrt(pooled_p * (1 - pooled_p) * (1/n1 + 1/n2))

    # Z-score
    z_score = (p1 - p2) / std_error

    # P-value (two-tailed test)
    p_value = 2 * (1 - norm.cdf(abs(z_score)))

    # print(f"Z-score: {z_score}")
    # print(f"P-value: {p_value}")

    |||||||||||||||||||||||||||| MEAN 2 IND. SAMPLES  |||||||||||||||||||||||||||||||||||||

    import numpy as np
    from scipy.stats import norm

    # Given sample means and standard deviations
    mean1, std_dev1, n1 = 100, 15, 1200  # Sample 1: mean, standard deviation, and size
    mean2, std_dev2, n2 = 105, 20, 1500  # Sample 2: mean, standard deviation, and size

    # Calculate the standard error of the difference in mean
    std_error_diff = np.sqrt((std_dev1**2 / n1) + (std_dev2**2 / n2))

    # Calculate the Z-score for the difference in means
    z_score = (mean1 - mean2) / std_error_diff

    # Calculate the p-value (two-tailed test)
    p_value = 2 * (1 - norm.cdf(abs(z_score)))

    print(f"Z-score: {z_score}")
    print(f"P-value: {p_value}")

    ||||||||||||||||||||||||||||| POP X SAMPLE ||||||||||||||||||||||||||||||||||||||||||||

    z_score = (sample_mean - population_mean) / (population_std_dev / np.sqrt(sample_size))
    p_value = 2 * (1 - norm.cdf(abs(z_score)))

    """,
    "t-test": """

    |||||||||||||||||||||||||||| MEAN 2 IND. SAMPLES |||||||||||||||||||||||||||||||||||||

    import numpy as np
    from scipy.stats import t

    # Given sample means, standard deviations, and sizes
    mean1, std_dev1, n1 = 100, 15, 30  # Sample 1: mean, standard deviation, and size
    mean2, std_dev2, n2 = 105, 20, 30  # Sample 2: mean, standard deviation, and size

    # Calculate the standard error of the difference in means
    std_error_diff = np.sqrt((std_dev1**2 / n1) + (std_dev2**2 / n2))

    # Calculate the degrees of freedom
    df = min(n1 - 1, n2 - 1)

    # Calculate the t-score for the difference in means
    t_score = (mean1 - mean2) / std_error_diff

    # Calculate the p-value (two-tailed test)
    p_value = 2 * (1 - t.cdf(abs(t_score), df))

    print(f"T-score: {t_score}")
    print(f"P-value: {p_value}")

    |||||||||||||||||||||||||||| MEAN PAIRED SAMPLES |||||||||||||||||||||||||||||||||||||

    import numpy as np
    from scipy.stats import t

    # Data for paired samples
    before = np.array([100, 105, 98, 87, 110])
    after = np.array([108, 110, 99, 89, 115])

    # Calculate the differences
    diffs = after - before

    # Calculate the mean of the differences
    mean_diff = np.mean(diffs)

    # Calculate the standard deviation of the differences
    std_diff = np.std(diffs, ddof=1)

    # Sample size
    n = len(diffs)

    # Degrees of freedom
    df = n - 1

    # Calculate the t-score
    t_score = mean_diff / (std_diff / np.sqrt(n))

    # Calculate the p-value (two-tailed test)
    p_value = 2 * (1 - t.cdf(abs(t_score), df))

    print(f"T-score: {t_score}")
    print(f"P-value: {p_value}")

    ||||||||||||||||||||||||||||| SINGLE SAMPLE VS POPULATION MEAN |||||||||||||||||||||||

    import numpy as np
    from scipy.stats import t

    # Sample data
    sample = np.array([100, 102, 104, 98, 96, 101, 99, 103, 97, 105])

    # Population mean
    pop_mean = 100

    # Sample mean and standard deviation
    sample_mean = np.mean(sample)
    sample_std = np.std(sample, ddof=1)

    # Sample size and degrees of freedom
    n = len(sample)
    df = n - 1

    # Calculate the t-score
    t_score = (sample_mean - pop_mean) / (sample_std / np.sqrt(n))

    # Calculate the p-value (two-tailed test)
    p_value = 2 * (1 - t.cdf(abs(t_score), df))

    print(f"T-score: {t_score}")
    print(f"P-value: {p_value}")

    """
}

tests = {
    "z-test": {
        "use-cases": ["Comparing proportions between two large independent samples", "Testing the difference between a sample mean and a population mean when the population standard deviation is known"],
        "description": "Used for testing hypotheses about proportions or means in large samples when the population standard deviation is known.",
        "examples": ["Comparing the proportion of voters favoring a candidate in two cities", "Comparing a sample mean to a known population mean in quality control"],
        "calculation-process": ["Calculate the test statistic based on sample and population parameters", "Determine the p-value from the z-distribution"],
        "formulas": ["Z = (X̄ - μ) / (σ / √n) for means", "Z = (p̂1 - p̂2) / √P(1-P)(1/n1 + 1/n2) for proportions"],
        "parametric": True,
        "summary": summaries.get('z-test'),
        "code_snippets": code_examples.get('z-test')
    },
    "t-test": {
        "use-cases": ["Comparing the means of two independent samples", "Comparing the mean of a single sample against a known mean", "Comparing the means of two paired samples"],
        "description": "Used to compare the means of two groups or a single group against a standard when the population standard deviation is unknown.",
        "examples": ["Testing the effectiveness of a new teaching method", "Measuring the impact of a drug on blood pressure"],
        "calculation-process": ["Calculate the t-score", "Determine the degrees of freedom", "Find the p-value from the t-distribution"],
        "formulas": ["t = (X̄ - μ) / (s / √n) for one-sample", "t = (X̄1 - X̄2) / √((s1²/n1) + (s2²/n2)) for independent samples"],
        "parametric": True,
        "summary": summaries.get('t-test'),
        "code_snippets": code_examples.get('t-test')
    },
    "f-test": {
        "use-cases": ["Comparing the variances of two populations", "Testing overall significance in regression analysis"],
        "description": "Used to compare the variances between two populations or to test the significance of predictors in a regression model.",
        "examples": ["Determining if two manufacturing processes have different variability", "Assessing the significance of predictors in a linear regression model"],
        "calculation-process": ["Calculate the variance of each group", "Divide the larger variance by the smaller variance to find the F-statistic", "Compare the F-statistic to the F-distribution"],
        "formulas": ["F = s1² / s2²"],
        "parametric": True,
        "summary": summaries.get('f-test')
    },
    "chi-square-test": {
        "use-cases": ["Testing independence between two categorical variables", "Goodness-of-fit testing against a known distribution"],
        "description": "Used to assess whether observed frequencies differ from expected frequencies in categorical data.",
        "examples": ["Analyzing the relationship between gender and voting preference", "Comparing observed dice rolls to a uniform distribution"],
        "calculation-process": ["Calculate the expected frequencies", "Compute the chi-square statistic using observed and expected frequencies", "Find the p-value from the chi-square distribution"],
        "formulas": ["χ² = Σ((O-E)² / E) where O is observed frequency and E is expected frequency"],
        "parametric": False,
        "summary": summaries.get('chi-square-test')
    },
    "wilcoxon-test": {
        "use-cases": ["Comparing two related samples", "Non-parametric alternative to the paired t-test"],
        "description": "A non-parametric test for assessing whether two paired samples come from the same distribution.",
        "examples": ["Evaluating pre- and post-treatment effects in a medical study", "Assessing learning outcomes before and after a training session"],
        "calculation-process": ["Rank the differences between pairs", "Sum the ranks for positive and negative differences", "Calculate the test statistic from the smaller of the two sums"],
        "formulas": ["W = min(W+, W-) where W+ is the sum of positive ranks and W- is the sum of negative ranks"],
        "parametric": False,
        "summary": summaries.get('wilcoxon-test'),
        "code_snippets": code_examples.get('wilcoxon'),
        "thorough_examples": thorough_examples.get('wilcoxon')
    },
    "mann-whitney-test": {
        "use-cases": ["Comparing ranks between two independent samples", "Non-parametric alternative to the independent samples t-test"],
        "description": "A non-parametric test used to determine if there is a statistically significant difference between the medians of two independent samples.",
        "examples": ["Evaluating the effect of two different teaching methods on student performance", "Assessing customer satisfaction levels between two service providers"],
        "calculation-process": ["Rank all observations across both groups", "Sum the ranks for each group", "Calculate the U statistic based on ranks"],
        "formulas": ["U = n1*n2 + (n1*(n1+1)/2) - R1", "where n1 and n2 are the sample sizes, and R1 is the sum of ranks for sample 1"],
        "parametric": False,
        "summary": summaries.get('mann-whitney-test'),
        "code_snippets": code_examples.get("mann-whitney-test"),
        "thorough_examples": thorough_examples.get('mann_whitney')
    },
    "anova": {
        "use-cases": ["Comparing means across three or more groups", "Testing the effect of a categorical variable on a continuous outcome"],
        "description": "Used to determine whether there are any statistically significant differences between the means of three or more independent (or related) groups.",
        "examples": ["Analyzing the impact of diet type (vegan, vegetarian, omnivore) on cholesterol levels", "Evaluating the performance of students across different education methods"],
        "calculation-process": ["Calculate group means and the overall mean", "Compute the between-group and within-group variances", "Calculate the F-statistic and find the p-value"],
        "formulas": ["F = (Between-group variance) / (Within-group variance)"],
        "parametric": True,
        "summary": summaries.get('anova')
    },
    "kruskal-wallis-test": {
        "use-cases": ["Comparing medians across three or more independent groups", "Non-parametric alternative to one-way ANOVA"],
        "description": "A non-parametric method for testing whether samples originate from the same distribution. It is used for comparing two or more independent samples of equal or different sample sizes.",
        "examples": ["Comparing customer satisfaction ratings across multiple store locations", "Assessing the effectiveness of different types of pain relief medication"],
        "calculation-process": ["Rank all data from all groups together", "Calculate the sum of ranks for each group", "Compute the test statistic using ranks"],
        "formulas": ["H = (12 / N(N+1)) * Σ(Ri²/ni) - 3(N+1)", "where N is the total number of observations, Ri is the sum of ranks for group i, and ni is the number of observations in group i"],
        "parametric": False,
        "summary": summaries.get('kruskal-wallis-test'),
        "code_snippets": code_examples.get("kruskall-wallis")
    },
    "friedman-test": {
        "use-cases": ["Comparing three or more paired groups", "Non-parametric alternative to repeated measures ANOVA"],
        "description": "Used for analyzing and comparing matched or paired samples across multiple test conditions. It assesses the differences in treatments across multiple test attempts.",
        "examples": ["Evaluating the performance of algorithms on different datasets", "Assessing the taste preference of a food product across different recipes"],
        "calculation-process": ["Rank each block (or subject) across all conditions", "Sum the ranks for each condition", "Compute the Friedman statistic"],
        "formulas": ["χ² = (12 / k(n+1)) ΣRi² - 3n(k+1)", "where k is the number of conditions, n is the number of blocks, and Ri is the sum of ranks for condition i"],
        "parametric": False,
        "summary": summaries.get('friedman-test'),
        "code_snippets": code_examples.get("friedman"),
        "thorough_examples": thorough_examples.get('friedman')
    },
    "correlation": {
        "use-cases": ["Measuring the strength and direction of association between two continuous variables"],
        "description": "Statistical technique to determine how strongly two variables are related to each other. Includes Pearson (linear relationship) and Spearman/Kendall (monotonic relationship) methods.",
        "examples": ["Exploring the relationship between height and weight", "Studying the association between education level and income"],
        "calculation-process": ["Calculate the covariance between variables", "Normalize by the product of the standard deviations", "Determine the correlation coefficient"],
        "formulas": ["Pearson's r = Σ((xi - x̄)(yi - ȳ)) / (nσxσy)", "Spearman's ρ = 1 - (6 Σd²) / (n(n² - 1))", "where d is the difference between ranks of each observation, n is the number of observations"],
        "parametric": True,
        "summary": summaries.get('correlation')
    },
    "regression": {
        "use-cases": ["Modeling the relationship between a dependent variable and one/more independent variables", "Predicting outcomes based on predictors"],
        "description": "A statistical approach to model and analyze the relationships between dependent and independent variables. It identifies the equation that best predicts the dependent variable.",
        "examples": ["Predicting house prices based on their size, location, and age", "Estimating future sales based on advertising spend"],
        "calculation-process": ["Estimate the model parameters (slopes, intercept) minimizing the error", "Calculate the coefficient of determination (R²)", "Assess the significance of predictors"],
        "formulas": ["y = β0 + β1x1 + ... + βnxn + ε", "where y is the dependent variable, x1,...,xn are independent variables, β0 is the intercept, β1,...,βn are the slopes, and ε is the error term"],
        "parametric": True,
        "summary": summaries.get('regression')
    },
    "time-series": {
        "use-cases": ["Analyzing trends, seasonal patterns, and cyclical fluctuations in data collected over time", "Forecasting future values based on past observations"],
        "description": "Statistical techniques for analyzing time-ordered data points to understand the underlying structure and function for forecasting.",
        "examples": ["Forecasting stock market prices", "Analyzing the trend of monthly sales data"],
        "calculation-process": ["Identify trend, seasonality, and residuals", "Model the time series using appropriate methods (ARIMA, Exponential Smoothing)", "Validate the model with diagnostics checks"],
        "formulas": ["Depends on the specific model: ARIMA(p,d,q), ETS, etc.", "where ARIMA is AutoRegressive Integrated Moving Average, ETS is Exponential Smoothing State Space model"],
        "parametric": True,
        "summary": summaries.get('time-series')
    },
    "kolmogorov-smirnov-test": {
        "use-cases": ["Comparing a sample distribution to a reference probability distribution (one-sample KS test)", "Comparing two sample distributions to determine if they are from the same distribution (two-sample KS test)"],
        "description": "The Kolmogorov-Smirnov (KS) test is a non-parametric test used to determine whether a sample comes from a specific distribution or to compare two samples. It is valuable for testing the goodness-of-fit between observed data and a known distribution or for comparing two empirical distributions without assuming them to follow a specific distribution.",
        "examples": ["Testing if reaction times in a cognitive psychology experiment follow a normal distribution", "Comparing the distribution of daily returns of two different stocks"],
        "calculation-process": ["Calculate the cumulative distribution function (CDF) for the reference distribution or for both samples in the case of the two-sample KS test", "Compute the maximum distance (D-statistic) between the CDFs", "Use the D-statistic to assess the hypothesis through the KS distribution"],
        "formulas": ["D = max|F1(x) - F2(x)| for the two-sample KS test", "D = max|Fn(x) - F(x)| for the one-sample KS test", "where F1(x) and F2(x) are the empirical distribution functions of the two samples, Fn(x) is the empirical distribution function of the sample, and F(x) is the CDF of the reference distribution"],
        "code_snippets": code_examples.get('kolmogorov'),
        "parametric": False,
        "summary": summaries.get('kolmogorov-smirnov-test'),
        "thorough_examples": thorough_examples.get('kolmogorov')
    }
}

