from . import thorough_examples
from .code_snippets.t_test import code_snippet as t_test_code_snippets
from .code_snippet.z_test import code_snippet as z_test_code_snippets
from .code_snippet.friedman import code_snippet as friedman_code_snippets
from .code_snippet.wilcoxon import code_snippet as wilcoxon_code_snippets
from .code_snippets.mann_whitney import code_snippet as mann_whitney_code_snippets
from .code_snippets.kruskall import code_snippet as kruskall_code_snippets
from .code_snippets.kolmogorov import code_snippet as kolmogorov_code_snippets

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
    "mann-whitney-test": mann_whitney_code_snippets,
    "wilcoxon": wilcoxon_code_snippets,
    "kruskall-wallis": kruskall_code_snippets,
    "kolmogorov": kolmogorov_code_snippets,
    "friedman": friedman_code_snippets,
    "z-test": z_test_code_snippets,
    "t-test": t_test_code_snippets
}

examples = {
    'z-test': [
        """
        \033[1mOne-sample z-test Example:\033[0m
        - \033[1mScenario:\033[0m Testing if light bulbs last 1200 hours on average. A sample of 50 bulbs
        has an average lifespan of 1180 hours.
        - \033[1mFormula:\033[0m z = (x̄ - μ) / (σ / √n)
        - \033[1mFormula Application:\033[0m z = (1180 - 1200) / (100 / sqrt(50)) = -1.41
        - \033[1mInterpretation:\033[0m Since z = -1.41, which is within the critical value of ±1.96,
        there's insufficient evidence to reject the claim that the bulbs last 1200 hours
        on average.

        \033[1mTwo-sample z-test Example:\033[0m
        - \033[1mScenario:\033[0m Comparing average test scores between two classes. Class A: avg=78, n=35.
        Class B: avg=82, n=40.
        - \033[1mFormula:\033[0m z = (x̄1 - x̄2) / √(σ1²/n1 + σ2²/n2)
        - \033[1mFormula Application:\033[0m z = (78 - 82) / sqrt(10^2/35 + 12^2/40) = -1.57
        - \033[1mInterpretation:\033[0m With z = -1.57, there's insufficient evidence to conclude a
        significant difference between the class scores.

        \033[1mProportion z-test Example:\033[0m
        - \033[1mScenario:\033[0m Testing if a website redesign increased the purchase rate from 15% to
        20% (40 out of 200 visitors).
        - \033[1mFormula:\033[0m z = (p̂ - p₀) / √(p₀(1-p₀)/n)
        - \033[1mFormula Application:\033[0m z = (0.20 - 0.15) / sqrt(0.15(1-0.15)/200) = 1.98
        - \033[1mInterpretation:\033[0m Since z = 1.98, which is slightly above 1.96, there's evidence
        suggesting a significant increase in the purchase rate after the website redesign.
        """
    ],
    't-test': [
        """
        \033[1mOne-sample t-test Example:\033[0m
        - \033[1mScenario:\033[0m Testing if the average height of a sample of students is different from the national average height of 170 cm.
        - \033[1mFormula:\033[0m t = (x̄ - μ) / (s / √n)
        - \033[1mFormula Application:\033[0m Hypothetically, t = (168 - 170) / (10 / sqrt(30)) = -1.09
        - \033[1mInterpretation:\033[0m Since t is within the critical value range, there's insufficient evidence to reject the null hypothesis that the sample mean is the same as the national average.

        \033[1mTwo-sample (independent) t-test Example:\033[0m
        - \033[1mScenario:\033[0m Comparing the average test scores of students from two different schools.
        - \033[1mFormula:\033[0m t = (x̄1 - x̄2) / √(s1²/n1 + s2²/n2)
        - \033[1mFormula Application:\033[0m Hypothetically, t = (75 - 80) / sqrt((15^2/50) + (20^2/50)) = -1.89
        - \033[1mInterpretation:\033[0m With t = -1.89, there's insufficient evidence to conclude a significant difference in average test scores.

        \033[1mPaired-sample t-test Example:\033[0m
        - \033[1mScenario:\033[0m Testing the effect of a study app on the scores of students by comparing their scores before and after using the app.
        - \033[1mFormula:\033[0m t = (d̄ - μd) / (sd / √n), where d̄ is the mean difference, and μd is the hypothesized mean difference (often 0).
        - \033[1mFormula Application:\033[0m Hypothetically, t = (5 - 0) / (10 / sqrt(30)) = 3.16
        - \033[1mInterpretation:\033[0m Since t = 3.16, which is beyond the critical value, there's evidence suggesting a significant effect of the study app on scores.
        """
    ],
    'wilcoxon': [
        """
        \033[1mWilcoxon Signed-Rank Test Example:\033[0m
        - \033[1mScenario:\033[0m Testing whether there's a significant difference in the median daily
        calorie intake before and after following a specific diet plan for a group of individuals.
        - \033[1mFormula:\033[0m The test statistic is W, which is the sum of the ranks of the positive
        differences between pairs. The calculation involves ranking the absolute differences, assigning
        signs based on the direction of the difference, and then summing the ranks for the positive differences.
        - \033[1mFormula Application:\033[0m Hypothetically, if the sum of ranks for the positive differences
        (after - before) is 120 and the number of pairs is 30, we would consult a Wilcoxon signed-rank
        table or use software to determine the significance based on W = 120.
        - \033[1mInterpretation:\033[0m Depending on the critical value for W from the Wilcoxon signed-rank
        table for n = 30 and a chosen significance level (e.g., α = 0.05), we determine if there's a
        significant difference in median calorie intake. If W is less than the critical value, we reject
        the null hypothesis, indicating a significant difference in medians before and after the diet.
        """
    ],
    'mann-whitney-test': [
        """
        \033[1mMann-Whitney U Test Example:\033[0m
        - \033[1mScenario:\033[0m Comparing the median waiting times at two different bus stops.
        - \033[1mFormula:\033[0m The test statistic is U, which is calculated based on the ranks of the
        observations from both groups. The U statistic is used to determine the significance of the
        difference in medians between the two groups.
        - \033[1mFormula Application:\033[0m Hypothetically, if the U statistic is 100 and the sample sizes
            are 40 and 50, we would consult a Mann-Whitney U table or use software to determine the
            significance based on U = 100.
        - \033[1mInterpretation:\033[0m Depending on the critical value for U from the Mann-Whitney U table
            for the given sample sizes and a chosen significance level (e.g., α = 0.05), we determine if there's
            a significant difference in median waiting times between the two bus stops.
        """
    ],
    'f-test': [
        """
        \033[1mF-Test Example:\033[0m
        - \033[1mScenario:\033[0m Comparing the variances of two different manufacturing processes.
        - \033[1mFormula:\033[0m The F-statistic is calculated as the ratio of the variances of the two groups.
            The F-test is used to determine if the variances are significantly different from each other.
        - \033[1mFormula Application:\033[0m Hypothetically, if the F-statistic is 1.5, we would consult an
            F-distribution table or use software to determine the significance based on F = 1.5.
        - \033[1mInterpretation:\033[0m Depending on the critical value for F from the F-distribution table for
            the given degrees of freedom and a chosen significance level (e.g., α = 0.05), we determine if there's
            a significant difference in variances between the two manufacturing processes.
        """
    ],
}

tests = {
    "z-test": {
        "use-cases": ["Comparing proportions between two large independent samples", "Testing the difference between a sample mean and a population mean when the population standard deviation is known"],
        "description": "Used for testing hypotheses about proportions or means in large samples when the population standard deviation is known.",
        "examples": examples.get('z-test'),
        "calculation-process": ["Calculate the test statistic based on sample and population parameters", "Determine the p-value from the z-distribution"],
        "formulas": ["Z = (X̄ - μ) / (σ / √n) for means", "Z = (p̂1 - p̂2) / √P(1-P)(1/n1 + 1/n2) for proportions"],
        "parametric": True,
        "summary": summaries.get('z-test'),
        "code_snippets": code_examples.get('z-test')
    },
    "t-test": {
        "use-cases": ["Comparing the means of two independent samples", "Comparing the mean of a single sample against a known mean", "Comparing the means of two paired samples"],
        "description": "Used to compare the means of two groups or a single group against a standard when the population standard deviation is unknown.",
        "examples": examples.get('z-test'),
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
        "examples": examples.get('wilcoxon'),
        "calculation-process": ["Rank the differences between pairs", "Sum the ranks for positive and negative differences", "Calculate the test statistic from the smaller of the two sums"],
        "formulas": ["W = min(W+, W-) where W+ is the sum of positive ranks and W- is the sum of negative ranks"],
        "parametric": False,
        "summary": summaries.get('wilcoxon-test'),
        "code_snippets": code_examples.get('wilcoxon'),
        "thorough_examples": thorough_examples.examples.get('wilcoxon')
    },
    "mann-whitney-test": {
        "use-cases": ["Comparing ranks between two independent samples", "Non-parametric alternative to the independent samples t-test"],
        "description": "A non-parametric test used to determine if there is a statistically significant difference between the medians of two independent samples.",
        "examples": examples.get('mann-whitney-test'),
        "calculation-process": ["Rank all observations across both groups", "Sum the ranks for each group", "Calculate the U statistic based on ranks"],
        "formulas": ["U = n1*n2 + (n1*(n1+1)/2) - R1", "where n1 and n2 are the sample sizes, and R1 is the sum of ranks for sample 1"],
        "parametric": False,
        "summary": summaries.get('mann-whitney-test'),
        "code_snippets": code_examples.get("mann-whitney-test"),
        "thorough_examples": thorough_examples.examples.get('mann_whitney')
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
        "thorough_examples": thorough_examples.examples.get('friedman')
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
        "thorough_examples": thorough_examples.examples.get('kolmogorov')
    }
}

