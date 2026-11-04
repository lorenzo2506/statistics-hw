📚 HMWK

Ex. 1 – Welford’s Recurrence Formulas
Derive Welford’s recurrence formulas for the sum of squares (SS) and the variance, following the same approach as Knuth’s recurrence formula for the mean.
Write the corresponding code. The core logic should require only about 3 lines.

Ex. 2 – Sampling Mean and Variance

Extend the previous homework by assigning each penetration attempt a score between 0 and 10, drawn from a fixed probability distribution of your choice. Let the theoretical mean and variance of this distribution be μ and σ², respectively.
Assume there are m hackers, each attempting to penetrate n systems. Each hacker receives one score for each system, resulting in a trajectory of n scores.
For each hacker, use Welford’s recurrence formulas to compute the sampling mean and sampling variance of their n scores.
Across the m hackers, compute the average of the sampling means and the average of the sampling variances. Compare these two averages with the theoretical mean μ and variance σ², respectively.
What conclusions can you draw from these comparisons? Why are these results useful in inferential statistics? Repeat the experiment for different values of n. How does increasing the number of systems attempted by each hacker affect the variability of the sampling means and sampling variances?