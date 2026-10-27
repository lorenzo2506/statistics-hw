

HMWK  ==>  Consider n systems and m attackers. Model the penetration process using a stochastic process with independent and identically distributed increments X₁, X₂, …, representing the penetration score at each step. Consider two possible models for the score: 
Bernoulli model: Xᵢ ∈ {0,1}, with P(Xᵢ = 1) = p.
Rademacher model: Xᵢ ∈ {−1,+1}, with P(Xᵢ = +1) = P(Xᵢ = −1) = 1/2. This second model is introduced to study the symmetric, zero-mean case.
For increasing values of n and suitable choices of m, simulate multiple sample trajectories and study the distributions of the following three quantities:
Total score: Sₙ = ∑ᵢ₌₁ⁿ Xᵢ
Average score: Sₙ/n = (1/n)∑ᵢ₌₁ⁿ Xᵢ
Normalized total score: Sₙ/√n = (1/√n)∑ᵢ₌₁ⁿ Xᵢ
Compare the behavior of these quantities as n increases. In particular, investigate the convergence in probability of Sₙ/n according to the Law of Large Numbers, the behavior of the variance under the different normalizations, and the convergence in distribution of Sₙ/√n to a standard normal distribution in the Rademacher case, according to the Central Limit Theorem. Interpret the simulated trajectories in the context of the n systems and m attackers, discussing how the system size and number of attackers affect the aggregate penetration score.

