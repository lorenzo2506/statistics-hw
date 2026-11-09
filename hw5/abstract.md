HOMEWORK
Es 1 – Riemann and Lebesgue Integration

Choose a collection of interesting functions f1, f2, ... and compute their integrals numerically. Approximate the integrals using approaches inspired by both the Riemann definition (partitioning the domain into subintervals) 
and the Lebesgue definition (partitioning the range of function values into measurable subsets and considering their preimages in the domain).
Compare the two numerical approaches in terms of accuracy, convergence, and computational efficiency. 
For functions that are Riemann-integrable, compare the numerical results and discuss whether the two approaches converge to the same value.
Include examples of functions for which the two approaches behave differently, if possible.

Examples:
f(x) = sin(x), x in [0, pi]
f(x) = x^2, x in [0, 1]
f(x) = 1 if x < 0.5, 3 if x >= 0.5, x in [0, 1]
f(x) = |x - 0.5|, x in [0, 1]
f(x) = 1 / sqrt(x), x in (0, 1]
f(x) = 1 if x < 0.5, 0 if x >= 0.5, x in [0, 1]
f(x) = sin(x) / x, x in [0, pi]

"Pathological" functions:
f(x) = 1 if x is rational, 0 if x is irrational, x in [0, 1] (Dirichlet function)
f(x) = 1/q if x = p/q in lowest terms, 0 if x is irrational, x in [0, 1] (Thomae's function)
...

--------------------------------------------------

Es 2 – Numerical Computation of the Normal CDF

Develop an efficient numerical method to compute the cumulative distribution function (CDF) of the standard normal distribution.
Implement the method and evaluate the CDF at several points, comparing the results with a trusted numerical reference found on the web. Investigate the numerical efficiency and accuracy of your implementation. Explain why the normal CDF cannot be expressed in terms of elementary functions, discussing the mathematical result behind this fact (Liouville's theorem).

Explain how the introduction of special functions, such as the error function (erf), provides an exact analytical representation of the CDF, 
and why numerical methods are still necessary for its practical computation.
Explore possible numerical approaches, such as numerical integration, series expansions, or approximations based on the error function.
https://faculty.etsu.edu/gardnerr/4047/notes-Hogg-McKean-Craig/Hogg-McKean-Craig-3-4.pdf

examples of approaches:
- Numerical integration: using the trapezoidal rule or Simpson's rule.
- Series expansion: using the power series expansion of the error function.
- Error function (erf): using the mathematical relationship between the normal CDF and the error function