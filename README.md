# Numerical-Analysis

1. Q1 is an attempt to compute the eigen spectrum of a covariance matrix of returns of an Index Fund. It uses Power Method to successively calculate eigenvalues and their eigenvectors. Tolerance defined as the ratio of eigenvalue to dominant eigenvalue, is used to terminate the program.

Usage python FE_Ass_Q1.py russell_cov.txt 0.01

2. Q2 uses the power method from Q1 on stock prices data. We first perform a missing value treatment using spline interpolation, calculate the returns matrix and its covariance and use Power Method with a tolerance value to output eigen spectrum;

Usage python FE_Ass_Q2.py missing.dat 0.01

3. Q3 analyses the effect of number of time samples used to construct a returns covariance matrix and also the change in eigen spectrum. If there are N assets and T time samples, then N-T projections of the true Covariance matrix cannot be explained (T < N). As T increases, our estimator will converge to a consistent estimator. We also look at change in the eigenvalues (risk associated with eigenportfolios) and the Tangent of the angle of the dominant eigenvectors to analyse how stable the covariance matrix is over time.

Usage python FE_Ass_Q3.py missing.dat

4. Q4 looks at an alternative approach to calcualte the eigen spectrum. We raise the covariance matrix to high power (repeated squaring) so that a random vector would align with most dominant eigenvector. We observe that for large matrices, this approach is slower because of matrix matricx multiplications despite it requiring on an average log(N) iterations where N = number of iterations required by Power Method.
