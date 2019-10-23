"""
Linear regression 23rd October 2019

Two quantitative (continuous) variables
  x = predictor, explanatory or independent
  y = target, response, outcome or dependent

Simple linear regression - one x. Multiple linear regression - more than one x.

Correlation (r) measures linear relationship between x and y. 0 = no relationship. -1/+1 neg/pos.

No need to fit data if weakly correlated or perfect straight line already (just work out gradient of straight line)
Only statistical relationships (non-perfect)

Square of correlation coefficient (r2) explains how much variation in y can be explained by x

Common method for linear regression "linear least squares" - must have at least 3 data points
Optimise by minimising sum of differences

In Python numpy or scikit-learn for calculating linear regression

Residuals measure error between data and predicted line - residual plot shows error for each datapoint x
- should be randomly distributed; if not then straight line not suitable

Loss function - non-linear optimisation Levenberg-Marquadt algorithm in Python

Can do simple linear regression in Prism, or can use website for multiple regression https://regression.structural-analyser.com

Limitations - only describes linear relationships, sensitive to outliers, data must be independent

Can do git clone and then add https:// link to pull down git repository (all files and folders) - have to do in git folder/repository
in terminal?

"""
# link to python notebook tutorial (also have exported to HTML) https://github.com/dwaithe/linear_regression_practical
