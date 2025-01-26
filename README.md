## Summary of the Project
For this project, I will be reviewing “Consistent high-dimensional Bayesian variable selection via penalized credible regions” (Bondell and Reich, 2012) and using “Priors for Bayesian shrinkage and high-dimension model selection” (Shin, 2017) as a supplementary
source. This paper studies a new Bayesian variable selection method in linear regression
models called penalized credible regions approach. The second paper provides some instruc-
tions on the dataset that I will be using and lists some methods of performance measurement.
I will be extensively summarizing the selected paper into four parts: introduction, penalized
credible regions, asymptotic theory, and simulations. Next, I will be implementing two types
of the proposed approach, the joint sets and the marginal sets approach, on Boston housing
dataset in varying settings and compare the results with those of the traditional Lasso.
Introduction of the Project
In the second half of the semester, we spent a great amount of time on Bayesian
linear regression and variable selection, which I am very interested in. With the increasing
dimensions of modern data and necessity to reduce dimensions, variable selection is a topic
gaining popularity. Statisticians have designed a variety of methods for variable selection,
including both frequentist and Bayesian approaches. Some frequentist methods include for-
ward selection, backward elimination, least absolute shrinkage and selection operator, fused
Lasso, adaptive Lasso and so on. On the other hand, the Bayesian statistics most frequently
approaches the variable selection problems with a probability-based MCMC sampling call
Stochastic Search Variable Selection (SSVS), which we focused a lot in class and assign-
ments. There are many papers discussing the pros and cons of each method with respect
to accuracy, computation time, whether predictors are correlated or not, dimension relative
to sample size and so on. Finding the most suitable variable selection method for a specific
dataset is always a debatable topic. The paper I selected is the same. This paper suggests
that the penalized credible regions approach is a more efficient and more accurate method
than other commonly used Bayesian and frequentist approach, and it exhibits selection con-
sistency in high dimensional settings. Therefore, I choose to review this paper to see how
the proposed approach differs from other else. The Boston housing dataset is then selected
because it is frequently used to study variable selection.
