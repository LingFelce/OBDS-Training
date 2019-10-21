"""
Introduction to machine learning 21st October 2019

Kaggle - website platform for predictive modelling and analystics competitions - can get scripts from website
scikit-learn - high quality python machine learning library, made machine learning more accessible https://scikit-learn.org/stable/

Supervised learning - regression, classification (know what you want to predict)
Unsupervised learning - clustering, dimensionality reduction (not known how many clusters etc)

Classification - supervised, want to derive set of rules to classify unknown data e.g single cell data where cell type known
  
  SVM - support vector machine - partition space to predict classes (draw lines of best fit to separate data - hyperplanes)
  
  Decision trees - sequence of decisions, go through each variable and choose cutoff, then continue to next variable.
  Train with data to find tree; tends to overfit   (classification only suitable for particular dataset, won't extrapolate other data)
  Can use DecisionTreeClassifier (part of scikit-learn) on dataset in Python, fit data and then use to predict what class data will be in.

Overfitting common problem in machine learning, prediction for new samples worse. Methods to minimise this - for decision trees a very
common variant to minimise overfitting is to use random forests.
Random forests - lots of decisions trees to run on random subsets of data (therefore forest!) - use subsets of variables when splitting
Consensus prediction from forest of trees, reduces overfitting, finds important variables
RandomForestClassifier from sklearn in Python - load and fit, then predict. Works quite quickly, so popular in machine learning.

Other classifiers and performance: https://scikit-learn.org/stable/auto_examples/classification/plot_classifier_comparison.html
Random Forest tends to be the most popular but should pick best method based on input data

Dimensionality reduction - very useful in single cell techniques
Simplified version of data as reduced 20,000 genes x 10,000 cells to 2 dimensions
  
  PCA - linear, get as many new dimensions (principal components) as you start with. New PCs ordered by how much information they have 
  (variance). Discard low information dimensions. Linear combination - can sum data
  
  tSNE - non-linear, more complex manipulation of data. Distance between cells uing Gaussian distribution in 20,000 dimensions.
  Find 2D arrangement that respects those distances. https://distill.pub/2016/misread-tsne/. tSNE not very good at distant neighbours
  
  UMAP - based on algebraic geometry, faster than tSNE, gives good results (distances are better mapped out that tSNE)
  
  Knn - K-nearest neighbour graph - can define trajectories or clusters - phenograph uses community detection. Calculate distance between
  cells, join every cell to its 5 closest cells (most similar)
  Knn can also be used for clustering - different names Louvain clustering and phenograph (Cytof data). Non-linear
  
Machine learning vs stats:
  Machine learning uses similar statistical tools to solved problems, very similar
  Stats - learn about system, moderate size data sets, quantify uncertainty, makes assumptions about data
  ML - focused on accurate prediction and speed, deal with huge datasets, less interested in understanding underlying data
  and how it was generated (fewer checks).
  
"""
