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
  
Machine learning exercise on single cell data done in Jupyter Lab - export to HTML.
Classify cell types in human PBMCs using machine learning

"""
!date
!pwd
!conda env list | grep '*'
!python --version

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import umap
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
# from sklearn.gaussian_process import GaussianProcessClassifier - don't bother, takes ages!
from sklearn.metrics import confusion_matrix

#makes plots display nicely in notebook
%matplotlib inline 

# Create labels vector
data_sets = os.listdir("Sep_data")
# Read in all files and create single list (all files joined together)
allcells = [pd.read_csv("Sep_data/"+i) for i in data_sets]
# Change to dataframe - rows are individual cells 10 cell types x 400 cells, columns are genes
df = pd.concat(allcells)

# Use umap with 20 neighbours to get lower dimensions for plotting
# want reduced representation of data
# Create a data frame with umap dimensions 
reducer = umap.UMAP(n_neighbors=20)
embedding = reducer.fit_transform(df)
embedding.shape
# should have 4000 samples, 2 feature columns as 2D

# Create vector containing names of cell types
label_names = np.repeat((data_sets),400)
print(label_names)

# Convert embedding to dataframe, add column describing cell type to table
embedding = pd.DataFrame(embedding)
embedding["cell_type"] = label_names
embedding.shape

# Plot Umap using seaborn scatterplot and colour by cell type
plt.figure(figsize=(10,10))
pbmc_umap = sns.scatterplot(x=0, y=1, hue="cell_type", data=embedding)

# use sklearn.model_selection.train_test_split 
# Split the data and the labels into 0.25 test data
# sklearn.model_selection.train_test_split(*arrays, **options)

df_train, df_test, names_train, names_test = train_test_split(df, label_names, test_size=0.25)

# training is 3000 cells, test is 1000 cells

# Train the Decision Tree Classifier
clf = DecisionTreeClassifier()
# Fit the training data with the classifier. We are using the default parameters
clf = clf.fit(df_train, names_train)

# score accuracy of fitted model
clf.score(df_train, names_train)
# 1.0

# see if model works for test data, score accuracy
clf.score(df_test, names_test)
# 0.682 - model doesn't work as well for test data

# train random forest
clf_forest = RandomForestClassifier() 
clf_forest = clf_forest.fit(df_train, names_train)

# score accuracy of fitted model
clf_forest.score(df_train, names_train)
# 0.997

# score accuracy of model for test data
clf_forest.score(df_test, names_test)
# 0.675 - worse than before!

# change number of trees used
clf_forest = RandomForestClassifier(n_estimators=250) 
clf_forest = clf_forest.fit(df_train, names_train)
# score accuracy of fitted model
clf_forest.score(df_train, names_train)

# score accuracy of model for test data
clf_forest.score(df_test, names_test)
# 0.83 - improved!

# random forest can calculate importance for each feature (genes) it uses to classify
feature_importances = pd.DataFrame(clf_forest.feature_importances_, index=df_train.columns, columns=['importance'])

feature_importances #2 columns with gene names and numerical value for importance, sorted alphabetically by gene name

# sort by importance, descending order
feature_importances.sort_values('importance', ascending=False)
# top genes are markers for cell types

clf_neighbors = KNeighborsClassifier(n_neighbors=1) # default n_neighbors=5; lowering no. of clusters increases accuracy
clf_neighbors = clf_neighbors.fit(df_train, names_train)
# score accuracy of fitting
clf_neighbors.score(df_train, names_train)
# 1.0

# score accuracy of test data
clf_neighbors.score(df_test, names_test)
# 0.703 (not as good as random forest)

# calculate confusion matrix (will be array) - run on test data
predicted_labels = clf_forest.predict(df_test)
cm = confusion_matrix(names_test, predicted_labels)
predicted_labels[:10] # can see is jumbled up compared to actual labels

# convert to pandas data frame to add column and row names
cm_df = pd.DataFrame(cm, index = data_sets, columns=data_sets)
cm_df
# see if predicted data matches up with correct names, should be diagonal line with high values
# if using all data then if matched perfectly should have 400 along diagonal
# algorithm does poorly with naive_cytotoxic and cd4_t_helper (and a few others)

# redo UMAP on df_test
reducer_test = umap.UMAP(n_neighbors=20)
embedding_test = reducer.fit_transform(df_test)
embedding_test.shape

# add names
embedding_test = pd.DataFrame(embedding_test)
embedding_test["cell_type"] = names_test
embedding_test.shape

# Umap for test data
plt.figure(figsize=(20,10))
ax1 = plt.subplot(1, 2, 1) # plot on axis 1; 1 row, 2 columns, 1st plot
ax2 = plt.subplot (1, 2, 2) # plot on axis 2; 1 row, 2 columns, 2nd plot
umap1 = sns.scatterplot(x=0, y=1, hue=names_test, data=embedding_test, ax = ax1, hue_order=data_sets).set_title("True Labels")
umap2 = sns.scatterplot(x=0, y=1, hue=predicted_labels, data=embedding_test, ax = ax2, hue_order=data_sets).set_title("Predicted Labels")

# plot Umaps for total data

all_predicted_labels = clf_forest.predict(df)

# Umap for all data
plt.figure(figsize=(20,10))
ax1 = plt.subplot(1, 2, 1) # plot on axis 1; 1 row, 2 columns, 1st plot
ax2 = plt.subplot (1, 2, 2) # plot on axis 2; 1 row, 2 columns, 2nd plot
umap1 = sns.scatterplot(x=0, y=1, hue="cell_type", data=embedding, ax = ax1, hue_order=data_sets).set_title("True Labels")
umap2 = sns.scatterplot(x=0, y=1, hue=all_predicted_labels, data=embedding, ax = ax2, hue_order=data_sets).set_title("Predicted Labels")

# 3/4 of data used for training, so true v predicted labels will look quite similar!

