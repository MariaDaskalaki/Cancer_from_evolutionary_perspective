#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 20:23:46 2021

@author: maria
"""

import itertools
import collections
import statsmodels.api
from scipy.stats import binom_test
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from sklearn import svm


#pathway="/home/maria/100_cancer_constant1/results"
pathway="/home/maria/no_weights/results"
#pathway="/home/maria/CA_100_cancer_constant1_second_try/results"
nruns=100

dfs_lista=[]

generations_lista=[]
for i in range(1,nruns+1):
    path=pathway+"{}/Frequencies_mutations_CA.txt".format(i)
    df=pd.read_csv(path,sep="\t",header=None)
    generations_lista.append(len(df)-1)
    #if i==2:
     #   break
print (generations_lista)
max_generation=max(generations_lista)
print (max_generation)

for i in range(1,nruns+1):
    path=pathway+"{}/Frequencies_mutations_CA.txt".format(i)
    df=pd.read_csv(path,sep="\t",header=None)
    df.drop(df.columns[-1], axis=1, inplace=True)
    df[0]=df[0].str.split().str.get(1).astype(float)
    selected_df=df
    #print (selected_df)
    transposed_df=selected_df.T
    #print (transposed_df)
    transposed_df.columns = transposed_df.iloc[0]
    transposed_df = transposed_df.reindex(transposed_df.index.drop(0)).reset_index(drop=True)
    transposed_df.columns.name = None
    transposed_df = transposed_df.apply(lambda x: x.str.split(r'\:').str.get(1).astype(float), axis=1)
    #print (transposed_df)
    new_column=pd.Series(transposed_df[transposed_df.columns[-1]])
    if max_generation!=len(selected_df)-1:
        extra_df=pd.concat([new_column]*(max_generation-len(selected_df)+1), axis=1)
        new_big_df=pd.concat([transposed_df, extra_df],ignore_index=True, axis=1)
    else:
        new_big_df=transposed_df
    #print (new_big_df)
    new_big_df['Fixed_freq']=new_big_df[new_big_df.columns[-1]]
    new_big_df['Max_freq'] = new_big_df.max(axis = 1)
    new_big_df['Label']=30*[1]+34*[0]
    dfs_lista.append(new_big_df)
    #if i==2:
     #   break
#print (dfs_lista)
big_df=pd.concat(dfs_lista, ignore_index=True)
print (big_df)

gens_to_keep=[i for i in range(0,max_generation+1,100)]
gens_to_keep.append(max_generation)
gens_to_keep.append('Max_freq')
gens_to_keep.append('Fixed_freq')
gens_to_keep.append('Label')
print (gens_to_keep)


features_df=big_df[gens_to_keep]
print (features_df)


features = list(features_df.columns.values)
x_feat = features_df.loc[:, features].values
#print (x_feat)
y_feat = features_df.loc[:,['Label']].values
#print (y_feat)
x_feat = StandardScaler().fit_transform(x_feat)
#print (x_feat)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x_feat)
#print (principalComponents)

principalDf = pd.DataFrame(data = principalComponents , columns = ['principal component 1', 'principal component 2'])

finalDf = pd.concat([principalDf, features_df[['Label']]], axis = 1)
#print (finalDf)

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0, 1]
colors = ['r', 'g']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['Label'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.show()
print ("The explained variance is",pca.explained_variance_ratio_)


X =features_df.drop('Label', axis=1)
print (X)
y = features_df['Label']
#print (y)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20)
svclassifier = SVC(kernel='poly',probability=True)
svclassifier.fit(X_train, y_train)
y_pred = svclassifier.predict(X_test)
print (y_pred)
print (len(y_pred))
print(confusion_matrix(y_test,y_pred))
print(classification_report(y_test,y_pred))
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
print("Precision:",metrics.precision_score(y_test, y_pred))
print("Recall:",metrics.recall_score(y_test, y_pred))
print (svclassifier.predict_proba(X_test))


cm = confusion_matrix(y_test, y_pred)

plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('SVM Linear Kernel Confusion Matrix')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.show()


def make_meshgrid(x, y, h=.02):
    x_min, x_max = x.min() - 1, x.max() + 1
    y_min, y_max = y.min() - 1, y.max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
    return xx, yy

def plot_contours(ax, clf, xx, yy, **params):
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    out = ax.contourf(xx, yy, Z, **params)
    return out

model = svm.SVC(kernel='poly')
clf = model.fit(principalComponents, y_feat)

fig, ax = plt.subplots()
# title for the plots
title = ('Decision surface of linear SVC ')
# Set-up grid for plotting.
X0, X1 = principalComponents[:, 0], principalComponents[:, 1]
xx, yy = make_meshgrid(X0, X1)

for target, color in zip(targets,colors):
    plot_contours(ax, clf, xx, yy, cmap=plt.cm.coolwarm, alpha=0.8)
    indicesToKeep = finalDf['Label'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
#plot_contours(ax, clf, xx, yy, cmap=plt.cm.coolwarm, alpha=0.8)
#ax.scatter(X0, X1, c=y, cmap=plt.cm.coolwarm, s=20, edgecolors='k')
ax.set_ylabel('PC2')
ax.set_xlabel('PC1')
ax.set_xticks(())
ax.set_yticks(())
ax.set_title('Decison surface using the PCA transformed/projected features')
ax.legend(targets)
plt.show()


