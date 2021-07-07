#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 00:26:33 2021

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
from sklearn.metrics import f1_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RandomizedSearchCV
import imblearn
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from sklearn.metrics import accuracy_score, recall_score
import seaborn as sns
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler



################## binomial test for recurrent mutations ############################

def read_fixation_probabilities(path, nruns):
    fixed_positions=[]
    for i in range(1,nruns+1):
        text=path+"{}".format(i)+"/Positions_that_matter.txt"

        with open (text,"r") as f:
            last_line=f.readlines()[-1].split(":")[1].split("\t")[:-1]
            fixed_positions.append(last_line)

    fixed_positions_flatten=list(itertools.chain(*fixed_positions))
    counter=dict(collections.Counter(fixed_positions_flatten))

    fix_probs={k : v/(nruns) for k,v in counter.items()}

    return fix_probs, counter ############### returns a dictionary position:fixation probability accross 100 runs.

neutral_probs=read_fixation_probabilities("/home/maria/Diplomatiki_parousiasi/machine_learning_neutral_new/results", 100)[0]
cancer_probs=read_fixation_probabilities("/home/maria/Diplomatiki_parousiasi/machine_learning_cancerous_new/results", 100)[1]
print ("The neutral probs are")
print (neutral_probs)
print ("The cancer counts are")
print (cancer_probs)

merged_dict = [neutral_probs, cancer_probs] ########## list of cancer and neutral dicitionaries
print ("The merged dict is")
print (merged_dict)

probs_dict = {}
for k in neutral_probs.keys():
    probs_dict[k] = tuple(probs_dict[k] for probs_dict in merged_dict)

print ("The probs dict is")    
print (probs_dict) ########### dictionary position: (neutral prob, cacner count of number of fixation accross 100 runs)

p_values={}
for k,v in probs_dict.items():
    test=binom_test(x=v[1], n=100, p=v[0], alternative='greater')
    p_values[int(k)]=test

print ("The p_values dict is")
print (p_values) ################# dictionary position:pvalues

pvalues_list=[v for v in p_values.values()]

stats_array=statsmodels.stats.multitest.multipletests(pvalues_list, alpha=0.01, method='holm', is_sorted=False, returnsorted=False)[0]
print ("The corrected p_values are")
print (stats_array)

pvals_adj_list=stats_array.tolist()
print (pvals_adj_list) ########### list with corrected pvals

merged_pvalues_pvals_adj=list(tuple(zip(pvalues_list,pvals_adj_list)))
print ("The merged adjusted pvals is")
print (merged_pvalues_pvals_adj) ########### list of tuples (pvalue before correction, pvalue after correction)

new_pvals_dict=dict(zip(p_values.keys(), merged_pvalues_pvals_adj))
new_pvals_dict=dict(sorted(new_pvals_dict.items()))
print ("The new_pvalues dict is")
print (new_pvals_dict) ########### dictionaty position:(pvalue before correction, pvalue after correction)

recurrent_mutations=[]
for k, v in new_pvals_dict.items():
    if v[1]==True: ############# true for hypothesis rejection
    #if v[1]==False:
        recurrent_mutations.append(int(k))
recurrent_mutations=sorted(recurrent_mutations)
print ("The possible drivers are")
print (recurrent_mutations)
print ("with length", len(recurrent_mutations)) ########## list of possible revurrent mutations keys=probs_dict.keys()


keys=probs_dict.keys()
print (keys)
values_list=[v[1] for k,v in probs_dict.items()]
print (values_list)
new_dictionary={}
new_dictionary=dict(zip(keys,values_list))
print (new_dictionary)
sorted_x = sorted(new_dictionary.items(), key=lambda x: x[1])
sorted_dict = dict(sorted_x)
print (sorted_dict)
colors=['b' for i in sorted_dict.keys()]
for i in recurrent_mutations:
    if str(i) in sorted_dict.keys():
        colors[i]='r'
#colors=['b' for i in sorted_dict.keys()]
#for i in sorted_dict.keys():
 #   if int(i) in recurrent_mutations:
  #      print ("yes")
        #colors[sorted_dict.keys()[i]]='r'
        
plt.bar(sorted_dict.keys(), sorted_dict.values(), color=colors)
plt.legend(["Possible drivers"], loc="upper left")
plt.xlabel("Genome positions")
plt.ylabel("Number of fixations in 100 cancerous simulation runs")
plt.title("Statistical significant recurrent mutations")
plt.show()
#####################################################################################################################


pathway="/home/maria/Diplomatiki_parousiasi/machine_learning_cancerous_new/results"

nruns=100

################## features dataframe ############################

dfs_lista=[]
recurrent_dfs_lista=[]

generations_lista=[]
for i in range(1,nruns+1):
    path=pathway+"{}/Frequencies_mutations_CA.txt".format(i)
    df=pd.read_csv(path,sep="\t",header=None)
    generations_lista.append(len(df)-1)
    #if i==2:
     #   break
#print (generations_lista)
max_generation=max(generations_lista)
print ("The max generation is", max_generation) #### max generation of fixation accross 100 runs

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
    new_big_df['Max_freq_time']=new_big_df.idxmax(axis=1)

    new_big_df['Label']=15*[1]+49*[0] ############## 1 for significant positions and 0 for neutral
    recurrent_df=new_big_df[new_big_df.index.isin(recurrent_mutations)]
    dfs_lista.append(new_big_df)
    recurrent_dfs_lista.append(recurrent_df)
    #if i==2:
     #   break
#print (dfs_lista)
big_df=pd.concat(dfs_lista, ignore_index=False) ########## dataframe with 64000 rows and max_generation columns with frequencies as elements
big_recurrent_df=pd.concat(recurrent_dfs_lista, ignore_index=False) ######### the recurrent mutations dataframe of frequencies accross generations.
#print ("The features df is")
#print (big_df)


gens_to_keep=[i for i in range(0,max_generation+1,100)]
gens_to_keep.append(max_generation)
gens_to_keep.append('Max_freq')
gens_to_keep.append('Fixed_freq')
gens_to_keep.append("Max_freq_time")
gens_to_keep.append('Label')

big_df=big_df[gens_to_keep]
big_recurrent_df=big_recurrent_df[gens_to_keep]

print ("The features df is")
print (big_df) ######### the final features dataframe of frequencies

print ("The recurrent features df is")
print (big_recurrent_df) ####### the final recurrent mutations features dataframe of frequencies

##################################################### PCA ###########################################

features = list(big_df.columns.values)
features=features[:-1]
print (features)

x_feat = big_df.loc[:, features].values
#print (x_feat)
y_feat = big_df.loc[:,['Label']].values
#print (y_feat)
x_feat = StandardScaler().fit_transform(x_feat)
#print (x_feat)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x_feat)
#print (principalComponents)

principalDf = pd.DataFrame(data = principalComponents , columns = ['principal component 1 (80%)', 'principal component 2 (0.09%)'])
print (principalDf)
print ("The explained variance is",pca.explained_variance_ratio_)

df_label=big_df[['Label']]
principalDf.reset_index(drop=True, inplace=True)
df_label.reset_index(drop=True, inplace=True)

finalDf = pd.concat([principalDf, df_label], axis = 1)

fig=plt.figure(figsize=(8,8))  
ax=fig.add_subplot(1,1,1)  
ax.set_xlabel('Principal Component 1',fontsize = 15)  
ax.set_ylabel('Principal Component 2',fontsize = 15)  
ax.set_title('2 PCA Components',fontsize=20)  
targets=[0,1] 
colors=['b','r']  
for target,color in zip(targets,colors):    
     indicesToKeep = finalDf['Label'] == target  
     ax.scatter(finalDf.loc[indicesToKeep,'principal component 1 (80%)'],
              finalDf.loc[indicesToKeep,'principal component 2 (0.09%)'],
             c=color,
             s=50)
ax.legend(targets)
plt.xlabel("Principal Component 1 - 80%")
plt.ylabel("Principal Component 2 - 0.09%")
plt.title("PCA analysis") 
ax.grid()
plt.show()



#####################################################################################

##################### Extraction of features and classes/labels ##############

X = big_df.drop(['Label'], axis=1)
X_recurrent=big_recurrent_df.drop(['Label'], axis=1)

target = big_df['Label']
recurrent_target=big_recurrent_df['Label']



###########################################################################


###################### Random Forest of raw data ###############################

X_train,X_test,y_train,y_test=train_test_split(X, target, test_size=0.20, stratify=target)
print('y_train class distribution of raw data')
print(y_train.value_counts(normalize=True))
print('y_test class distribution of raw data')
print(y_test.value_counts(normalize=True))


scaler=StandardScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
rfc=RandomForestClassifier(random_state=34)
rfc.fit(X_train,y_train)

rfc_predict=rfc.predict(X_test)
#rfc_cv_score = cross_val_score(rfc, selected_df, target, cv=10, scoring="roc_auc")
rfc_cv_score = cross_val_score(rfc, X, target, cv=StratifiedKFold(10), scoring="accuracy")


print("=== Confusion Matrix of raw data ===")
print(confusion_matrix(y_test, rfc_predict))
print('\n')
print("=== Classification Report of raw data ===")
print(classification_report(y_test, rfc_predict))
print('\n')
print("=== All AUC Scores of raw data ===")
print(rfc_cv_score)
print('\n')
print('Accuracy of raw data: ',metrics.accuracy_score(y_test, rfc_predict))
print('Recall of raw data: ',metrics.recall_score(y_test, rfc_predict))
print("=== Mean AUC Score of raw data ===")
print("Mean AUC Score - Random Forest: ", rfc_cv_score.mean())

cm = confusion_matrix(y_test, rfc_predict)

plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('Raw Random Forest Confusion Matrix')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.savefig('/home/maria/machine_learning/Confusion_Matrix_raw.png')
plt.show()


#############################################################################

##################### Recursive Feature Elimination for feature selection ###########################


rfc = RandomForestClassifier(random_state=34)
min_features_to_select = 3
rfecv = RFECV(estimator=rfc, step=1, cv=StratifiedKFold(10), scoring='accuracy',min_features_to_select=min_features_to_select)
rfecv.fit(X, target)

print('Optimal number of features: {}'.format(rfecv.n_features_))
plt.figure(figsize=(16, 9))
plt.title('Recursive Feature Elimination with Cross-Validation', fontsize=18, fontweight='bold', pad=20)
plt.xlabel('Number of features selected', fontsize=14, labelpad=20)
plt.ylabel('% Correct Classification', fontsize=14, labelpad=20)
plt.plot(range(1, len(rfecv.grid_scores_) + 1), rfecv.grid_scores_, color='#303F9F', linewidth=3)
         
plt.savefig('/home/maria/machine_learning/Number of best Features.png')
plt.show()
print("The non important features are", np.where(rfecv.support_ == False)[0])
non_important_lista=np.where(rfecv.support_ == False)[0]

selected_df=X.drop(X.columns[non_important_lista], axis = 1)
selected_recurrent_df=X_recurrent.drop(X_recurrent.columns[non_important_lista], axis=1)
print ("The final dataframe of best features is")
print (selected_df)
print ("The final recurrent df of best features is")
print (selected_recurrent_df)

feat_importances = pd.Series(rfecv.estimator_.feature_importances_, index=selected_df.columns)
feat_importances.nlargest(rfecv.n_features_).plot(kind='barh')
plt.ylabel("Features")
plt.xlabel("Contributing percentage in accuracy")
plt.title("Best features")
plt.show()

####################################################################################

############################ Random Forest after feature selection #############################################

X_train,X_test,y_train,y_test=train_test_split(selected_df, target, test_size=0.20, stratify=target)


scaler=StandardScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
rfc=RandomForestClassifier(random_state=34)
rfc.fit(X_train,y_train)

rfc_predict=rfc.predict(X_test)
#rfc_cv_score = cross_val_score(rfc, selected_df, target, cv=10, scoring="roc_auc")
rfc_cv_score = cross_val_score(rfc, selected_df, target, cv=StratifiedKFold(10), scoring="accuracy")


print("=== Confusion Matrix after feature selection ===")
print(confusion_matrix(y_test, rfc_predict))
print('\n')
print("=== Classification Report after feature selection ===")
print(classification_report(y_test, rfc_predict))
print('\n')
print("=== All AUC Scores after feature selection ===")
print(rfc_cv_score)
print('\n')
print('Accuracy after feature selection: ',metrics.accuracy_score(y_test, rfc_predict))
print('Recall after feature selection: ',metrics.recall_score(y_test, rfc_predict))
print("=== Mean AUC Score ===")
print("Mean AUC Score after feature selection - Random Forest: ", rfc_cv_score.mean())


cm = confusion_matrix(y_test, rfc_predict)

plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('Raw Random Forest Confusion Matrix after feature selection')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.savefig('/home/maria/machine_learning/Confusion_Matrix_after_feature_selection.png')
plt.show()

#################################################################

############################ Random Forest parameters improvement ############################################

n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# number of features at every split
max_features = ["auto", "sqrt"]

# max depth
max_depth = [int(x) for x in np.linspace(100, 500, num = 11)]
max_depth.append(None)
# create random grid
random_grid = {
 "n_estimators": n_estimators,
 "max_features": max_features,
 "max_depth": max_depth
 }
# Random search of parameters
rfc_random = RandomizedSearchCV(estimator = rfc, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=34, n_jobs = -1)
# Fit the model
rfc_random.fit(X_train, y_train)
# print results
print ("The best params for random forest are")
print(rfc_random.best_params_)
best_params=rfc_random.best_params_ ############# best parameters for random forest


rfc = RandomForestClassifier(n_estimators=best_params['n_estimators'], max_depth=best_params['max_depth'], max_features=best_params['max_features'], random_state=34) ####### random forest model classifier with the best params
rfc.fit(X_train,y_train)
rfc_predict = rfc.predict(X_test)
rfc_cv_score = cross_val_score(rfc, selected_df, target, cv=StratifiedKFold(10), scoring='accuracy')
print("=== Confusion Matrix after feature selection and params improvement ===")
print(confusion_matrix(y_test, rfc_predict))
print('\n')
print("=== Classification Report after feature selection and params improvement ===")
print(classification_report(y_test, rfc_predict))
print('\n')
print("=== All AUC Scores after feature selection and params improvement ===")
print(rfc_cv_score)
print('\n')
print('Accuracy after parameters improvement: ',metrics.accuracy_score(y_test, rfc_predict))
print('Recall after parameters improvement: ',metrics.recall_score(y_test, rfc_predict))
print("=== Mean AUC Score after parameters improvement ===")
print("Mean AUC Score - Random Forest: ", rfc_cv_score.mean())

cm = confusion_matrix(y_test, rfc_predict)

plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('Random Forest Confusion Matrix after parameters improvement')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.savefig('/home/maria/machine_learning/Confusion_Matrix_after_parameters_improvement.png')
plt.show()

#######################################################################################

################## SMOTE for rebalancing the dataset ################################

sm = SMOTE(random_state=34)
X_sm, y_sm = sm.fit_resample(selected_df, target) ####### SMOTE resampling for features datafeame after features selection
X_sm_recurrent, y_sm_recurrent=sm.fit_resample(selected_recurrent_df, recurrent_target) #### SMOTE resamplinf for recurrent features df

print(f'''Shape of X before SMOTE: {X.shape}
Shape of X after SMOTE: {X_sm.shape}''')

print('\nBalance of positive and negative classes after SMOTE (%):')
print (y_sm.value_counts(normalize=True) * 100)

#print ("geiaaaaa")
#print (pd.concat([X_sm,y_sm],ignore_index=False, axis=1))

########################################################################################

################ Random Forest after SMOTE ###################################

X_train, X_test, y_train, y_test = train_test_split(X_sm, y_sm, test_size=0.20, random_state=34)

scaler=StandardScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)


model = RandomForestClassifier(n_estimators=best_params['n_estimators'], max_depth=best_params['max_depth'], max_features=best_params['max_features'], random_state=34)
model.fit(X_train, y_train)
preds = model.predict(X_test)
model_cv_score = cross_val_score(model, X_sm, y_sm, cv=StratifiedKFold(n_splits=10), scoring='accuracy')
metrics.plot_roc_curve(model, X_test, y_test)
plt.title("Random Forest ROC-AUC curve")
plt.show()  


print("=== Confusion Matrix after SMOTE  ===")
print(confusion_matrix(y_test, preds))
print('\n')
print("=== Classification Report after SMOTE ===")
print(classification_report(y_test, preds))
print('\n')

print(f'Accuracy after SMOTE = {accuracy_score(y_test, preds):.2f}\nRecall after SMOTE = {recall_score(y_test, preds):.2f}\n')
print("=== Mean AUC Score after selection and SMOTE ===")
print("Mean AUC Score - Random Forest: ", model_cv_score.mean())
cm = confusion_matrix(y_test, preds)
plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('Random Forest Confusion Matrix with feature selection and SMOTE')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.savefig('/home/maria/machine_learning/Confusion_Matrix_with_feature_selection_and_SMOTE.png')
plt.show()
#metrics.plot_roc_curve(model, X_test, y_test)
#plt.show()  

model.fit(selected_df, target)
total_preds=model.predict(selected_df)
print ("the type is",type(total_preds))
print (collections.Counter(total_preds))
predictions_df=pd.concat([selected_df,target], axis=1)
predictions_df["Predicted_Label"]=total_preds
print ("The predictions dataframe is")
print (predictions_df)
predictions_df.to_csv("data.txt", sep='\t', index=True, header=True)

predictions_index_array=predictions_df[predictions_df['Predicted_Label']==1].index.values
print ("The indexes of possible drivers are")
print(predictions_index_array)
predictions_dictionary=collections.Counter(predictions_index_array)
print ("The occurences of possible drivers are")

predictions_str_dictionary={str(k):v for k, v in predictions_dictionary.items()}
sorted_x_rf = sorted(predictions_str_dictionary.items(), key=lambda x: x[1], reverse=True)
sorted_dict_rf = dict(sorted_x_rf)
print ("The sorted dict is")
print (sorted_dict_rf)
#print (sorted_dict)
#print (predictions_str_dictionary)
str_list=[str(i) for i in range(15)]
colors_new=['b' for i in sorted_dict_rf.keys()]
for i in predictions_dictionary.keys():
    if str(i) in str_list:
        colors_new[i]='r'
        '''
for i in sorted_dict_rf.keys():
    if str(i) in str_list:
        colors_new[i]='r'
        '''
plt.bar(sorted_dict_rf.keys(), sorted_dict_rf.values(), color=colors_new)
plt.legend(["Class 1"], loc="upper left")
plt.xlabel("Genome Positions")
plt.ylabel("Frequencies of positive class")
plt.title("Classification of class 1 for each genome position")
#plt.bar(*zip(*predictions_dictionary.items()))
plt.show()


##############################################################################

######################### Random Forest of recurrent after SmOTE ##########################################

X_train, X_test, y_train, y_test = train_test_split(X_sm_recurrent, y_sm_recurrent, test_size=0.20, random_state=34)

scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)


model = RandomForestClassifier(n_estimators=best_params['n_estimators'], max_depth=best_params['max_depth'], max_features=best_params['max_features'],random_state=34)
model.fit(X_train, y_train)
preds = model.predict(X_test)
model_cv_score = cross_val_score(model, X_sm, y_sm, cv=StratifiedKFold(10), scoring='accuracy')

print("=== Confusion Matrix of recurrent after SMOTE  ===")
print(confusion_matrix(y_test, preds))
print('\n')
print("=== Classification Report of recurrent after SMOTE ===")
print(classification_report(y_test, preds))
print('\n')


print(f'Accuracy of recurrent after SMOTE = {accuracy_score(y_test, preds):.2f}\nRecall of recurrent after SMOTE = {recall_score(y_test, preds):.2f}\n')
print("=== Mean AUC Score of recurrent ===")
print("Mean AUC Score - Random Forest of recurrent: ", model_cv_score.mean())
cm = confusion_matrix(y_test, preds)
plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('Random Forest Confusion Matrix with SMOTE of recurrent')
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


model.fit(selected_recurrent_df, recurrent_target)
total_preds=model.predict(selected_recurrent_df)
print (collections.Counter(total_preds))
predictions_recurrent_df=pd.concat([selected_recurrent_df, recurrent_target], axis=1)
predictions_recurrent_df["Predicted_Label of recurrent"]=total_preds
print ("The predictions of recurrent dataframe is")
print (predictions_recurrent_df)

predictions_recurrent_index_array=predictions_recurrent_df[predictions_recurrent_df['Predicted_Label of recurrent']==1].index.values
print ("The indexes of possible drivers from recurrent are")
print(predictions_recurrent_index_array)
predictions_recurrent_dictionary=collections.Counter(predictions_recurrent_index_array)
print ("The occurences of possible drivers from recurrent are")
print (predictions_recurrent_dictionary)
############################################################################################

####################################### SVM ########################################

X =big_df.drop('Label', axis=1)
print (X)
y = big_df['Label']
#print (y)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.20, stratify=y)
scaler=StandardScaler()
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)
svclassifier = SVC(kernel='poly',probability=True) ######### polynomial kernel due to unknown distribution of the features data.
svclassifier.fit(X_train, y_train)
y_pred = svclassifier.predict(X_test)
#print (y_pred)
#print (len(y_pred))
#print ("geiaaa",set(y_test) - set(y_pred))
print(confusion_matrix(y_test,y_pred))
print(classification_report(y_test,y_pred))
#print("geiaaaa",metrics.f1_score(y_test, y_pred, labels=np.unique(y_pred)))
print("Accuracy of SVM raw:",metrics.accuracy_score(y_test, y_pred))
print("Precision of SVM raw:",metrics.precision_score(y_test, y_pred))
print("Recall of SVM raw:",metrics.recall_score(y_test, y_pred))
print("F1 of SVM raw:",metrics.f1_score(y_test, y_pred))
#print (svclassifier.predict_proba(X_test))


cm = confusion_matrix(y_test, y_pred)

plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('SVM Poly Kernel Confusion Matrix of Raw')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.savefig('/home/maria/machine_learning/SVM_Poly_confusion_matrix_Raw.png')
plt.show()

################################################################################################

############################## SVM after SMOTE ############################################

X =big_df.drop('Label', axis=1)
y = big_df['Label']
X_sm, y_sm = sm.fit_resample(X,y)
X_train, X_test, y_train, y_test = train_test_split(X_sm, y_sm, test_size = 0.20, random_state=34)
scaler.fit(X_train)
X_train = scaler.transform(X_train)
X_test = scaler.transform(X_test)


model = SVC(kernel='poly')
model.fit(X_train, y_train)
preds = model.predict(X_test)
model_cv_score = cross_val_score(model, X_sm, y_sm, cv=StratifiedKFold(10), scoring='accuracy')

print(f'Accuracy of SVM after SMOTE = {accuracy_score(y_test, preds):.2f}\nRecall of SVM after SMOTE = {recall_score(y_test, preds):.2f}\n')
print ('Precision of SVM after SMOTE', metrics.precision_score(y_test, preds))
print ('F1 of SVM after SMOTE', metrics.f1_score(y_test, preds))
cm = confusion_matrix(y_test, preds)
print (cm)
plt.imshow(cm, interpolation='nearest', cmap=plt.cm.Wistia)
classNames = ['Negative','Positive']
plt.title('Poly SVM Confusion Matrix with SMOTE')
plt.ylabel('True label')
plt.xlabel('Predicted label')
tick_marks = np.arange(len(classNames))
plt.xticks(tick_marks, classNames, rotation=45)
plt.yticks(tick_marks, classNames)
s = [['TN','FP'], ['FN', 'TP']]
 
for i in range(2):
    for j in range(2):
        plt.text(j,i, str(s[i][j])+" = "+str(cm[i][j]))

plt.savefig('/home/maria/machine_learning/Poly_SVM_Confusion_Matrix_after_SMOTE.png')
plt.show()







