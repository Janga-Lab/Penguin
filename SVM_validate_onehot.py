# To set seed random number in order to reproducable results in keras
from numpy.random import seed
seed(4)
#import tensorflow
#tensorflow.random.set_seed(1234)
########################################
import pandas as pd
from pandas import *
import numpy as np
import random
import pickle
import joblib
from sklearn import svm
classifier =svm.SVC(gamma='scale',C=1,probability=True)
import plot_learning_curves as plc
from sklearn.preprocessing import MinMaxScaler #For feature normalization
scaler = MinMaxScaler()

df1 = pd.read_csv("Pseu_Modification_coors_ns_hek.txt",sep=' ',skiprows=(0),header=(0))
df2 = pd.read_csv("xaa.txt",sep='\t',skiprows=(0),header=(0))
df3 = pd.read_csv("Pseu_Modification_coors_hela.txt",sep=' ',skiprows=(0),header=(0))
df4 = pd.read_csv("hela-reads-ref.eventalign.txt",sep='\t',skiprows=(0),header=(0))

print(df2.shape)
print("&&&&&&&&")
print(df1.head())
print("***********************")
print(df2.head())
print("######################")
print(df3.head())
print("######################")
print(df4.head())
#get the 2nd column of 
#column_2 = df1.iloc[:, 1]

#print(df1['position'].iloc[0:5])
print(df1.iloc[0:5, 1])
print("@@@@@@@@@@@@@@@")
#print(df2['position'].iloc[0:5])
#print(df2.iloc[0:5, 1])
print(df2.iloc[0:5, 9])
print("######################")
model_kmer_list=list(df2.iloc[:, 9]) #10 for model-kmer that
model_kmer_list_test=list(df4.iloc[:, 9]) #10 for model-kmer that
print("333333333333333333", type(model_kmer_list))
print("333333333333333333", type(model_kmer_list_test))
print(model_kmer_list[5])
print(model_kmer_list[5][2])
print(model_kmer_list_test[5])
print(model_kmer_list_test[5][2])
U_kmer_list=[]
for i in model_kmer_list:
    #print(i)
    if i[2]=='T':
        U_kmer_list.append(i)

U_kmer_list_test=[]
for i in model_kmer_list_test:
    #print(i)
    if i[2]=='T':
        U_kmer_list_test.append(i)
print("length of U_kmer_list",len(U_kmer_list))
print(U_kmer_list[0:50])

print("length of U_kmer_list_test",len(U_kmer_list_test))
print(U_kmer_list_test[0:50])

df=df2[df2['model_kmer'].isin(U_kmer_list)]
print(df.shape)
print(df.head())

np.savetxt('filtered_df_U_kmer', df,fmt='%s')


df_test=df4[df4['model_kmer'].isin(U_kmer_list_test)]
print(df_test.shape)
print(df_test.head())

np.savetxt('filtered_df_U_kmer_test', df_test,fmt='%s')



#label the data
x=list(set(df1.iloc[:,1]).intersection(set(df.iloc[:,1])))
print("length of intersection list",len(x))

x_test=list(set(df3.iloc[:,1]).intersection(set(df_test.iloc[:,1])))
print("length of intersection list",len(x_test))



df_pseu=df[df['position'].isin(x)]
listofones = [1] * len(df_pseu.index)
# Using DataFrame.insert() to add a column 
df_pseu.insert(13, "label", listofones, True)

df_pseu_test=df_test[df_test['position'].isin(x_test)]
listofones = [1] * len(df_pseu_test.index)
# Using DataFrame.insert() to add a column 
df_pseu_test.insert(13, "label", listofones, True)


df_U=df[~df['position'].isin(x)]
listofzeros=[0]*len(df_U.index)
df_U.insert(13, "label", listofzeros, True)
print(df_pseu.shape)
print(df_pseu.head())
print(df_U.shape)
print(df_U.head())


df_U_test=df_test[~df_test['position'].isin(x_test)]
listofzeros=[0]*len(df_U_test.index)
df_U_test.insert(13, "label", listofzeros, True)
print(df_pseu_test.shape)
print(df_pseu_test.head())
print(df_U_test.shape)
print(df_U_test.head())
#np.savetxt('pseu_samples.txt', df_pseu,fmt='%s')
#np.savetxt('U_samples.txt', df_U,fmt='%s')
##########prepare training datast       
df_U = df_U.sample(n=len(df_pseu), replace=False) #try replace=false

# Create DataFrame from positive and negative examples
dataset = df_U.append(df_pseu, ignore_index=True)
#dataset['label'] = dataset['label'].astype('category')
columns=['event_level_mean','event_stdv','event_length']
#columns=['event_level_mean']
#columns=['event_stdv']
#columns=['event_length']


##################################
#prepare testing datast       
df_U_test = df_U_test.sample(n=len(df_pseu_test), replace=False) #try replace=false

# Create DataFrame from positive and negative examples
dataset_test = df_U_test.append(df_pseu_test, ignore_index=True)
#dataset_test['label'] = dataset_test['label'].astype('category')

#shuffle the test and train datasets
from sklearn.utils import shuffle
dataset = shuffle(dataset)
dataset_test=shuffle(dataset_test)

#combine onehot_encoding of train and test
union_reference_kmer_set=set(dataset.iloc[:, 2]).union(set(dataset_test.iloc[:, 2]))
union=list(union_reference_kmer_set)
print(len(union))
dataset['reference_kmer']=pd.Categorical(dataset['reference_kmer'], categories=list(union))
dataset_test['reference_kmer']=pd.Categorical(dataset_test['reference_kmer'], categories=list(union))


X_train = dataset[columns]
#insert onehot encoding of reference-kmer in train data
Onehot=pd.get_dummies(dataset['reference_kmer'], prefix='reference_kmer')
X_train= pd.concat([X_train,Onehot],axis=1)
#X_train=Onehot
print("#############",X_train.shape)
print(X_train.head())
#scale training data
X_train= scaler.fit_transform(X_train)
y_train = dataset['label'] 
print(",,,,,,,,",X_train.shape)

X_test = dataset_test[columns]

#insert onehot encoding of reference-kmer in test data
Onehot=pd.get_dummies(dataset_test['reference_kmer'], prefix='reference_kmer')
X_test= pd.concat([X_test,Onehot],axis=1)
#X_test= Onehot

print("#############",X_test.shape)
print(X_test.head())
#scale training data
X_test= scaler.fit_transform(X_test)
y_test = dataset_test['label'] 
print(",,,,,,,,",X_test.shape)

###################################

from sklearn.model_selection import train_test_split

#train, test = train_test_split(df, test_size=0.2)   
#X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=0)

#X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.3) for unblanced dataset

#clf = classifier.fit(X_train,y_train)
clf = classifier.fit(X_train,y_train.ravel())

# Evaluate the model: Model Accuracy, how often is the classifier correct
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
from sklearn.metrics import classification_report #for classifier evaluation
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score # for printing AUC
from sklearn.metrics import confusion_matrix

'''
#save the ML model to test on unseen dataset
filename = 'finalized_model.sav'
pickle.dump(classifier, open(filename, 'wb'))
 
# some time later...
 
# load the model from disk
loaded_model = pickle.load(open(filename, 'rb'))
result = loaded_model.score(X_test, y_test)
print(result)
'''
y_pred = classifier.predict(X_test)
y_prob = classifier.predict_proba(X_test)
y_prob = y_prob[:,1]

print("Accuracy:",metrics.accuracy_score(y_test, y_pred)*100)
 
print(classification_report(y_test, y_pred))
auc=roc_auc_score(y_test.round(),y_pred)
auc = float("{0:.3f}".format(auc))
print("AUC=",auc)
#true negatives c00, false negatives C10, true positives C11, and false positives C01 
#tn c00, fpC01, fnC10, tpC11 
print('CF=',confusion_matrix(y_test, y_pred))
l=confusion_matrix(y_test, y_pred)#https://towardsdatascience.com/accuracy-precision-recall-or-f1-331fb37c5cb9
print('TN=',l.item((0, 0)))
print('FP=',l.item((0, 1)))
print('FN=',l.item((1, 0)))
print('TP=',l.item((1, 1)))
#print(type(X_train), type(y_train))


#plot learning curve: works with all classifier and all features except x(padded signal) as it leads to error with SVM 
#References:https://medium.com/@datalesdatales/why-you-should-be-plotting-learning-curves-in-your-next-machine-learning-project-221bae60c53
import matplotlib.pyplot as plt


plc. plot_learning_curves(classifier, X_train, y_train, X_test, y_test)

# Create plot
#plt.title("Learning Curve")
plt.xlabel("Training Set Size"), plt.ylabel("Accuracy Score"), plt.legend(loc="best")
plt.tight_layout()
plt.savefig('SVM_validation_onehot_LC.png',dpi=300)
plt.savefig('SVM_validation_onehot_LC.svg')
plt.close()


#plot ROC curve: https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
from sklearn import metrics
import numpy as np
import matplotlib.pyplot as plt

fpr, tpr, thresholds = metrics.roc_curve(y_test, y_prob)


# Print ROC curve
plt.plot(fpr,tpr)
#plt.title("ROC Curve")
# axis labels
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.savefig('SVM_validation_onehot_ROC.png',dpi=300)
plt.savefig('SVM_validation_onehot_ROC.svg')

plt.close() 


#############################################
#old code to plot learning curve: works only with RandomForest
#Reference: https://www.dataquest.io/blog/learning-curves-machine-learning/
##################
