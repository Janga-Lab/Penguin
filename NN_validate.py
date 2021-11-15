# To set seed random number in order to reproducable results in keras
from numpy.random import seed
seed(4)
import tensorflow
tensorflow.random.set_seed(1234)
########################################
import pandas as pd
from pandas import *
import numpy as np
import random
import pickle
import joblib
import matplotlib.pyplot as plt
import plot_learning_curves as plc
from sklearn.preprocessing import MinMaxScaler #For feature normalization
scaler = MinMaxScaler()


dataset = pd.read_csv('all_hek_data.csv', index_col=0) #all_hek_dat= Hek293 benchmark dataset
dataset_test = pd.read_csv('all_hela_data.csv', index_col=0)    #all_hela_dat= Hela benchmark dataset 
print("===============")
     
#dataset['label'] = dataset['label'].astype('category')
columns=['event_level_mean','event_stdv','event_length']
#columns=['event_level_mean']
#columns=['event_stdv']
#columns=['event_length']

##################################

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
a,b=X_train.shape

X_test = dataset_test[columns]

#insert onehot encoding of reference-kmer in test data
Onehot=pd.get_dummies(dataset_test['reference_kmer'], prefix='reference_kmer')
X_test= pd.concat([X_test,Onehot],axis=1)
#X_test=Onehot
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


###################Build NN model

from keras.models import Sequential
from keras.layers import Dense 
from keras.optimizers import SGD
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
# Evaluate the model: Model Accuracy, how often is the classifier correct
from sklearn import metrics #Import scikit-learn metrics module for accuracy calculation
from sklearn.metrics import classification_report #for classifier evaluation
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score # for printing AUC
from sklearn.metrics import confusion_matrix



model = Sequential()
model.add(Dense(12, input_dim=b, activation='relu'))
model.add(Dense(8, activation='relu'))
model.add(Dense(1, activation='sigmoid'))

# compile the keras model
#model.compile(loss='binary_crossentropy', optimizer=SGD(lr=0.05, momentum=0.99), metrics=['accuracy'])
model.compile(loss='binary_crossentropy', optimizer='Adam', metrics=['accuracy'])

# Fit the model                            #epochs= ﬁxed number of iterations through the dataset called epochs
#model.fit(X_train, y_train, validation_split=0.2, epochs=150, batch_size=16) #batch_size=the number of instances that are evaluated before a weight update
classifier=model.fit(X_train, y_train, validation_split=0.2, epochs=150, batch_size=len(X_train))

# evaluate the keras model on the same dataset, but the dataset can be devided into training and testing, then fit the model on the training and evalaute the model on testing
_, accuracy_train = model.evaluate(X_train, y_train, verbose=1)
_, accuracy_test = model.evaluate(X_test, y_test, verbose=1)
print('Accuracy on training: %.2f' % (accuracy_train*100))
print('Accuracy on testing: %.2f' % (accuracy_test*100))

#needed for plotting learning curve
train_loss = classifier.history['loss']
val_loss   = classifier.history['val_loss']
train_acc  = classifier.history['accuracy']
val_acc    = classifier.history['val_accuracy']
xc         = range(50)


#plot epochs versus accuracy curve.

plt.plot(train_acc, label="Training accuracy")
plt.plot(val_acc, label="Validation accuracy")
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.legend(loc="best")
#plt.show()
plt.savefig('NN_validate_onehot_LC.png',dpi=300)
plt.savefig('NN_validate_onehot_LC.svg')

plt.close()


plt.plot(train_loss, label="Training loss")
plt.plot(val_loss, label="Validation loss")
plt.xlabel('Epochs')
plt.ylabel('loss')
plt.legend(loc="best")
#plt.show()
plt.savefig('NN_onehot_LL.png',dpi=300)
plt.savefig('NN_onehot_LL.svg')
plt.close()

# make class predictions with the model
predictions = model.predict_classes(X_test)
y_pred = model.predict_classes(X_test)
y_prob = model.predict_proba(X_test)
#y_prob = y_prob[:,1]

print(classification_report(y_test, y_pred))
auc=roc_auc_score(y_test.round(),y_pred)
auc = float("{0:.3f}".format(auc))
print("AUC=",auc)

print('CF=',confusion_matrix(y_test, y_pred))
l=confusion_matrix(y_test, y_pred)#https://towardsdatascience.com/accuracy-precision-recall-or-f1-331fb37c5cb9
print('TN=',l.item((0, 0)))
print('FP=',l.item((0, 1)))
print('FN=',l.item((1, 0)))
print('TP=',l.item((1, 1)))

fpr, tpr, thresholds = metrics.roc_curve(y_test, y_prob)


# Print ROC curve
plt.plot(fpr,tpr)
#plt.title("ROC Curve")
# axis labels
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
#plt.show()
plt.savefig('NN_onehot_validate_ROC.png',dpi=300)
plt.savefig('NN_onehot_validate_ROC.svg')

plt.close() 


