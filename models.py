import numpy as np
import csv
import os
import logomaker
from Bio import motifs
from Bio.Seq import Seq
import math
from scipy.stats import binom
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from sklearn.svm import SVR
import parameters as para
import utils
from scipy import optimize


from tensorflow import keras
from keras import layers
from keras.models import Sequential
from keras.layers.core import Dense
from keras.models import load_model


from keras.layers.core import Activation
from keras.layers.convolutional import Conv2D, MaxPooling2D
from keras.layers.core import Flatten
# from tensorflow.keras import optimizers
from keras import optimizers
from keras.callbacks import EarlyStopping,ModelCheckpoint
from scipy.stats import pearsonr

import pandas as pd
import matplotlib.pyplot as plt

def svr(X, y):
    l = int(len(X)*0.8)
    X_train = X[0:l]
    y_train = y[0:l]
    X_test = X[l:len(X)]
    y_test = y[l:len(X)]
  ###############################################################################
    # Fit regression model
    svr_rbf = SVR(kernel='rbf', C=1e3, gamma=0.1).fit(X_train, y_train)
    svr_lin = SVR(kernel='linear', C=1e3).fit(X_train, y_train)
    svr_poly = SVR(kernel='poly', C=1e3, degree=2).fit(X_train, y_train)
    svr_sig = SVR(kernel='sigmoid').fit(X_train, y_train)
    y_sig = svr_sig.predict(X_test)
    y_sig_train = svr_sig.predict(X_train)
    y_rbf = svr_rbf.predict(X_test)
    y_rbf_train = svr_rbf.predict(X_train)
    y_lin = svr_lin.predict(X_test)
    y_poly = svr_poly.predict(X_test)
    ###############################################################################
    # look at the results
    lw = 2
    plt.scatter(y_train, y_rbf_train, color='navy', lw=lw, label='train')
    # plt.plot(x)
    # plt.scatter(y_test, y_rbf, color='c', lw=lw, label='test')
    # plt.scatter(y_test, y_lin, color='c', lw=lw, label='Linear model')
    # plt.scatter(y_test, y_poly, color='cornflowerblue', lw=lw, label='Polynomial model')
    plt.xlabel('data')
    plt.ylabel('target')
    plt.title('Support Vector Regression')
    plt.legend()
    plt.show()
    plt.close('all')
    ###############################################################################
    LR = LinearRegression()
    reg_f = LR.fit(X, y)
    reg_s = LR.score(X,y, sample_weight=None)
    sf = reg_f.score(X, y)
    coef = reg_f.coef_
    inter = reg_f.intercept_
    a = 0


def fully_connect(x,y,x_a,y_a,writer,i,lors,dataform):
    x = np.array(x)
    y = np.log(np.array(y)+0.0001).reshape(-1)
    x_a = np.array(x_a)
    y_a = np.log(np.array(y_a)*max(y)+0.0001).reshape(-1)
    l = int(len(y) * 0.8)
    input_d = len(x[0])
    coxy = []

    for m in range(1):
        ## generate training and testing dataset
        indices = np.random.permutation(x.shape[0])
        training_idx, test_idx = indices[:l], indices[l:]
        print('Building model...')
        ## Fully Connected layer (FC) 200 nodes
        ## parameters not optimized

        model = Sequential()
        model.add(Dense(200, input_dim=input_d, activation='tanh'))
        model.add(Dense(1, activation='relu'))

        print('Compiling model...')
        model.compile(loss='mean_squared_error',
                    optimizer='adam',
                    metrics=['accuracy'])

        print(model.summary())
        model.fit(x[training_idx], y[training_idx], epochs=100, batch_size = 20)
        print('train completed')
        y_test = model.predict(x[test_idx])
        y_test2 = model.predict(x[training_idx])
        y_test3 = model.predict(x_a)
        

        correlation_matrix = np.corrcoef(y[test_idx], y_test.reshape(-1))
        correlation_xy = correlation_matrix[0,1]
        print("testing set correlations ", correlation_xy)
        coxy.append(correlation_xy)

        correlation_matrix = np.corrcoef(y[training_idx], y_test2.reshape(-1))
        correlation_xy = correlation_matrix[0,1]
        print("training set correlations ", correlation_xy)
        coxy.append(correlation_xy)

        correlation_matrix = np.corrcoef(y_a, y_test3.reshape(-1))
        correlation_xy = correlation_matrix[0,1]
        print("anderson set correlations ", correlation_xy)
        coxy.append(correlation_xy)
    print('round ' + str(i) + ' ' + lors+dataform)
    writer.writerow({'model':'CNN', 'data':para.data_kind[i],'lors':lors,'dataform':dataform,'train':coxy[1], 'test':coxy[0],'anderson':coxy[2]})
    return writer

def f_1(x, a, b):
    return a*x + b

def plt_model(y1, y2, title, r):
    a1, b1 = optimize.curve_fit(f_1, y1, y2)[0]
    e = math.ceil(max(max(y1),max(y2)))
    s = int(min(min(y1),min(y2)))

    x = np.linspace(s,e,e-s)
    plt.plot(x,a1*x+b1,'red',label = 'y = '+str(round(a1,2))+'*x + '+str(round(b1,2)))
    plt.scatter(y1, y2)
    plt.title(title+', r = ' + str(r))
    plt.legend()
    plt.savefig(para.FIG_DIR+'model/'+title+'.png')
    plt.show()
    plt.close()


def model_out(x, model, title, y, reshape):
    yp = model.predict(x)
    if reshape:
        yp = yp.reshape(-1)
    # yp = yp[:,0]
    co = pearsonr(yp, y)[0]
    # co = correlation_xy(y,yp)
    print(title+" set correlations ", co)
    plt_model(y, yp, title, co)
    return co

class PREDICT():   #Predict the sequence expression
    
    def __init__(self):  
        self.model_weight = para.TRAIN_DIR+'/weight_CNN.h5'
        self.shuffle_flag = 2
        
    def CNN_model(self,promoter_length):
        model = Sequential()
        model.add(
                Conv2D(100, (6, 1),
                padding='same',
                input_shape=(promoter_length, 1, 4))
                )
        model.add(Activation('relu'))
        model.add(MaxPooling2D(pool_size=(2, 1)))
        model.add(Conv2D(200, (5, 1),padding='same'))
        model.add(Activation('relu'))
        model.add(Conv2D(200, (5, 1),padding='same'))
        model.add(Activation('relu'))
        model.add(MaxPooling2D(pool_size=(2, 1)))
        model.add(Flatten())
        model.add(Dense(1024))
        model.add(Activation('relu'))
        model.add(Dense(1))
        return model

    def CNN_predict(self,seq_onehot):  #Using trained CNN model to predict the expression of promoters.
        model = self.CNN_model(len(seq_onehot[0]))
        # model.load_weights(self.model_weight)
        batch_exp = model.predict(seq_onehot,verbose=0)
        return batch_exp   #return the expression of promoters in this batch

    def random_perm(self,seq,exp,shuffle_flag):  #random perm the sequence and expression data
        indices = np.arange(seq.shape[0])
        np.random.seed(shuffle_flag)
        np.random.shuffle(indices)
        seq = seq[indices]
        exp = exp[indices]
        return seq,exp
        
        
    def SVR_train(self,data,expression):
        data,expression = self.random_perm(data,expression,self.shuffle_flag)
        data = data.reshape(len(data),data.shape[1])
        from sklearn.svm import SVR  
        from sklearn.model_selection import GridSearchCV 
        svr = GridSearchCV(SVR(kernel='rbf', gamma=0.1), cv=5,  
                           param_grid={"C": [1e0, 1e1, 1e2, 1e3],  
                                       "gamma": np.logspace(-2, 2, 5)})  
        svr.fit(data, expression)  
        return svr
    
    def CNN_train(self,data,expression,x_a_cnn,y_a, title):        
        # expression_new = np.zeros((len(expression),))
        # i = 0
        # while i < len(expression):
        #     expression_new[i] = float(expression[i])
        #     i = i + 1

        data,expression = self.random_perm(data,expression,self.shuffle_flag + 1) # Used different shuffle flag from SVR_train
        # Split training/validation and testing set
        r = int(len(expression)*0.8)
        train_feature = data[0:r]
        test_feature = data[r:]
        train_label = expression[0:r]
        test_label = expression[r:]
        # construct CNN model and training 
        model = self.CNN_model(len(data[0]))
        # sgd = optimizers.SGD(lr=0.0005, decay=1e-6, momentum=0.9, nesterov=True)
        sgd = optimizers.gradient_descent_v2.SGD(lr=0.0005, decay=1e-6, momentum=0.9, nesterov=True)

        model.compile(loss='mean_squared_error', optimizer=sgd)
        model.fit(train_feature,train_label, epochs = 1000, batch_size = 128, validation_split=0.1, 
                callbacks=[EarlyStopping(patience=10)],shuffle=True)
        # model.load_weights(self.model_weight)

        co_train = model_out(train_feature, model, title + 'CNN_training', train_label,True)
        co_test = model_out(test_feature, model, title + 'CNN_test', test_label,True)
        co_anderson = model_out(x_a_cnn, model, title + 'CNN_anderson', y_a,True)
        
        # result = model.predict(test_feature, verbose=0)
        # result = result[:,0]
        # cor_test = pearsonr(test_label,result)
        # co_test = utils.correlation_xy(test_label,result)
        # print("testing set correlations ", cor_test)

        # result = model.predict(train_feature, verbose=0)
        # result = result[:,0]
        # cor_train = pearsonr(train_label,result)
        # co_train = utils.correlation_xy(train_label,result)
        # print("training set correlations ", cor_train)

        # result = model.predict(x_a_cnn, verbose=0)
        # result = result[:,0]
        # cor_anderson = pearsonr(y_a,result)
        # co_anderson = utils.correlation_xy(y_a,result)
        # print("testing set correlations ", cor_anderson)

        return co_train, co_test, co_anderson

    def fully_conn(self,x,y,x_a,y_a,title):
        l = int(len(y) * 0.8)
        input_d = len(x[0])

        for m in range(1):
            ## generate training and testing dataset
            indices = np.random.permutation(x.shape[0])
            training_idx, test_idx = indices[:l], indices[l:]
            print('Building model...')
            ## Fully Connected layer (FC) 200 nodes
            ## parameters not optimized

            model = Sequential()
            model.add(Dense(200, input_dim=input_d, activation='tanh'))
            model.add(Dense(1, activation='relu'))

            print('Compiling model...')
            model.compile(loss='mean_squared_error',
                        optimizer='adam',
                        metrics=['accuracy'])

            print(model.summary())
            model.fit(x[training_idx], y[training_idx], epochs=100, batch_size = 20)
            print('train completed')
            co_train = model_out(x[training_idx], model, title + 'NN_training', y[training_idx],True)
            co_test = model_out(x[test_idx], model, title + 'NN_test', y[test_idx],True)
            co_anderson = model_out(x_a, model, title + 'NN_anderson', y_a,True)
            
            # y_test = model.predict(x[test_idx])
            # y_test2 = model.predict(x[training_idx])
            # y_test3 = model.predict(x_a)
            
            # correlation_matrix = np.corrcoef(y[test_idx], y_test.reshape(-1))
            # co_test = correlation_matrix[0,1]
            # print("testing set correlations ", co_test)

            # correlation_matrix = np.corrcoef(y[training_idx], y_test2.reshape(-1))
            # co_train = correlation_matrix[0,1]
            # print("training set correlations ", co_train)

            # y_max = max(y_test3)
            # correlation_matrix = np.corrcoef(y_a, y_test3.reshape(-1))
            # co_anderson = correlation_matrix[0,1]
            # print("anderson set correlations ", co_anderson)
        
        return co_train, co_test, co_anderson

    def predict_svr(self,x,y,x_a,y_a, title):
        l = int(len(y) * 0.8)
        svr = self.SVR_train(x[0:l],y[0:l])
        # seq_svr = x.reshape(len(x),x.shape[1] * 4)
        co_train = model_out(x[:l], svr, title + 'SVR_training', y[:l],False)
        co_test = model_out(x[l:], svr, title + 'SVR_test', y[l:],False)
        co_anderson = model_out(x_a, svr, title + 'SVR_anderson', y_a,False)
        # yp_train = svr.predict(x[:l])
        # yp_test = svr.predict(x[l:])
        # yp_anderson = svr.predict(x_a)
        # co_train = utils.correlation_xy(y[0:l],yp_train)
        # co_test = utils.correlation_xy(y[l:],yp_test)
        # co_anderson = utils.correlation_xy(y_a,yp_anderson)
        
        # print("testing set correlations ", co_test)
        # print("training set correlations ", co_train)
        # print("anderson set correlations ", co_anderson)
        # plt_model(y1, y2, title, r)
        # plt_model(y1, y2, title, r)
        # plt_model(y1, y2, title, r)
        
        return co_train, co_test, co_anderson

    def predict(self,x,y,x_a,y_a,writer,i,lors,dataform):
        print('round ' + str(i) + ' ' + lors+dataform)
        
        x = np.array(x)
        y = np.log(np.array(y)+0.0001).reshape(-1)
        x_a = np.array(x_a)
        y_a = np.log(np.array(y_a)*max(y)+0.0001).reshape(-1)

        # Predict by SVR
        co_train, co_test, co_anderson = self.predict_svr(x,y,x_a,y_a,para.data_kind[i]+dataform)
        writer.writerow({'model':'SVR', 'data':para.data_kind[i],'lors':lors,'dataform':dataform,
                        'train':co_train, 'test':co_test,'anderson':co_anderson})
        
        co_train, co_test, co_anderson = self.fully_conn(x,y,x_a,y_a,para.data_kind[i]+dataform)
        writer.writerow({'model':'nn', 'data':para.data_kind[i],'lors':lors,'dataform':dataform,
                        'train':co_train, 'test':co_test,'anderson':co_anderson})
        # Predict by CNN
        n = len(x)
        line = len(x[0])
        x_cnn = x.reshape(n,int(line/4),1,4)
        x_a_cnn = x_a.reshape(len(x_a),int(line/4),1,4)
        co_train, co_test, co_anderson = self.CNN_train(x_cnn,y,x_a_cnn,y_a,para.data_kind[i]+dataform)
        writer.writerow({'model':'CNN', 'data':para.data_kind[i],'lors':lors,'dataform':dataform,
                        'train':co_train, 'test':co_test,'anderson':co_anderson})

        return writer