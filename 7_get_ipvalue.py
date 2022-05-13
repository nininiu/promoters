# from curses.ascii import SP
from math import log2
from os import close, startfile
from turtle import st
from typing import Iterator
import numpy as np
import argparse
import matplotlib.pyplot as plt
from numpy.lib.function_base import average
from numpy.lib.utils import deprecate, deprecate_with_doc
import pandas as pd
import pdb
from pandas.core.construction import is_empty_data
import scipy
from scipy.stats import binom
import xlwt
import time
import info
import math
from Bio import motifs
from Bio.Seq import Seq
from tensorflow import keras
import os

import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


# 967:73-87
inputf = info.inp
inputs = info.srand
inputf4 = info.inp4

tf_keras = False

def read_result(fin):
    pv = {}
    s = 0
    f = open(fin,'r')
    for item in f:
        if 'Start' in item:
            s = 1
            continue
        if s > 1:
            if '------' in item:
                break
            item = item.split()
            keylist = item[0].split('.')
            key = keylist[0]
            srand = keylist[1]
            if srand not in pv:
                pv[srand] = {}
            # if key == '1598476':
            #     a = 0
            pv[srand][key] = info.get_pvalue(item[2])
        if s > 0:
            s+=1
            continue
    f.close()
    return pv

def at_num(seq):
    res = 0
    l = len(seq)
    i = 0
    res = 0
    while i < l:
        cnt = 0
        if seq[i] not in 'AT':
            i+=1
            continue
        while i < l and seq[i] in 'AT':
            cnt+=1
            i+=1
        if cnt > res:
            res = cnt
    return res

def split_sequence(sequence, sw_width, n_features):
    # '''
    # 这个简单的示例，通过for循环实现有重叠截取数据，滑动步长为1，滑动窗口宽度为sw_width。
    # 以后的文章，会介绍使用yield方法来实现特定滑动步长的滑动窗口的实例。
    # '''
    X, y = [], []
    
    for i in range(len(sequence)):
        # 获取单个样本中最后一个元素的索引，因为python切片前闭后开，索引从0开始，所以不需要-1
        end_element_index = i + sw_width
        # 如果样本最后一个元素的索引超过了序列索引的最大长度，说明不满足样本元素个数，则这个样本丢弃
        if end_element_index > len(sequence) - 1:
            break
        # 通过切片实现步长为1的滑动窗口截取数据组成样本的效果
        seq_x, seq_y = sequence[i:end_element_index], sequence[end_element_index]
        
        X.append(seq_x)
        y.append(seq_y)
        
        process_X, process_y = np.array(X), np.array(y)
        process_X = process_X.reshape((process_X.shape[0], process_X.shape[1], n_features))
    
    print('split_sequence:\nX:\n{}\ny:\n{}\n'.format(np.array(X), np.array(y)))
    print('X_shape:{},y_shape:{}\n'.format(np.array(X).shape, np.array(y).shape))
    print('train_X:\n{}\ntrain_y:\n{}\n'.format(process_X, process_y))
    print('train_X.shape:{},trian_y.shape:{}\n'.format(process_X.shape, process_y.shape))
    return process_X, process_y




if __name__ == '__main__':
    gdict = info.open_gdict_npy()
    ddict = info.open_npy(info.indid+'ddict.npy')
    promoter_dict = info.open_npy(info.PROMOTER_DIR+inputf+'promoter_dict.npy')
    spacers = info.open_npy_list(info.MOTIF_DIR+inputf+'spacers.npy')
    m10ilist = info.open_npy_list(info.MOTIF_DIR+inputf+'m10ilist.npy')
    mextilist = info.open_npy_list(info.MOTIF_DIR+inputf+'mextilist.npy')
    m35ilist = info.open_npy_list(info.MOTIF_DIR+inputf+'m35ilist.npy')

    pssm10 = info.get_pssm(m10ilist)
    pssmext = info.get_pssm(mextilist)
    pssm35 = info.get_pssm(m35ilist)

    cntspacers = pd.value_counts(spacers)
    print(cntspacers)
    lspacers = info.log_ods(cntspacers)

    pv10 = read_result(info.RESULT_DIR+inputf+'10i')
    pvext = read_result(info.RESULT_DIR+inputf+'exti')
    pv35 = read_result(info.RESULT_DIR+inputf+'35i')

    xpvalue = []
    xpssm = []
    y = []

    fp = open(info.PROMOTER_DIR+'promoters.txt','w')
    p_npy = np.zeros((len(promoter_dict['+'])+len(promoter_dict['-']), 80*4))
    dic = {'A':0, 'T':1, 'C':2, 'G':3}
    pi = 0
    for srand, plist in promoter_dict.items():
        for key, iseq in plist.items():
            tss = iseq.correcttss
            promoter = gdict[srand][tss-60:tss+20]
            fp.write(promoter+'\n')
            for index in range(len(promoter)):
                p_npy[pi][4*index+dic[promoter[index]]] = 1
            pi+=1
    fp.close()
    np.save(info.PROMOTER_DIR+'promoters.npy',p_npy)


    for srand, plist in promoter_dict.items():
        for key, iseq in plist.items():
            pvalueitem = []
            pssmitem = []
            
            promoter_dict[srand][key].m35pvalue = pv35[srand][key]
            promoter_dict[srand][key].m10pvalue = pv10[srand][key]
            promoter_dict[srand][key].mext10pvalue = pvext[srand][key]
            promoter_dict[srand][key].mspacerscore = lspacers[iseq.spacer]
            promoter_dict[srand][key].m35p = promoter_dict[srand][key].m10p-promoter_dict[srand][key].spacer
            promoter_dict[srand][key].mupp = promoter_dict[srand][key].m35p-7
            promoter_dict[srand][key].mupseq = gdict[srand][promoter_dict[srand][key].mupp-20:promoter_dict[srand][key].mupp]
            promoter_dict[srand][key].mupscore = at_num(iseq.mupseq)
            cpssm10 = pssm10.calculate(Seq(iseq.m10seq))
            cpssmext = pssmext.calculate(Seq(iseq.mext10seq))
            cpssm35 = pssm35.calculate(Seq(iseq.m35seq))
            pvalueitem = [-math.log2(iseq.m10pvalue), -math.log2(iseq.mext10pvalue),-math.log2(iseq.m35pvalue), iseq.mspacerscore, iseq.mupscore]
            pssmitem = [cpssm10,cpssmext,cpssm35,iseq.mspacerscore,iseq.mupscore]
            xpvalue.append(pvalueitem)
            xpssm.append(pssmitem)
            # y.append(math.log2(iseq.strength))
            y.append([np.log10(iseq.strength)])
    
    # y = np.log10(y)
    np.save(info.MOTIF_DIR+'xpssm.npy',xpssm)
    np.save(info.MOTIF_DIR+'xpvalue.npy',xpvalue)
    np.save(info.MOTIF_DIR+'y.npy',y)


    if tf_keras:
        x = np.array(xpssm)
        train_x, train_y = x[:1050], y[:1050]
        test_x, test_y = x[1050:], y[1050:]
        model = keras.models.Sequential([
            keras.layers.Dense(10, activation=keras.activations.relu),
            keras.layers.Dense(1),
        ])
        model.compile(
            optimizer=keras.optimizers.SGD(0.01),
            loss=keras.losses.MeanSquaredError(),
            metrics=[keras.metrics.MeanSquaredError()],
        )

        model.fit(train_x, train_y, batch_size=32, epochs=3, validation_split=0.2, shuffle=True)
        model.evaluate(test_x, test_y, verbose=1)
        print(model.summary())
        model.fit(train_x, train_y, epochs=15, batch_size = 20)
        print('train completed')
        y_test = model.predict(test_x)
        y_test2 = model.predict(train_x)

        correlation_matrix = np.corrcoef(test_y, y_test.reshape(-1))
        correlation_xy = correlation_matrix[0,1]
        print("testing set correlations ", correlation_xy)

        correlation_matrix = np.corrcoef(train_y, y_test2.reshape(-1))
        correlation_xy = correlation_matrix[0,1]
        print("training set correlations ", correlation_xy)

    LR = LinearRegression()
    reg_f = LR.fit(xpvalue, y)
    reg_s = LR.score(xpvalue,y, sample_weight=None)
    sf = reg_f.score(xpvalue, y)
    coef = reg_f.coef_
    inter = reg_f.intercept_
    r2 = r2_score(y,LR.predict(xpvalue))
    print(r2)
    print(coef)

    LRpssm = LinearRegression()
    reg_fpssm = LRpssm.fit(xpssm, y)
    reg_spssm = LRpssm.score(xpssm,y, sample_weight=None)
    sfpssm = reg_fpssm.score(xpssm, y)
    coefpssm = reg_fpssm.coef_
    interpssm = reg_fpssm.intercept_
    r2pssm = r2_score(y,LR.predict(xpssm))    
    print(r2pssm)
    print(coefpssm)
    a=0
