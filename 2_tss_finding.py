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
# from gtfparse import read_gtf
from pandas.core.construction import is_empty_data
import scipy
from scipy.stats import binom
import xlwt
import time
import info
import seqlogo

inputf = info.inp
inputs = info.srand

def potential_promoter_cdf(ddict, gdict):
    fp= open(info.PROMOTER_DIR+inputf+'promoter.txt','w')
    data=np.load(info.GENE_DIR+inputf+'gene_protein.npy',allow_pickle=True)
    promoter_dict= {'+':{},'-':{}} #tss:iseq
    tss_dict = {'+':{},'-':{}} #tss:tss
    gene=data.tolist()
    for srand, gitem in gene.items():
        depth = ddict[srand]
        genome = gdict[srand]
        allreads = sum(depth)

        for item in gitem:
            start = item.start
            ratiom = 0
            for i in range(start-100, start):
                high = np.mean(depth[i:i+10])
                low = np.mean(depth[i-10:i])
                ratio = high/(low+10)
                if (info.if_start(i, depth) or (low == 0 and depth[i]> 10) )and ratio>ratiom:
                    item.tss = i
                    ratiom = ratio
            tss = item.tss
            if tss >= start - 100:
                item.promoter = genome[tss-100:tss]
                high = high = np.mean(depth[tss:tss+10])
                low = np.mean(depth[tss-10:tss])
                item.strength = (high - low)/info.copy_num(tss)
                item.avedepth = np.mean(depth[item.start:item.end+1])
                item.fpkm = pow(10,9) / allreads * item.avedepth
                tss_dict[srand][str(tss)] = tss
                promoter_dict[srand][str(tss)] = item
                fp.write(item.promoter+'\n')
    fp.close()
    np.save(info.TSS_DIR+inputf+'tss.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+'gene_after_tss.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+'promoter_dict.npy',promoter_dict)

def potential_promoter_cdf_single(ddict, gdict,srand):
    fp= open(info.PROMOTER_DIR+inputf+srand+'promoter.txt','w')
    data=np.load(info.GENE_DIR+inputf+'gene_protein.npy',allow_pickle=True)
    gene=data.tolist()
    depth = ddict[srand]
    genome = gdict[srand]
    allreads = sum(depth)

    promoter_dict= {} #tss:iseq
    tss_dict = {} #tss:tss
    
    for item in gene[srand]:
        start = item.start
        ratiom = 0
        for i in range(start-100, start):
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            ratio = high/(low+10)
            if (info.if_start(i, depth) or (low == 0 and depth[i]> 10) )and ratio>ratiom:
                item.tss = i
                ratiom = ratio
        tss = item.tss
        if tss >= start - 100:
            item.promoter = genome[tss-100:tss]
            high = high = np.mean(depth[tss:tss+10])
            low = np.mean(depth[tss-10:tss])
            item.strength = (high - low)/info.copy_num(tss)
            item.avedepth = np.mean(depth[item.start:item.end+1])
            item.fpkm = pow(10,9) / allreads * item.avedepth
            tss_dict[str(tss)] = tss
            promoter_dict[str(tss)] = item
            fp.write(item.promoter+'\n')
    fp.close()
    np.save(info.TSS_DIR+inputf+srand+'tss.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+srand+'gene_after_tss.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+srand+'promoter_dict.npy',promoter_dict)

def potential_promoters(ddict, gdict,srand):
    fp= open(info.PROMOTER_DIR+inputf+srand+'promoter.txt','w')
    data=np.load(info.GENE_DIR+inputf+'gene_protein.npy',allow_pickle=True)
    gene=data.tolist()
    depth = ddict[srand]
    genome = gdict[srand]
    allreads = sum(depth)

    promoter_dict= {} #tss:iseq
    tss_dict = {} #tss:tss
    
    for item in gene[srand]:
        start = item.start
        ratiom = 0
        for i in range(start-100, start):
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            ratio = high/(low+10)
            if (binom.cdf(low, low + high, 0.5)<0.01 and high/(low+0.01)>2) and low < 2 and ratio>ratiom:
            # if (info.if_start(i, depth) or (low == 0 and depth[i]> 1) )and ratio>ratiom:
                item.tss = i
                ratiom = ratio
        tss = item.tss
        # start = item.start
        # # info.plot_test(depth,start-100,start)
        # ratiom = 0
        # for i in range(start-100, start):
        #     high = np.mean(depth[i:i+10])
        #     low = np.mean(depth[i-10:i])
        #     ratio = high/(low+10)
        #     if info.if_start(i, depth) and ratio>ratiom:
        #         item.tss = i
        #         ratiom = ratio
        #     # high = np.mean(depth[i:i+10])
        #     # low = max(1,np.mean(depth[i-10:i]))
        #     # ratio = high/low
        #     # if info.if_start(i, depth) and ratio>ratiomin:
        #     #     item.tss = i
        #     #     ratiomin = ratio
        # tss = item.tss
        if tss >= start - 100:
            item.promoter = genome[tss-100:tss]
            high = high = np.mean(depth[tss:tss+10])
            low = np.mean(depth[tss-10:tss])
            item.strength = (high - low)/info.copy_num(tss)
            item.avedepth = np.mean(depth[item.start:item.end+1])
            item.fpkm = pow(10,9) / allreads * item.avedepth
            tss_dict[str(tss)] = tss
            promoter_dict[str(tss)] = item
            fp.write(item.promoter+'\n')

    fp.close()
    np.save(info.TSS_DIR+inputf+srand+'tss.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+srand+'gene_after_tss.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+srand+'promoter_dict.npy',promoter_dict)

def potential_promoters_shaobin(ddict, gdict,srand):
    fp= open(info.PROMOTER_DIR+inputf+srand+'promoter_s.txt','w')
    data=np.load(info.GENE_DIR+inputf+'gene_protein.npy',allow_pickle=True)
    gene=data.tolist()
    depth = ddict[srand]
    genome = gdict[srand]
    allreads = sum(depth)

    promoter_dict= {} #tss:iseq
    tss_dict = {} #tss:tss
    
    for item in gene[srand]:
        start = item.start
        ratiom = 0
        for i in range(start-100, start):
    #         np.mean(Rwig[i:i+10,1])- np.mean(Rwig[i-10:i,1]) > 20*copy_number and \
    # np.mean(Rwig[i:i+10,1])/(np.mean(Rwig[i-10:i,1])+10) > 3 and \
    # np.mean(Rwig[i:i+50,1])>20*copy_number
            high = np.mean(depth[i:i+10])
            high50 = np.mean(depth[i:i+50])
            low = np.mean(depth[i-10:i])
            ratio = high/(low+10)
            copy = info.copy_num(i)
            if high - low > 20*copy and ratio>3 and high50>20*copy and ratio>ratiom:
                item.tss = i
                ratiom = ratio
        tss = item.tss
        if tss >= start - 100:
            item.promoter = genome[tss-100:tss]
            high = high = np.mean(depth[tss:tss+10])
            low = np.mean(depth[tss-10:tss])
            item.strength = (high - low)/info.copy_num(tss)
            item.avedepth = np.mean(depth[item.start:item.end+1])
            item.fpkm = pow(10,9) / allreads * item.avedepth
            tss_dict[str(tss)] = tss
            promoter_dict[str(tss)] = item
            fp.write(item.promoter+'\n')

    fp.close()
    np.save(info.TSS_DIR+inputf+srand+'tss_s.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+srand+'gene_after_tss_s.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+srand+'promoter_dict_s.npy',promoter_dict)

def potential_promoters_5(ddict, gdict,srand):
    fp= open(info.PROMOTER_DIR+inputf+srand+'promoter_5.txt','w')
    data=np.load(info.GENE_DIR+inputf+'gene_protein.npy',allow_pickle=True)
    gene=data.tolist()
    depth = ddict[srand]
    genome = gdict[srand]
    allreads = sum(depth)

    promoter_dict= {} #tss:iseq
    tss_dict = {} #tss:tss
    
    for item in gene[srand]:
        start = item.start
        ratiom = 0
        for i in range(start-100, start):
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            ratio = high/(low+10)
            if ((low == 0 and depth[i]> 1) or high/low >= 5)and ratio>ratiom:
                item.tss = i
                ratiom = ratio
        tss = item.tss
        # start = item.start
        # # info.plot_test(depth,start-100,start)
        # ratiom = 0
        # for i in range(start-100, start):
        #     high = np.mean(depth[i:i+10])
        #     low = np.mean(depth[i-10:i])
        #     ratio = high/(low+10)
        #     if info.if_start(i, depth) and ratio>ratiom:
        #         item.tss = i
        #         ratiom = ratio
        #     # high = np.mean(depth[i:i+10])
        #     # low = max(1,np.mean(depth[i-10:i]))
        #     # ratio = high/low
        #     # if info.if_start(i, depth) and ratio>ratiomin:
        #     #     item.tss = i
        #     #     ratiomin = ratio
        # tss = item.tss
        if tss >= start - 100:
            item.promoter = genome[tss-100:tss]
            high = high = np.mean(depth[tss:tss+10])
            low = np.mean(depth[tss-10:tss])
            item.strength = (high - low)/info.copy_num(tss)
            item.avedepth = np.mean(depth[item.start:item.end+1])
            item.fpkm = pow(10,9) / allreads * item.avedepth
            tss_dict[str(tss)] = tss
            promoter_dict[str(tss)] = item
            fp.write(item.promoter+'\n')

    fp.close()
    np.save(info.TSS_DIR+inputf+srand+'tss_5.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+srand+'gene_after_tss_5.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+srand+'promoter_dict_5.npy',promoter_dict)

def potential_promoters_2(ddict, gdict,srand):
    fp= open(info.PROMOTER_DIR+inputf+srand+'promoter_2.txt','w')
    data=np.load(info.GENE_DIR+inputf+'gene_protein.npy',allow_pickle=True)
    gene=data.tolist()
    depth = ddict[srand]
    genome = gdict[srand]
    allreads = sum(depth)

    promoter_dict= {} #tss:iseq
    tss_dict = {} #tss:tss
    
    for item in gene[srand]:
        start = item.start
        ratiom = 0
        for i in range(start-100, start):
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            ratio = high/(low+10)
            if info.if_start(i, depth) and ratio>ratiom:
                item.tss = i
                ratiom = ratio
        tss = item.tss
        # start = item.start
        # # info.plot_test(depth,start-100,start)
        # ratiom = 0
        # for i in range(start-100, start):
        #     high = np.mean(depth[i:i+10])
        #     low = np.mean(depth[i-10:i])
        #     ratio = high/(low+10)
        #     if info.if_start(i, depth) and ratio>ratiom:
        #         item.tss = i
        #         ratiom = ratio
        #     # high = np.mean(depth[i:i+10])
        #     # low = max(1,np.mean(depth[i-10:i]))
        #     # ratio = high/low
        #     # if info.if_start(i, depth) and ratio>ratiomin:
        #     #     item.tss = i
        #     #     ratiomin = ratio
        # tss = item.tss
        if tss >= start - 100:
            item.promoter = genome[tss-100:tss]
            high = high = np.mean(depth[tss:tss+10])
            low = np.mean(depth[tss-10:tss])
            item.strength = (high - low)/info.copy_num(tss)
            item.avedepth = np.mean(depth[item.start:item.end+1])
            item.fpkm = pow(10,9) / allreads * item.avedepth
            tss_dict[str(tss)] = tss
            promoter_dict[str(tss)] = item
            fp.write(item.promoter+'\n')

    fp.close()
    np.save(info.TSS_DIR+inputf+srand+'tss_2.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+srand+'gene_after_tss_2.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+srand+'promoter_dict_2.npy',promoter_dict)


def open_gene(depth, genome,srand):
    fp= open(info.PROMOTER_DIR+inputf+srand+'promoter.txt','w')
    data=np.load(info.GENE_DIR+info.inp[0:4]+srand+'gene_protein.npy',allow_pickle=True)
    gene=data.tolist()
    allreads = sum(depth)
    
    promoter_dict= {}
    tss_dict = {}

    for item in gene:
        start = item.start
        ratiom = 0
        for i in range(start-100, start):
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            ratio = high/(low+10)
            if info.if_start(i, depth) and ratio>ratiom:
                item.tss = i
                ratiom = ratio
        tss = item.tss
        if tss >= start - 100:
            item.promoter = genome[tss-100:tss]
            high = 0
            for i in range(0,10):
                if np.mean(depth[tss+i:tss+10+i]) > high:
                    high = np.mean(depth[tss+i:tss+10+i])
            item.strength = (high - np.mean(depth[tss-10:tss]))/info.copy_num(tss)
            item.avedepth = np.mean(depth[item.start:item.end+1])
            item.fpkm = pow(10,9) / allreads * item.avedepth
            tss_dict[str(tss)] = tss
            promoter_dict[str(tss)] = item
            fp.write(item.promoter+'\n')
            

    now = time.strftime("%Y-%m-%d-%H_%M_%S",time.localtime(time.time()))
    fp.close()
    np.save(info.TSS_DIR+inputf+srand+'tss.npy',tss_dict)
    np.save(info.GENE_DIR+inputf+srand+'gene_after_tss.npy',gene)
    np.save(info.PROMOTER_DIR+inputf+srand+'promoter_dict.npy',promoter_dict)

if __name__ == '__main__':
    gdict = info.open_gdict_npy()
    g_len = len(gdict['+'])
    ddict = info.open_npy(info.indid+'ddict.npy')
    potential_promoter_cdf(ddict,gdict)
    info.gen_fasta(inputf)
    info.gen_logo(inputf)
    
    