import numpy as np
import csv
import os
import logomaker
from Bio import motifs
from Bio.Seq import Seq
import math
from scipy.stats import binom
import matplotlib.pyplot as plt

import parameters as para

def open_gdict_npy():
    return open_npy(para.GTF_DIR+'gdict.npy')

def open_npy(input_file):
    return np.load(input_file,allow_pickle=True).item()

def open_npy_list(input_file):
    data=np.load(input_file,allow_pickle=True)
    data = data.tolist()
    return data

def parse_fpkm(id):
    fpkm_dict = {}
    tpm_dict = {}
    cov_dict = {}
    gtf_file = para.STIE_DIR + id + '.gtf'
    if not os.path.isfile(gtf_file):
        return fpkm_dict, tpm_dict, cov_dict
    data = np.loadtxt(gtf_file, dtype=np.str_, delimiter='\t').tolist()
    for item in data:
        if len(item) > 8 and item[2] == 'transcript' and 'reference_id' in item[8]:
            att_list = item[8].split()
            i = 1
            while (i < len(att_list)):
                if 'gene_id' or 'gene' in att_list[i]:
                    break
                i+=1
            gid = att_list[i+1][1:-2].split('-')[-1]
            if gid not in cov_dict:
                cov_dict[gid] = 0 
                fpkm_dict[gid] = 0
                tpm_dict[gid] = 0
            for atti in range(0, len(att_list)):
                if att_list[atti] == 'cov':
                    cov_dict[gid] += float(att_list[atti + 1][1:-2])
                elif att_list[atti] == 'FPKM':
                    fpkm_dict[gid] += float(att_list[atti + 1][1:-2])
                elif att_list[atti] == 'TPM':
                    tpm_dict[gid] += float(att_list[atti + 1][1:-2])
    return fpkm_dict, tpm_dict, cov_dict

def csv_writer(file, fieldnames):
    cf = open(file, 'w', newline='')
    writer = csv.DictWriter(cf, fieldnames=fieldnames)
    writer.writeheader()
    return writer

def gen_fasta_psi(s, e, id, kind, promoters):
    fa = open(para.FASTA_DIR + id + kind + 'raw.fasta','w+')
    print('generating fasta file:' + id + kind)
    motifs = []
    for gid, pitemlist in promoters.items():
        for i in range(len(pitemlist)):
            pitem = pitemlist[i]
            srand = pitem.srand
            if len(pitem.promoter[s:e]) < e - s:
                continue
            fa.write('>'+gid+'.'+str(i)+'\n')
            fa.write(pitem.promoter[s:e]+'\n')
            motifs.append(pitem.promoter[s:e])
    fa.close()
    plot_logo(motifs, id + kind + 'raw')

def gen_fasta(s, e, id, kind, promoters):
    fa = open(para.FASTA_DIR + id + kind + 'raw.fasta','w+')
    print('generating fasta file:' + id + kind)
    motifs = []
    i = 0
    for gid, pallitem in promoters.items():
        if pallitem['p'][0].tss > 0:
            i += 1
            pitem = pallitem['p'][0]
            srand = pitem.srand
            if len(pitem.promoter[s:e]) < e - s:
                continue
            fa.write('>'+gid+'.'+srand+'\n')
            fa.write(pitem.promoter[s:e]+'\n')
            motifs.append(pitem.promoter[s:e])
    fa.close()
    plot_logo(motifs, id + kind + 'raw')

def plot_logo(promoters, logo_file):
    print('plot '+logo_file+' logo:'+str(len(promoters))+' promoters/motifs')
    if len(promoters) == 0:
        return
    maxLength = max(map(lambda X: len(X), promoters))
    promoters = list(filter(lambda X: len(X) == maxLength, promoters))
    paraMatrix = logomaker.alignment_to_matrix(promoters, to_type = 'information')
    plotLogo = logomaker.Logo(paraMatrix)
    plotLogo.fig.savefig(para.LOGO_DIR+logo_file+'.png', format = 'png')

def get_pssm(seq_list):
    if len(seq_list) == 0:
        return []
    i = []
    maxLength = max(map(lambda X: len(X), seq_list))
    seq_list = list(filter(lambda X: len(X) == maxLength, seq_list))
    for item in seq_list:
        i.append(Seq(item))
    m = motifs.create(i)
    print(m.counts)
    pwm = m.counts.normalize(pseudocounts=0.5)
    print(pwm)
    
    pssm = pwm.log_odds()
    for key,value in pssm.items():
        for i in range(0, len(value)):
            if value[i] <-2.37:
                pssm[key][i] = -2.37
    print(pssm)
    # input()
    return pssm


def compare(i1,i2):
    l1 = get_pwm(i1)
    l2 = get_pwm(i2)
    for key, value in l1.items():
        for i in range(0,len(value)):
            if abs(l1[key][i] - l2[key][i]) > 0.01:
                return False
    return True
    
def get_pwm(seq_list):
    i = []
    maxLength = max(map(lambda X: len(X), seq_list))
    seq_list = list(filter(lambda X: len(X) == maxLength, seq_list))
    for item in seq_list:
        if item == '':
            continue
        i.append(Seq(item))
    m = motifs.create(i)
    pwm = m.counts.normalize(pseudocounts=0.5)
    return pwm

def get_pwm_score(pwm, seq):
    score = 0
    for i in range(0, len(seq)):
        score += pwm[seq[i]][i]
    return score

def log_ods(list):
    res = {}
    s = sum(list)
    l = float(len(list))
    for key,value in list.items():
        p = value/s
        res[key] = -p*math.log2(p)
    return res

def if_start(index, depth):
    # original version
    # r = 1
    # readsi=depth[index]
    # readsi_1= depth[index-1*r]
    # readsi_2= depth[index-2*r]
    # readsi_10 = np.mean(depth[max(index-10,0):index])
    # cdf2 = binom.cdf(readsi_2, readsi_2 + readsi, 0.5)
    # cdf1 = binom.cdf(readsi_1, readsi_1 + readsi, 0.5)
    # return ((cdf1 < 0.01 and readsi / (readsi_1 + 0.0001) > 1.5) 
    #     and (cdf2 < 0.01 and readsi / (readsi_2 + 0.0001) > 1.5)) or (readsi_10 == 0 and readsi >= 5)

    readsi=np.mean(depth[index:min(len(depth), index+5)])
    readsi_1= np.mean(depth[max(index-5,0):index])
    cdf = binom.cdf(readsi_1, readsi_1 + readsi, 0.5)
    return (cdf < 0.01 and readsi / (readsi_1 + 0.0001) > 1.1)  or (readsi_1 == 0 and readsi >= 5)

def copy_num(i):
    copy_number = 1
    if i < para.origin:
        copy_number = 1+3*abs(i-para.terminus)/(para.origin-para.terminus)
    else:
        copy_number = 4-3*(i-para.origin)/(para.origin-para.terminus)
    return copy_number

def get_pvalue(str):
    plist = str.split('e')
    num = float(plist[0])
    e = int(plist[1])
    pvalue = pow(num,e)
    return pvalue

def str2list(str):
    str = str[1:-1]
    if len(str) == 0:
        return []
    return list(map(int,str.split(',')))

def localmax(tss, ddict, g_len):
    inter = para.str_search_len
    te = min(g_len, tss + inter)
                
    depth = ddict[tss:te]
    window = []
    for i in range(0, inter-20):
        window.append(np.mean(depth[i:i+20]))
    derivative = np.gradient(window)
    # x = np.linspace(0,inter-20,inter-20)
    # plt.plot(x,window,'blue',label = 'window')
    # plt.plot(x,depth[0:inter-20],'green',label = 'depth')
    # plt.plot(x,derivative,'orange', label = 'derivative')
    maxw = max(window)
    maxi = np.where(window == maxw)[0][0]
    # plt.plot(maxi,maxw,'o')
    # plt.legend()
    # plt.show()
    # plt.close('all')
    if maxi >= inter-21:
        i = tss+maxi+1
        while np.mean(ddict[i-20:i]) > np.mean(ddict[i-1-20:i-1]):
            i+=1
        return max(np.mean(ddict[i-1-20:i-1]),maxw)
    # for i in range(0, len(window) - 1):
    #     if derivative[i] >= 0 and derivative[i+1] < 0 and window[i] == max(window[0:i+1]):
    #         plt.plot(i,window[i],'o')
    #         plt.legend()
    #         plt.show()
    #         plt.close('all')
    #         return window[i]
    # plt.legend()
    # plt.show()
    # plt.close('all')
    return maxw

def seq2onehot(seq):
    dic = {'A':[1,0,0,0], 'T':[0,1,0,0], 'C':[0,0,1,0], 'G':[0,0,0,1]}
    onehot = []
    for i in seq:
        onehot+= dic[i]
    return onehot

def correlation_xy(x,y):
    return np.corrcoef(x, y)[0,1]

def plt_depth_se_tss(s,e,tss,depth,start, name, algtss):
    if e<s:
        return
    x = np.linspace(s,e,e-s)
    plt.plot(x,depth[s:e],'blue',label = 'depth')
    plt.plot(tss,depth[tss],'o',label = 'dbtss')
    plt.plot(start,depth[start],'*',label = 'start')
    plt.plot(algtss,depth[algtss],'+',label = 'algtss')
    plt.legend()
    plt.savefig(para.FIG_DIR+name+'.png')
    plt.show()
    plt.close()
