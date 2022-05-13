import numpy as np
import math
import scipy
from scipy.stats import binom
from Bio.Seq import Seq
from Bio import motifs
import matplotlib.pyplot as plt
from weblogo import *
import csv
import time
from Bio import motifs
from Bio.Seq import Seq

GTF_DIR = './00_genome/'
DEPTH_DIR = './01_depth/'
MEME_DIR = './02_meme/'
SCORE_DIR = './03_score/'
TSS_DIR = './04_tss/'
PROMOTER_DIR = './05_promoterseq/'
FASTA_DIR = './06_fasta/'
RESULT_DIR = './07_result/'
MOTIF_DIR = './09_motif/'
SIGMA_DIR = './10_sigma/'
GENE_DIR = './08_gene_data/'
LOGO_DIR = './11_logo/'


SIGMA_NUM = '70'
SRR967='SRR967'
SRR9651=['SRR1-965','+','./01_depth/SRR1173965_f.depth.txt']
SRR9652=['SRR1-965','-','./01_depth/SRR1173965_r.depth.txt']
SRR965S1=['SRR1-965S','+','./01_depth/SRR1173965_f.depth.txt']
SRR965S2=['SRR1-965S','-','./01_depth/SRR1173965_r.depth.txt']
SRR965S=['SRR1-965S','-','./01_depth/SRR1173965_f.depth.txt','./01_depth/SRR1173965_r.depth.txt']
GSE53767F='GSE53767F'
GSE53767R='GSE53767R'
SRR1='SRR1'


CDS2 = ['cds_protein','cds_protein_50','cds_protein_100','cds_protein_200','cds_protein_500','cds_protein_1000','cds_protein_max']
CDS = ['cds_protein','cds_protein1','cds_protein2','cds_protein3','cds_protein4']
FKPM1 = 10.87
FKPM2 = 23.23
FKPM3 = 49.77

terminus = 1600000
origin = 3900000
#cds fkpm 10.87, 23.23,49.77
# SRR = './03_srr/'
# TSS = './04_tss/'
# SEQ = './05_promoterseq/'

int_arr = SRR965S
inp = int_arr[0]
inp4 = inp[0:4]
ind = int_arr[2]
srand = int_arr[1]
indf = int_arr[2]
indr = int_arr[3]
indid = indf.split('_f')[0]
# inp = SRR967

gtf_file = './00_genome/GCF_000005845.2_ASM584v2_genomic.gtf'
g_dic = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
protein_start = ['ATG','GTG','TTG']
protein_end = ['TAA','TAG','TGA']

INF_MIN = float('-inf')
class GENESEQ():

    def __init__(self):
        #from gtf
        self.seqid = ''
        self.source = ''
        self.feature = ''
        self.start = -1
        self.end = -1
        self.len = 0
        self.srand = '+'
        self.phase = -1
        self.score = 0
        self.attributes = ''
        self.geneid = ''
        self.transcriptid = ''
        self.proteinid = ''
        self.genename=''
        #genome
        self.seq = ''
        #compute
        self.avedepth = 0
        self.strength = 0
        self.fpkm = 0
        #algorithm
        self.tss = -1
        self.sigma = 0
        self.correcttss = 0
        #motif position
        self.m10p = 0
        self.m35p = 0
        self.spacer = 0
        self.mext10p = 0
        #p-value
        self.m10pvalue = 0.0
        self.m35pvalue = 0.0
        self.mext10pvalue = 0.0
        #score
        self.m10score = 0.0
        self.m35score = 0.0
        self.mext10score = 0.0
        self.mspacerscore = 0.0
        self.mupscore = 0.0
        self.mscore = 0
        #motif seq
        self.m10seq = ''
        self.m35seq = ''
        self.mext10seq = ''
        self.mupseq = ''
        self.promoter = ''

    def set_setid(self, id):
        self.seqid = id
    
    def set_source(self, source):
        self.source = source

    def set_feature(self, feature):
        self.feature = feature

    def set_start(self, start, end, g_len):
        if self.srand == '+':
            self.start = start - 1
        else:
            self.start = g_len - end

    def set_end(self, start, end, g_len):
        if self.srand == '+':
            self.end = end - 1
        else:
            self.end = g_len - start
        if self.feature == 'CDS':
            self.end+=3

    def set_len(self):
        self.len = abs(self.end - self.start) + 1

    def set_score(self, score)    :
        if score != '.':
            self.score = float(score)
        else:
            self.score = 0
    
    def set_srand(self, srand):
        self.srand = srand

    def set_phase(self, phase):
        if self.feature == 'CDS':
            self.phase = int(phase)
        else:
            self.phase = phase
    
    def set_attributes(self, attributes):
        self.attributes = attributes
    
    def set_seq(self,gdict):
        self.seq = gdict[self.srand][self.start:self.end+1]
        # if self.srand == '+':
        #     self.seq = seq
        # else:
        #     self.seq = str(Seq(seq).reverse_complement())
        # self.seq = gdict[self.srand][self.start:self.end+1]

    def set_avedepth(self, depthf,depthr):
        if self.srand == 1:
            self.avedepth = np.mean(depthf[self.start:self.end+1])
        else:
            self.avedepth = np.mean(depthr[self.start:self.end+1])
    
    def set_geneid(self, geneid):
        self.geneid = geneid

    def set_transcriptid(self, transcriptid):
        self.transcriptid = transcriptid

    def set_proteinid(self, proteinid):
        self.proteinid = proteinid

    def set_fpkm(self,allreadsr,allreadsf):
        if self.srand == 1:
            # reads = np.sum(depthf[self.start:self.end+1])
            self.fpkm = pow(10,9) / allreadsf * self.avedepth
            # pow(10,9) / allreadsf * reads / self.len 
        else:
            # reads = np.sum(depthr[self.start:self.end+1])
            self.fpkm = pow(10,9) / allreadsr * self.avedepth
            # pow(10,9) / allreadsr * reads / self.len 
        # a = 0

def open_npy(input_file):
    data=np.load(input_file,allow_pickle=True)
    data = data.item()
    return data

def open_npy_list(input_file):
    data=np.load(input_file,allow_pickle=True)
    data = data.tolist()
    return data

def open_genome():
    genome = ''
    fg = './00_genome/GCF_000005845.2_ASM584v2_genomic.fna'
    f= open(fg,'r')
    i = 0
    for item in f:
        if i > 0:
            genome = genome + item[0:-1]
        i = i + 1
    f.close()
    print('genome len:', len(genome), end='\n')
    genome_seq = Seq(genome)
    # np.save(input_file+'.npy',genome)
    return genome_seq

def open_gdict_npy():
    return open_npy(GTF_DIR+'gdict.npy')

def open_genome_npy():
    return open_npy('./00_genome/genome')
    # if srand == '+':
    #     return open_npy('./00_genome/genome')
    # else:
    #     return open_npy('./00_genome/genomer')

def open_genomer_npy():
    return open_npy('./00_genome/genomer.npy')

def open_depth(input_file, g_len):
    dth = np.zeros((g_len))
    f= open(input_file,'r')
    for item in f:
        item = item.split()
        dth[int(item[1])-1]=float(item[2])
    f.close()
    np.save(input_file+'.npy',dth)
    return dth

def if_start(index, depth):
    r = 1
    readsi=depth[index]
    readsi_1= depth[index-1*r]
    readsi_2= depth[index-2*r]
    readsi_10 = np.mean(depth[max(index-10,0):index])
    cdf2 = binom.cdf(readsi_2, readsi_2 + readsi, 0.5)
    cdf1 = binom.cdf(readsi_1, readsi_1 + readsi, 0.5)
    return ((cdf1 < 0.01 and readsi / (readsi_1 + 0.01) > 2) 
        and (cdf2 < 0.01 and (readsi + 0.01)/ readsi_2 > 2)) or (readsi_10 == 0 and readsi > 10)

def if_start2(i, depth,copy_number,high,low):
    # high = np.mean(depth[i:i+10])
    # low = np.mean(depth[i-10:i])
    return high - low >20*copy_number and (high/low > 3) and (high>20*copy_number)

def copy_num(i):
    copy_number = 1
    if i < origin:
        copy_number = 1+3*abs(i-terminus)/(origin-terminus)
    else:
        copy_number = 4-3*(i-origin)/(origin-terminus)
    return copy_number

def plot_test(dth,plt_s,plt_e):
    x = np.linspace(plt_s,plt_e,plt_e-plt_s)
    plt.plot(x,dth[plt_s:plt_e])
    plt.show()
    plt.close("all")

def gen_fasta(id):
    fp = open(PROMOTER_DIR+id+'promoter.txt','r')
    fa = open(FASTA_DIR+id+'promoter.fa','w')
    i = 0
    for item in fp:
        fa.write('>'+str(i)+'\n'+item)
        i+=1
    fa.close()

def gen_logo(id):
    fin = open(FASTA_DIR+id+'.fasta')
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)

    options = LogoOptions() 
    options.fineprint = 'Python plot'   # 生成下角标
    format = LogoFormat(data, options)

    eps = eps_formatter(data, format)

    with open(LOGO_DIR+id+'.eps','wb') as f:
        f.write(eps)

def gen_logo_logo(id):
    fin = open(FASTA_DIR+id+'logo.fasta')
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)

    options = LogoOptions() 
    options.fineprint = 'Python plot'   # 生成下角标
    format = LogoFormat(data, options)

    eps = eps_formatter(data, format)

    with open(LOGO_DIR+id+'.eps','wb') as f:
        f.write(eps)
    
def gen_fasta_s_e(s,e,id,kind):
    data=np.load(PROMOTER_DIR+id+'promoter_dict.npy',allow_pickle=True)
    promoters = data.item()
    fa = open(FASTA_DIR+id+kind+'.fasta','w+')
    for srand, plist in promoters.items():
        for i in plist.keys():    
            fa.write('>'+i+'.'+srand+'\n')
            fa.write(plist[i].promoter[s:e]+'\n')
    fa.close()

def get_pvalue(str):
    plist = str.split('e')
    num = float(plist[0])
    e = int(plist[1])
    pvalue = pow(num,e)
    return pvalue

def get_pssm(list):
    i = []
    for item in list:
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
    return pssm

def get_motif(list):
    i = []
    for item in list:
        i.append(Seq(item))
    m = motifs.create(i)
    return m

def motif2pssm(m):
    pwm = m.counts.normalize(pseudocounts=0.5)
    print(pwm)
    pssm = pwm.log_odds()
    for key,value in pssm.items():
        for i in range(0, len(value)):
            if value[i] <-2.37:
                pssm[key][i] = -2.37
    print(pssm)
    return pssm
    
def get_pwm(list):
    i = []
    for item in list:
        i.append(Seq(item))
    m = motifs.create(i)
    pwm = m.counts.normalize(pseudocounts=0.5)
    # print(pwm)
    # pssm = pwm.log_odds()
    # print(pssm)
    # pwm_log = {}
    # for key,value in pwm.items():
    #     pwm_log[key]=[]
    #     for i in range(0,len(value)):
    #         pwm_log[key].append(math.log2(value[i]/0.25))
    
    # print(pwm_log)
    return pwm

def cal_pwmscore(pwm, seq):
    score = 0
    for i in range(0, pwm):
        score+=math.log2(pwm[seq[i]][i]/0.25)
    return score

def cal_score(pssm,seq):
    score = 0
    x=pssm.calculate(Seq(seq))
    for i in range(0,len(seq)):
        score+=pssm[seq[i]][i]
    return score

def plot_bin(data):
    #以batch为间隔统计，并且x轴与箱子对齐
    batch = 1
    bands = int((max(data)-min(data))/batch)+1
    bins = np.arange(0,int(max(data))+1,batch)
    print(bins)
    plt.hist(data,bins,rwidth=0.5)
    x_ticks=bins
    print(x_ticks)
    plt.xticks(x_ticks)
    plt.xlim(0,int(max(data))+1)
    plt.show()

def normalize(list):
    res = list
    s = sum(list)
    for i in range(0, len(list)):
        res[i] = list[i]/s
    return res

def log_ods(list):
    res = {}
    s = sum(list)
    l = float(len(list))
    for key,value in list.items():
        p = value/s
        res[key] = math.log2(p*l)
    return res