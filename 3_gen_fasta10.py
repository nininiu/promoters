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
from weblogo import *



# 967:73-87
inputf = info.inp
inputf4 = info.inp4



def gen_fasta10(srand):
	# input_file=info.PROMOTER_DIR+input+'promoter.npy'
	data=np.load(info.PROMOTER_DIR+inputf+srand+'promoter_dict.npy',allow_pickle=True)
	promoters = data.item()
	fa = open(info.FASTA_DIR+inputf+srand+'m10+ext.fasta','w+')
	for i in promoters.keys():	
		fa.write('>'+i+'\n')
		tmp = promoters[i]
		if inputf == info.SIGMA_NUM:
			fa.write(promoters[i][46:58]+'\n')
		if inputf == info.SRR967:
			fa.write(promoters[i][73:98]+'\n')
		if inputf4 == info.SRR1:
			fa.write(promoters[i].promoter[70:95]+'\n')
	fa.close()

def open_sigmasrr():
	data=np.load(info.TSS_DIR+inputf+'corecttss.npy',allow_pickle=True)
	data = data.item()
	return data

def gen_fasta35(genome):
	# input_file=info.PROMOTER_DIR+input+'promoter.npy'
	# data=np.load(info.PROMOTER_DIR+input+'promoter.npy',allow_pickle=True)
	promoters = info.open_npy(info.TSS_DIR+inputf+'cdscorecttss')
	# open_sigmasrr()
	fa = open(info.FASTA_DIR+inputf+'-cds35.fasta','w+')
	for i in promoters.keys():	
		
		tmp = promoters[i]
		if inputf == info.SIGMA_NUM:
			fa.write(promoters[i][46:58]+'\n')
		if inputf == info.SRR967:
			fa.write(promoters[i][73:98]+'\n')
		if inputf == info.SRR1:
			tss=promoters[i].correcttss
			if tss !=0:
				fa.write('>'+i+'\n')
				fa.write(genome[tss-37:tss-28]+'\n')
	fa.close()

def gen_fastafimo(genome):
	# input_file=info.PROMOTER_DIR+input+'promoter.npy'
	# data=np.load(info.PROMOTER_DIR+input+'promoter.npy',allow_pickle=True)
	promoters = info.open_npy(info.TSS_DIR+inputf+'cdscorecttss')
	# open_sigmasrr()
	fa = open(info.FASTA_DIR+inputf+'-fimo.fasta','w+')
	for i in promoters.keys():	
		tss=promoters[i].tss
		if tss !=0:
			fa.write('>'+i+'\n')
			fa.write(genome[tss-100:tss]+'\n')		
	fa.close()

def open_motif(n):
	input_file = info.MOTIF_DIR+inputf+'-'+str(n)+'.txt'
	f= open(input_file,'r')
	P = []
	PWM = []
	for item in f:
		p = []
		if '#' not in item:
			item = item.split()
			i = 1
			while i < 5:
				p.append(int(item[i]))
				i+=1
			P.append(p)
			s = sum(p)

			pwm = np.array(p)/s
			PWM.append(pwm)
	f.close()
	np.save(info.MOTIF_DIR+inputf+'-'+str(n)+'.npy',PWM)

	return PWM

def open_srr():
	data=np.load(info.TSS_DIR+inputf+'corecttss.npy',allow_pickle=True)
	data = data.item()
	# tmp = data[10:30]
	# l = len(data)
	return data
	
if __name__ == '__main__':
	# genome = info.open_genome()
	info.gen_fasta_s_e(78,95,inputf,'10')
	info.gen_logo(inputf+'10')
	# os.system("ping 192.168.1.101")
