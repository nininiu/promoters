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
import csv
import time
import info
from Bio import motifs
from Bio.Seq import Seq
inputf = info.inp
inputs = info.srand
inputf4 = info.inp4

def open_result10(id):
	tss_dict=info.open_npy(info.TSS_DIR+inputf+'tss.npy')
	gdict = info.open_gdict_npy()
	promoter_dict = info.open_npy(info.PROMOTER_DIR+inputf+'promoter_dict.npy')

	f= open(info.RESULT_DIR+inputf+id,'r')
	fp= open(info.MOTIF_DIR+inputf+id+'.txt','w')
	fe= open(info.MOTIF_DIR+inputf+'ext.txt','w')
	f10f= open(info.FASTA_DIR+inputf+id+'logo.fasta','w')
	fef= open(info.FASTA_DIR+inputf+'extlogo.fasta','w')
	f35f = open(info.FASTA_DIR+inputf+'35.fasta','w')

	cf = open(info.TSS_DIR+inputf+id+'.csv', 'w', newline='')
	fieldnames = ['geneid','srand','tss','m10','cm10','mext']
	writer = csv.DictWriter(cf, fieldnames=fieldnames)
	writer.writeheader()

	m10list= []
	mextlist = []

	s = 0
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
			tss_original = int(keylist[0])
			srand = keylist[1]
			m10 = item[4]
			inter = int(item[1])
			promoter_dict[srand][key].correcttss = tss_original + inter - 11
			tss = promoter_dict[srand][key].correcttss
			promoter_dict[srand][key].m10p = tss-12
			promoter_dict[srand][key].m10pvalue = info.get_pvalue(item[2])
			promoter_dict[srand][key].m10seq = m10
			promoter_dict[srand][key].mext10p = tss - 17
			promoter_dict[srand][key].mext10seq = gdict[srand][promoter_dict[srand][key].mext10p:promoter_dict[srand][key].mext10p+4]

			m10list.append(m10)
			mextlist.append(promoter_dict[srand][key].mext10seq)
			writer.writerow({'geneid':key,'srand':srand,'tss':tss,'m10':m10,'cm10':gdict[srand][promoter_dict[srand][key].m10p:
				promoter_dict[srand][key].m10p+6],'mext':promoter_dict[srand][key].mext10seq})
			fp.write(item[4]+'\n')
			f10f.write('>'+str(item[0])+'\n'+item[4]+'\n')
			fe.write(promoter_dict[srand][key].mext10seq+'\n')
			fef.write('>'+str(item[0])+'\n'+'AA'+ promoter_dict[srand][key].mext10seq+'\n')
			p35 = tss-12-19-7
			m35 = gdict[srand][p35:p35+11]
			f35f.write('>'+str(item[0])+'\n'+m35+'\n')
		if s > 0:
			s+=1
			continue
	np.save(info.MOTIF_DIR+inputf+'m10list.npy',m10list)
	np.save(info.MOTIF_DIR+inputf+'mextlist.npy',mextlist)
	np.save(info.PROMOTER_DIR+inputf+'promoter_dict.npy',promoter_dict)
	info.get_pssm(m10list)
	info.get_pssm(mextlist)
	# m10m = motifs.create(m10i)
	# mextm = motifs.create(mexti)
	# print(m10m.counts)
	# print(mextm.counts)
	# m10pwm = m10m.counts.normalize(pseudocounts=0.5)
	# mextpwm = mextm.counts.normalize(pseudocounts=0.5)
	# print(m10pwm)
	# print(mextpwm)
	# m10m.weblogo("m10m.png")
	# m10pssm = m10pwm.log_odds()
	# print(print(m10pssm))
	# test_seq=Seq("TATAAT")
	# x=m10pssm.calculate(test_seq)
	info.gen_logo(inputf+id)
	info.gen_logo(inputf+'ext')
	fp.close()
	f.close()

if __name__ == '__main__':
	open_result10('10')
	

	