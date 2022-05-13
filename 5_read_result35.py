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


def open_result():
	tss_dict=info.open_npy(info.TSS_DIR+inputf+'tss.npy')
	gdict = info.open_gdict_npy()
	promoter_dict = info.open_npy(info.PROMOTER_DIR+inputf+'promoter_dict.npy')

	fext= open(info.RESULT_DIR+inputf+'ext','r')
	f35= open(info.RESULT_DIR+inputf+'35','r')
	f35p= open(info.MOTIF_DIR+inputf+'35.txt','w')
	f35f = open(info.FASTA_DIR+inputf+'35logo.fasta','w')

	cf = open(info.TSS_DIR+inputf+'35.csv', 'w', newline='')
	fieldnames = ['geneid','srand','tss','m35','cm35']
	writer = csv.DictWriter(cf, fieldnames=fieldnames)
	writer.writeheader()

	m35list = []
	spacers=[]

	s = 0
	for item in fext:
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
			promoter_dict[srand][key].mextpvalue = info.get_pvalue(item[2])
		if s > 0:
			s+=1
			continue
	s=0
	for item in f35:
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
			spacer = 20-int(item[1])
			m35 = item[4]
			t = promoter_dict[srand][key].correcttss
			upp=t-12-spacer-7-20
			m35p=t-12-spacer-7
			up = gdict[srand][upp:upp+20]

			promoter_dict[srand][key].m35p = m35p
			promoter_dict[srand][key].m35seq = m35
			promoter_dict[srand][key].spacer = spacer
			promoter_dict[srand][key].mupseq = up
			promoter_dict[srand][key].m35pvalue = info.get_pvalue(item[2])

			m35list.append(m35)
			spacers.append(spacer)
			writer.writerow({'geneid':key,'srand':srand,'m35':m35,'cm35':gdict[srand][promoter_dict[srand][key].m35p:
				promoter_dict[srand][key].m35p+7]})
			f35p.write(item[4]+'\n')
			f35f.write('>'+item[0]+'\n'+m35+'\n')
		if s > 0:
			s+=1
			continue
	np.save(info.MOTIF_DIR+inputf+'m35list.npy',m35list)
	np.save(info.MOTIF_DIR+inputf+'spacers.npy',spacers)
	np.save(info.PROMOTER_DIR+inputf+'promoter_dict.npy',promoter_dict)
	result = pd.value_counts(spacers)
	print(result)
	sm=np.mean(spacers)
	info.get_pssm(m35list)
	info.gen_logo(inputf+'35')
	fext.close()
	f35.close()
	f35p.close()
	f35f.close()

if __name__ == '__main__':
	open_result()
