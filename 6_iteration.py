from calendar import prmonth
from os import close, startfile
from socket import SocketIO
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

# 967:73-87

def parse_args():
    parser = argparse.ArgumentParser(description='read promoters',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fasta', type=str, default='./promoterdata/PromoterSigma', help='sigma promoter file')
    parser.add_argument('--num', type=int, default=70, help='sigma number')
    # parser.add_argument('--meme', type=str, default=info.SIGMA70, help='meme result file')
    parser.add_argument('--genome_file', type=str, default='./00_genome/GCF_000005845.2_ASM584v2_genomic.fna', help='genome file')
    parser.add_argument('--gtf_file', type=str, default='./00_genome/GCF_000005845.2_ASM584v2_genomic.gtf', help='genome file')	
    parser.add_argument('--depth_file', type=str, default='./01_depth/xxx.txt', help='depth file')
    # params
    args = parser.parse_args()
    return args

def open_memeresultext(tss,genome,promoters):
	f= open(info.RESULT_DIR+inputf+inputs+'m10+extresult','r')
	# open(info.RESULT_DIR+inputf+'cdsresult2','r')
	fx=xlwt.Workbook()
	sheet1 = fx.add_sheet('all',cell_overwrite_ok=True)
	row = [u'geneid',u'corrected tss',u'motif10',u'corrected motif10','ext10']
	for j in  range(0,len(row)):
		sheet1.write(0,j,row[j])
	
	now = time.strftime("%Y-%m-%d-%H_%M_%S",time.localtime(time.time())) 
	
	s = 0
	ti=1
	for item in f:
		if 'Start' in item:
			s = 1
			continue
		if s > 1:
			if '------' in item:
				break
			item = item.split()
			tss_original = int(item[0])
			key = item[0]
			ext = item[4][0:4]
			m10 = item[4][5:11]
			spacer = int(item[1])
			if inputf == info.SIGMA_NUM:
				tss[item[0]]=tss[item[0]]+int(item[1])-4
			elif inputf == info.SRR967:
				tss[item[0]]=tss[item[0]]+int(item[1])-16
			elif inputf4 == info.SRR1:
				promoters[key].correcttss = tss_original+spacer-14
				t = promoters[key].correcttss
				promoters[key].m10p = t-12
				promoters[key].m10seq = m10
				promoters[key].mext10p = t - 17
				promoters[key].mext10seq = ext
				if ext == 'CTGT':
					a = 0
			row = [tss_original, t,m10,genome[t-12:t-6],ext]
			for j in  range(0,len(row)):
				sheet1.write(ti,j,row[j])
			ti+=1
		if s > 0:
			s+=1
			continue
	now = time.strftime("%Y-%m-%d-%H_%M_%S",time.localtime(time.time())) 
	fx.save(info.TSS_DIR+now+inputf+inputs+'corectedtss.xls')
	np.save(info.TSS_DIR+inputf+info.srand+'tss_correct10.npy',tss)
	np.save(info.PROMOTER_DIR+inputf+inputs+'promoter_dict_corecttss.npy',promoters)
	return tss,promoters

def compare(i1,i2):
	l1 = info.get_pwm(i1)
	l2 = info.get_pwm(i2)
	for key, value in l1.items():
		for i in range(0,len(value)):
			if abs(l1[key][i] - l2[key][i]) > 0.02:
				return False
	return True

if __name__ == '__main__':
	gdict = info.open_gdict_npy()
	promoter_dict = info.open_npy(info.PROMOTER_DIR+inputf+'promoter_dict.npy')
	m10list = info.open_npy_list(info.MOTIF_DIR+inputf+'m10list.npy')
	mextlist = info.open_npy_list(info.MOTIF_DIR+inputf+'mextlist.npy')
	m35list = info.open_npy_list(info.MOTIF_DIR+inputf+'m35list.npy')
	pssm10 = info.get_pssm(m10list)
	pssmext = info.get_pssm(mextlist)
	pssm35 = info.get_pssm(m35list)
	# pwm = info.get_pwm(m10list)
	
	cnt = 0
	spacers = info.open_npy_list(info.MOTIF_DIR+inputf+'spacers.npy')
	while 1:
		cntspacers = pd.value_counts(spacers)
		print(cntspacers)
		lspacers = info.log_ods(cntspacers)
		cf = open(info.PROMOTER_DIR+'all.csv', 'w', newline='')
		fieldnames = ['algtss','tss','10','ext','spacer','35','up','m10score','mextscore','m35score','spacerscore','promoter']
		writer = csv.DictWriter(cf, fieldnames=fieldnames)
		writer.writeheader()
		# for srand, plist in promoter_dict.items():
		# 	for key, iseq in plist.items():
		# 		tss = int(key)
		# 		ctss = iseq.correcttss
		# 		writer.writerow({'algtss':key,'tss':ctss,'10':iseq.m10seq,'ext':iseq.mext10seq,'spacer':gdict[srand][ctss-12-iseq.spacer:ctss-12],
		# 			'35':iseq.m35seq,'up':gdict[srand][ctss-12-iseq.spacer-20:ctss-12-iseq.spacer],'m10score':iseq.m10score,'mextscore':
		# 			iseq.mext10score,'m35score':iseq.m35score,'spacerscore':lspacers[iseq.spacer],'promoter':gdict[srand][tss-100:tss]})

		
		
		m10ilist = []
		mextilist = []
		m35ilist = []
		spacers = []
		pcnt = 0
		for srand, plist in promoter_dict.items():
			for key, iseq in plist.items():
				pcnt+=1
				tss = int(key)
				ti = 0
				tj = 0
				mscore = info.INF_MIN
				for i in range(5,20):
					ctss = tss-i+6
					if ctss == 0:
						a = 0
					m10 = gdict[srand][ctss-12:ctss-6]
					mext = gdict[srand][ctss-17:ctss-13]
					m10score = pssm10.calculate(Seq(m10))
					mextscore = pssmext.calculate(Seq(mext))
					m35sscore = info.INF_MIN
					tm35=''
					tspacer = 0
					tm35score=info.INF_MIN
					tsscore=info.INF_MIN
					for j in range(15,20):
						m35 = gdict[srand][ctss-12-j-7:ctss-12-j]
						c35score = pssm35.calculate(Seq(m35))
						csscore = lspacers[j]
						c35sscore = c35score+csscore
						# pvalue = motif.instances[0].pvalue
						if c35sscore>m35sscore:
							tj = j
							m35sscore=c35sscore
							tm35 = m35
							tspacer = j
							tm35score = c35score
							tsscore = csscore
					ctotal = float(m10score) + float(mextscore) + float(tm35score)+float(tsscore)
					if ctotal>mscore: #and (m10score>m35score or promoter_dict[srand][key].correcttss ==0):
						ti = i
						mscore = ctotal
						up = gdict[srand][ctss-12-iseq.spacer-7-20:ctss-12-iseq.spacer-7]

						promoter_dict[srand][key].m10seq = m10
						promoter_dict[srand][key].mext10seq = mext
						promoter_dict[srand][key].correcttss = ctss
						promoter_dict[srand][key].m10p = ctss-12
						promoter_dict[srand][key].mextp = ctss-17
						promoter_dict[srand][key].m10score = m10score
						promoter_dict[srand][key].mext10score = mextscore
						promoter_dict[srand][key].m35seq = tm35
						promoter_dict[srand][key].spacer = tspacer
						promoter_dict[srand][key].m35score = tm35score
						promoter_dict[srand][key].mspacerscore = tsscore
						promoter_dict[srand][key].mupseq = up

						iseq = promoter_dict[srand][key]
						m10 = iseq.m10seq
						m35 = iseq.m35seq
						mspacer = gdict[srand][ctss-12-iseq.spacer:ctss-12]
						m35s10 = m35+mspacer+m10
						target = gdict[srand][ctss-12-iseq.spacer-7:ctss-6]
						if m35s10 != target:
							a=0
						
		f10f= open(info.FASTA_DIR+inputf+'10ilogo.fasta','w')
		fef= open(info.FASTA_DIR+inputf+'extilogo.fasta','w')
		fef2= open(info.FASTA_DIR+inputf+'exti.fasta','w')
		f35f = open(info.FASTA_DIR+inputf+'35ilogo.fasta','w')
		fupf = open(info.FASTA_DIR+inputf+'uplogo.fasta','w')

		

		pcnt = 0
		for srand, plist in promoter_dict.items():
			for key, iseq in plist.items():
				tss = int(key)
				ctss = iseq.correcttss
				# up = gdict[srand][ctss-12-iseq.spacer-7-20:ctss-12-iseq.spacer-7]
				f10f.write('>'+key+'.'+srand+'\n'+iseq.m10seq+'\n')
				fef.write('>'+key+'.'+srand+'\n'+iseq.mext10seq+'\n')
				fef2.write('>'+key+'.'+srand+'\n'+'AA'+iseq.mext10seq+'\n')
				f35f.write('>'+key+'.'+srand+'\n'+iseq.m35seq+'\n')
				fupf.write('>'+key+'.'+srand+'\n'+iseq.mupseq+'\n')
				
				m10ilist.append(iseq.m10seq)
				mextilist.append(iseq.mext10seq)
				m35ilist.append(iseq.m35seq)
				spacers.append(iseq.spacer)

				m10 = iseq.m10seq
				m35 = iseq.m35seq
				mspacer = gdict[srand][ctss-12-iseq.spacer:ctss-12]
				m35s10 = m35+mspacer+m10
				target = gdict[srand][ctss-12-iseq.spacer-7:ctss-6]
				if m35s10 != target:
					a=0
				
				writer.writerow({'algtss':key,'tss':ctss,'10':iseq.m10seq,'ext':iseq.mext10seq,'spacer':gdict[srand][ctss-12-iseq.spacer:ctss-12],
					'35':iseq.m35seq,'up':iseq.mupseq,'m10score':iseq.m10score,'mextscore':iseq.mext10score,'m35score':iseq.m35score,
					'spacerscore':iseq.mspacerscore,'promoter':gdict[srand][tss-100:tss]})
				pcnt+=1
		info.gen_logo_logo(inputf+'10i')
		info.gen_logo_logo(inputf+'exti')
		info.gen_logo_logo(inputf+'35i')

		# motif10 = info.get_motif(m10ilist)
		# pvalue = motif10.instances[0].pvalue
		# motifext = info.get_motif(mextilist)
		# motif35 = info.get_motif(m35ilist)
		pssm10 = info.get_pssm(m10ilist)
		pssmext = info.get_pssm(mextilist)
		pssm35 = info.get_pssm(m35ilist)
		print(np.mean(spacers))
		info.plot_bin(spacers)
		f10f.close()
		fef.close()
		fef2.close()
		f35f.close()
		fupf.close()
		if compare(m10ilist, m10list) and compare(mextilist, mextlist) and compare(m35ilist, m35list):
			np.save(info.MOTIF_DIR+inputf+'m10ilist.npy',m10ilist)
			np.save(info.MOTIF_DIR+inputf+'mextilist.npy',mextilist)
			np.save(info.MOTIF_DIR+inputf+'m35ilist.npy',m35ilist)
			np.save(info.PROMOTER_DIR+inputf+'promoter_dict.npy',promoter_dict)
			np.save(info.MOTIF_DIR+inputf+'spacers.npy',spacers)
			break
		m10list = m10ilist
		mextlist = mextilist
		m35list = m35ilist	
		cnt+=1
#todo:提取pvalue，pssm，以及pwm计算多重不同的方式。做regression