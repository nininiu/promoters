from cmath import inf
from os import close, startfile
from turtle import st
from typing import Iterator
from xml.sax.xmlreader import AttributesImpl
import numpy as np
import argparse
import matplotlib.pyplot as plt
from numpy.lib.function_base import average
from numpy.lib.utils import deprecate, deprecate_with_doc
import pandas as pd
import pdb
from gtfparse import read_gtf
from pandas.core.construction import is_empty_data
import scipy
from scipy.stats import binom
import xlwt
import time
import info

inputf = info.inp

def open_gtf(gtf_file,gdict):
	g_len = len(gdict['+'])
	cds = []
	gene = []
	cds_protein = {'+':[],'-':[]}
	gene_protein = {'+':[],'-':[]}

	f=xlwt.Workbook()
	sall = f.add_sheet('all',cell_overwrite_ok=True)
	sCDS = f.add_sheet('CDS',cell_overwrite_ok=True)
	sgene = f.add_sheet('gene',cell_overwrite_ok=True)
	sCDS_protein = f.add_sheet('CDS_protein',cell_overwrite_ok=True)
	sgene_protein = f.add_sheet('gene_protein',cell_overwrite_ok=True)
	stranscript = f.add_sheet('transcript',cell_overwrite_ok=True)
	sexon = f.add_sheet('exon',cell_overwrite_ok=True)
	stRNA = f.add_sheet('tRNA',cell_overwrite_ok=True)
	srRNA = f.add_sheet('rRNA',cell_overwrite_ok=True)
	rowall=['feature','start','end','score','srand','phase','seq','len',
	'geneid','transcriptid','proteinid','attributes']
	for i in range(0,len(rowall)):
		sall.write(0,i,rowall[i])  # 顺序为x行x列写入第x个元素
		sCDS.write(0,i,rowall[i]) 
		sgene.write(0,i,rowall[i])
		sCDS_protein.write(0,i,rowall[i])
		sgene_protein.write(0,i,rowall[i])
		stranscript.write(0,i,rowall[i])
		stranscript.write(0,i,rowall[i])
		sexon.write(0,i,rowall[i])
		stRNA.write(0,i,rowall[i])
		srRNA.write(0,i,rowall[i])
	rall = rCDS = rgene = rCDS_protein = rgene_protein = rtranscript = rexon = rtRNA = rrRNA= 1
	i=-1
	stop_codon = ''
	start_codon = ''

	data = np.loadtxt(gtf_file, dtype=np.str, delimiter='\t').tolist()
	for item in data:
		if len(item) > 8:
			feature= item[2]
			start = int(item[3])
			end = int(item[4])
			score = item[5]
			srand = item[6]
			if srand == '-':
				a = 0
			phase = item[7]
			attributes = item[8]
			att_list = item[8].split()

			iseq = info.GENESEQ()
			iseq.set_srand(srand)
			iseq.set_feature(feature)
			iseq.set_setid(item[0])
			iseq.set_source(item[1])
			iseq.set_start(start,end,g_len)
			iseq.set_end(start,end,g_len)
			iseq.set_len()
			iseq.set_score(score)
			iseq.set_phase(phase)
			iseq.set_attributes(attributes)
			iseq.set_seq(gdict)
			for atti in range(0, len(att_list)):
				if att_list[atti] == 'gene_id':
					iseq.set_geneid(att_list[atti+1][1:-2])
				elif att_list[atti] == 'transcript_id':
					iseq.set_transcriptid(att_list[atti+1][1:-2])
				elif att_list[atti] == 'protein_id':
					iseq.set_proteinid(att_list[atti+1][1:-2])
			rowall=[iseq.feature,iseq.start,iseq.end,iseq.score,iseq.srand,iseq.phase,iseq.seq,iseq.len,iseq.avedepth,iseq.fpkm,iseq.geneid,iseq.transcriptid,iseq.proteinid,iseq.attributes]
			if iseq.feature not in ['start_codon','stop_codon']:
				for j in  range(0,len(rowall)):
					sall.write(rall,j,rowall[j])
				rall+=1

			if (iseq.feature == 'CDS'):
				for j in  range(0,len(rowall)):
					sCDS.write(rCDS,j,rowall[j])
				rCDS+=1
				cds.append(iseq)
				
			if (iseq.feature == 'gene'):
				for j in  range(0,len(rowall)):
					sgene.write(rgene,j,rowall[j])  
				rgene+=1
				gene.append(iseq)

			if (iseq.feature == 'transcript'):
				for j in  range(0,len(rowall)):
					stranscript.write(rtranscript,j,rowall[j])  
				rtranscript+=1

			if (iseq.feature == 'exon'):
				for j in  range(0,len(rowall)):
					sexon.write(rexon,j,rowall[j])  
				rexon+=1

			if ('tRNA' in iseq.attributes):
				for j in  range(0,len(rowall)):
					stRNA.write(rtRNA,j,rowall[j])  # 顺序为x行x列写入第x个元素
				rtRNA+=1

			if ('rRNA' in iseq.attributes):
				for j in  range(0,len(rowall)):
					srRNA.write(rrRNA,j,rowall[j])  # 顺序为x行x列写入第x个元素
				rrRNA+=1

			if iseq.feature == 'start_codon':
				start_codon = iseq.seq
			if iseq.feature == 'stop_codon':
				stop_codon = iseq.seq
				if srand == '-':
					a = 0
				if start_codon in info.protein_start and stop_codon in info.protein_end:
					lastgene = gene[-1]
					lastcds = cds[-1]
					rowgene=[lastgene.feature,lastgene.start,lastgene.end,lastgene.score,lastgene.srand,lastgene.phase,lastgene.seq,lastgene.len,lastgene.avedepth,lastgene.fpkm,lastgene.geneid,
						lastgene.transcriptid,lastgene.proteinid,lastgene.attributes]
					rowcds=[lastcds.feature,lastcds.start,lastcds.end,lastcds.score,lastcds.srand,lastcds.phase,lastcds.seq,lastcds.len,lastcds.avedepth,lastcds.fpkm,lastcds.geneid,
						lastcds.transcriptid,lastcds.proteinid,lastcds.attributes]
					for j in  range(0,len(rowall)):
						sCDS_protein.write(rCDS_protein,j,rowcds[j])
						sgene_protein.write(rgene_protein,j,rowgene[j])
					rCDS_protein+=1
					rgene_protein+=1
					cds_protein[srand].append(lastcds)
					gene_protein[srand].append(lastgene)
	now = time.strftime("%Y-%m-%d-%H",time.localtime(time.time())) 
	f.save(info.GENE_DIR+now+inputf+'allinfo.xls')
	np.save(info.GENE_DIR+inputf+'gene_protein.npy',gene_protein)
	return gene_protein

if __name__ == '__main__':
	genome = info.open_genome()
	genomer = genome.reverse_complement()
	l = len(genome)

	np.save(info.GTF_DIR+'genome.npy',str(genome))
	np.save(info.GTF_DIR+'genomer.npy',str(genomer))
	gdict = {'+':str(genome),'-':str(genomer)}
	np.save(info.GTF_DIR+'gdict.npy',gdict)
	tmp = info.open_npy(info.GTF_DIR+'gdict.npy')

	open_gtf(info.gtf_file,gdict)

	