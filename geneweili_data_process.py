from cmath import inf, pi
import imp
from sys import api_version
from typing import Iterator
from xml.sax.xmlreader import AttributesImpl
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.function_base import average
from numpy.lib.utils import deprecate, deprecate_with_doc
import pandas as pd
from pandas.core.construction import is_empty_data
from scipy.stats import binom
import csv
import logomaker
import os
import time
from Bio import motifs
from Bio.Seq import Seq

import utils
import input_process as ip
import promoter_process as pp
import parameters as para

if __name__ == '__main__':
    gdict = utils.open_gdict_npy()
    g_len = len(gdict['+'])
    id = para.sid
    gene_dict = ip.gen_gen_dict(id, gdict)
    paper_dict = ip.gen_paper_dict(gene_dict)
    ddict = ip.gen_ddict(id, g_len)
    gene_dict, pdict = pp.potential_promoter_cdf(ddict, gdict, gene_dict, id, g_len)
    # debug
    # pdict = utils.open_npy(para.PROMOTER_DIR + id + '-raw-pdict.npy')
    # gene_dict = utils.open_npy(para.GENE_DIR + id +'gene_after_tss.npy')
    # plist = open(para.PROMOTER_DIR + id + 'ppromoter.txt',"r").read().split()
    # slist = open(para.PROMOTER_DIR + id + 'spromoter.txt',"r").read().split()
    # ilist = open(para.PROMOTER_DIR + id + 'ipromoter.txt',"r").read().split()
    # alist = open(para.PROMOTER_DIR + id + 'apromoter.txt',"r").read().split()
    # utils.plot_logo(plist, id+'p')
    # utils.plot_logo(slist, id+'s')
    # utils.plot_logo(ilist, id+'i')
    # utils.plot_logo(alist, id+'a')

    # fieldnames = ['gid', 'srand','start', 'end', 'ptss', 'stss', 'itss', 'astss']
    # writer = utils.csv_writer(para.TSS_DIR+id+'tss.csv', fieldnames)
    # for gid, iseq in gene_dict.items():
    #     writer.writerow({'gid':gid, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'ptss':iseq.ptss, 'stss':iseq.stss,
    #                      'itss':iseq.itss, 'astss':iseq.astss})

    utils.gen_fasta(para.s10, para.e10, id, '10', pdict)
    pp.meme(id, 6, '10raw')
    pdict, m10list, mextlist = pp.open_result10(gdict, pdict, id, paper_dict)
    pp.meme(id, 6, 'extmeme')
    pp.meme(id, 7, '35raw')
    pdict, m35list, spacers = pp.open_result35(gdict, pdict, id)
    pdict, m10list, mextlist, m35list, spacers = pp.iteration_base(gdict, pdict, id, 'p', m10list, mextlist, m35list, spacers)
    # debug
    # m10list = utils.open_npy_list(para.MOTIF_DIR+id+'m10ilist.npy')
    # mextlist = utils.open_npy_list(para.MOTIF_DIR+id+'mextilist.npy')
    # m35list =  utils.open_npy_list(para.MOTIF_DIR+id+'m35ilist.npy')
    # pdict = utils.open_npy(para.PROMOTER_DIR+id+'ppdict-iteration.npy')
    # spacers = utils.open_npy_list(para.MOTIF_DIR+id+'spacersi.npy')
    for kind in para.klist:
        pdict = pp.pwm_tss(gdict, pdict, id, kind, m10list, mextlist, m35list, spacers)
    
    #debug
    # pdict = utils.open_npy(para.PROMOTER_DIR+id+'apdict-pwmtss.npy')
    pltdict = {'p': 'o', 'i': 'x', 's':'*'}
    
    for gid, pallitem in pdict.items():
        if len(pallitem['p']) <= 0:
            continue
        paper_item = paper_dict[gid]
        iseq = gene_dict[gid]

        srand = iseq.srand
        start = iseq.start
        end = iseq.end
        depth = ddict[srand]
        plots = max(0, start - 300)
        ptss = pallitem['p'][0].ctss
        depth_plot = depth[ptss : start + 100].astype(np.int).tolist()
        maxd = max(depth_plot)
        maxi = depth_plot.index(max(depth_plot)) + ptss
        plote = maxi + 50
        x = np.linspace(plots,plote,plote-plots)
        plt.plot(x,depth[plots:plote])
        plt.plot()
        for kind in para.plot_klist:
            for pitem in pallitem[kind]:
                ctss = pitem.ctss
                dtss = depth[ctss]
                plt.plot(ctss,depth[ctss],pltdict[kind])
                plt.text(ctss, dtss+3, kind+'('+str(ctss)+','+str(dtss)+')', fontsize=10)
            for pitem in paper_item[kind]:
                ctss = pitem.ctss
                dtss = depth[ctss]
                plt.plot(ctss,depth[ctss],pltdict[kind])
                plt.text(ctss, dtss+3, kind+'paper('+str(ctss)+','+str(dtss)+')', fontsize=10)
        plt.plot(maxi,maxd,'+')
        plt.text(maxi, maxd+3, 'max('+str(maxi)+','+str(maxd)+')', fontsize=10)
        plt.savefig(para.FIG_DIR+'/'+id+'/'+gid+'.png')
        # plt.show()
        plt.close("all")

    print('end')