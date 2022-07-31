from cProfile import label
from cmath import inf, pi
import imp
from operator import itemgetter
from sys import api_version
from typing import Iterator
from xml.sax import make_parser
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
import out_process as op
import models as ml

if __name__ == '__main__':
    id = para.sid
    gdict = ip.get_gdict()
    g_len = len(gdict['+'])
    short_id = para.short_id
    gene_dict = ip.parse_gff(id, gdict)
    paper_dict = ip.gen_paper_dict(gene_dict)

    ddict = ip.gen_ddict(id, g_len)
    sddcit = ip.gen_ddict(short_id, g_len, suffix=para.short_suffix, inter=para.short_dp)

    # #debug
    # ddict = utils.open_npy(para.DEPTH_DIR + id +'ddict.npy')
    # sddcit = utils.open_npy(para.DEPTH_DIR + short_id +'ddict.npy')
    gene_dict, pdict = pp.potential_promoter_cdf(ddict, gdict, gene_dict, id, g_len)

    # debug
    # pdict = utils.open_npy(para.PROMOTER_DIR + id + '-raw-pdict.npy')
    # gene_dict = utils.open_npy(para.GENE_DIR + id +'gene_after_tss.npy')
    # # plist = open(para.PROMOTER_DIR + id + 'ppromoter.txt',"r").read().split()
    # # slist = open(para.PROMOTER_DIR + id + 'spromoter.txt',"r").read().split()
    # # ilist = open(para.PROMOTER_DIR + id + 'ipromoter.txt',"r").read().split()
    # # alist = open(para.PROMOTER_DIR + id + 'apromoter.txt',"r").read().split()
    # # utils.plot_logo(plist, id+'p')
    # # utils.plot_logo(slist, id+'s')
    # # utils.plot_logo(ilist, id+'i')
    # # utils.plot_logo(alist, id+'a')

    # # fieldnames = ['gid', 'srand','start', 'end', 'ptss', 'stss', 'itss', 'astss']
    # # writer = utils.csv_writer(para.TSS_DIR+id+'tss.csv', fieldnames)
    # # for gid, iseq in gene_dict.items():
    # #     writer.writerow({'gid':gid, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'ptss':iseq.ptss, 'stss':iseq.stss,
    # #                      'itss':iseq.itss, 'astss':iseq.astss})

    utils.gen_fasta(75, 95, id, '10', pdict)
    pp.meme(id, 6, '10raw')
    pdict, m10list, mextlist = pp.open_result10(gdict, pdict, id, paper_dict)
    pp.meme(id, 6, 'extmeme')
    pp.meme(id, 7, '35raw')
    pdict, m35list, spacers = pp.open_result35(gdict, pdict, id)
    pdict, m10list, mextlist, m35list, spacers = pp.iteration_base(gdict, pdict, id, 'p', m10list, mextlist, m35list, spacers)

    # # debug
    # # m10list =  utils.open_npy_list(para.MOTIF_DIR+id+'m10ilist.npy')
    # # mextlist = utils.open_npy_list(para.MOTIF_DIR+id+'mextilist.npy')
    # # m35list =  utils.open_npy_list(para.MOTIF_DIR+id+'m35ilist.npy')
    # # pdict = utils.open_npy(para.PROMOTER_DIR+id+'ppdict-iteration.npy')
    # # spacers = utils.open_npy_list(para.MOTIF_DIR+id+'spacersi.npy')

    for kind in para.klist:
        pdict = pp.pwm_tss(gdict, pdict, id, kind, m10list, mextlist, m35list, spacers)
    
    ###############################################################################
    # leaning model
    # debug
    # pdict = utils.open_npy(para.PROMOTER_DIR+id+'apdict-pwmtss.npy')

    anderson_pmt, anderson_motif, anderson_y = ip.anderson_process()

    # pdict, allmotifs = op.get_motif_list(pdict, gdict, id, ddict,sddcit, g_len)
    # pdict = op.get_score(pdict, allmotifs, id)
    # train_pmt,train_motifs,trains_scores,ly,sy,anderson_scores = op.gen_train_data(pdict, gdict, id)
    
    

    #debug
    train_pmt = utils.open_npy_list(para.TRAIN_DIR+id+'train_pmt.npy')
    train_motifs = utils.open_npy_list(para.TRAIN_DIR+id+'train_motif.npy')
    trains_scores = utils.open_npy_list(para.TRAIN_DIR+id+'trains_scores.npy')
    ly = np.load(para.TRAIN_DIR+id+'ly.npy', allow_pickle=True)
    sy = np.load(para.TRAIN_DIR+id+'sy.npy', allow_pickle=True)
    anderson_scores = utils.open_npy_list(para.TRAIN_DIR+id+'anderson_scores.npy')

    predictor = ml.PREDICT()
    
    

    fieldnames = ['model', 'data','lors','dataform','test', 'train','anderson']
    writer = utils.csv_writer(para.TRAIN_DIR + id+'model.csv', fieldnames)
    

    for i in range(3):
        writer = predictor.predict(train_pmt[i],ly[i],anderson_pmt,anderson_y,writer,i,'long','pmt')
        writer = predictor.predict(train_motifs[i],ly[i],anderson_motif,anderson_y,writer,i,'long','motif')
        writer = predictor.predict(train_pmt[i],sy[i],anderson_pmt,anderson_y,writer,i,'short','pmt')
        writer = predictor.predict(train_motifs[i],sy[i],anderson_motif,anderson_y,writer,i,'short','motif')

        # # predict only for svr
        # writer = predictor.predict(trains_scores[i],ly[i],anderson_scores[i],anderson_y,writer,i,'long','scores')
        # writer = predictor.predict(train_pmt[i],ly[i],anderson_pmt,anderson_y,writer,i,'long','pmt')
        # writer = predictor.predict(train_motifs[i],ly[i],anderson_motif,anderson_y,writer,i,'long','motif')

        # writer = predictor.predict(trains_scores[i],sy[i],anderson_scores[i],anderson_y,writer,i,'short','scores')
        # writer = predictor.predict(train_pmt[i],sy[i],anderson_pmt,anderson_y,writer,i,'short','pmt')
        # writer = predictor.predict(train_motifs[i],sy[i],anderson_motif,anderson_y,writer,i,'short','motif')

        # writer = ml.fully_connect(trains_scores[i],ly[i],anderson_scores[i],anderson_y,writer,i,'long','scores')
        # writer = ml.fully_connect(trains_scores[i],sy[i],anderson_scores[i],anderson_y,writer,i,'short','scores')
        # writer = ml.fully_connect(train_pmt[i],ly[i],anderson_pmt,anderson_y,writer,i,'long','pmt')
        # writer = ml.fully_connect(train_pmt[i],sy[i],anderson_pmt,anderson_y,writer,i,'short','pmt')
        # writer = ml.fully_connect(train_motifs[i],ly[i],anderson_motif,anderson_y,writer,i,'long','motif')
        # writer = ml.fully_connect(train_motifs[i],sy[i],anderson_motif,anderson_y,writer,i,'short','motif')

        # exp_CNN,exp_SVR = predictor.predict(np.array(train_pmt[i]),np.log(np.array(ly[i])+0.0001),
        #                     anderson_pmt,np.log(np.array(anderson_y)*max(ly[i])+0.0001))
        
        # ml.svr(trains_scores[i],np.log(np.array(ly[i])+0.0001))
        # ml.svr(trains_scores[i],np.log(np.array(sy[i])+0.0001))
        # ml.cnn_model(trains_scores[i],np.log(np.array(ly[i])+0.0001),i,'l')
        # ml.cnn_model(trains_scores[i],np.log(np.array(sy[i])+0.0001),i,'s')



        # writer = ml.fully_connect(train_pmt[i],np.log(np.array(ly[i])+0.0001),i,'lpmt',anderson_pmt,np.log(np.array(anderson_y)*max(ly[i])+0.0001))
        # writer.writerow({'model':'CNN', 'data':data_kind[i],'lors':'long','dataform':'-100~20','train':coxy[1], 'test':coxy[0],'anderson':coxy[2]})
        # writer = ml.fully_connect(train_pmt[i],np.log(np.array(sy[i])+0.0001),i,'spmt',anderson_pmt,np.log(np.array(anderson_y)*max(sy[i])+0.0001))
        # writer.writerow({'model':'CNN', 'data':data_kind[i],'lors':'short','dataform':'-100~20','train':coxy[1], 'test':coxy[0],'anderson':coxy[2]})
        # # train_motifs_line = []
        # # for item in train_motifs[i]:
        # #     train_motifs+=item
        # writer = ml.fully_connect(train_motifs[i],np.log(np.array(ly[i])+0.0001),i,'lmotif',anderson_motif,np.log(np.array(anderson_y)*max(ly[i])+0.0001))
        # writer.writerow({'model':'CNN', 'data':data_kind[i],'lors':'long','dataform':'pmt+motif','train':coxy[1], 'test':coxy[0],'anderson':coxy[2]})
        # writer = ml.fully_connect(train_motifs[i],np.log(np.array(sy[i])+0.0001),i,'smotif',anderson_motif,np.log(np.array(anderson_y)*max(sy[i])+0.0001))
        # writer.writerow({'model':'CNN', 'data':data_kind[i],'lors':'short','dataform':'pmt+motif','train':coxy[1], 'test':coxy[0],'anderson':coxy[2]})



        # ml.svr(trains_scores[i],list(map(lambda x:np.log(x+1), sy)))
        # list(map(lambda x:np.log(x+1), ly))
    # strength_dict = {} #gid:{} tss:[],kind,index,strength
    # for gid, pallitem in pdict.items():
    #     strength_dict[gid] = []
    #     sitem = strength_dict[gid]
    #     for kind, pitemlist in pallitem.items():
    #         for i in range(0, len(pitemlist)):
                
    #             pitem = pitemlist[i]
    #             tss = pitem.ctss
    #             # sitem[tss] = [kind, i]
    #             if len(sitem) == 0 or tss>sitem[-1][0]:
    #                 sitem.append([tss, kind, i]) 
    #     sitem = sorted(sitem, key=itemgetter(0))
    #     a = 0 


# cntspacers = pd.value_counts(spacersi)
#         lspacers = utils.log_ods(cntspacers)
    # a = 0
    # if not os.path.exists(para.FIG_DIR+id):
    #     os.mkdir(para.FIG_DIR+id)
    # for gid, pallitem in pdict.items():
    #     for kind, pitemlist in pallitem.items():
    #         for i in range(0, len(pitemlist)):
    #             pitem = pitemlist[i]
    #             tss = pitem.ctss
    #             inter = 200
    #             te = tss+inter
    #             srand = pitem.srand
    #             depth = sddcit[srand][tss:te]
    #             ldepth = ddict[srand][tss:te]
    #             depth_window = []
    #             for i in range(tss, tss+inter):
    #                 depth_window.append(np.mean(depth[i-tss:i+20-tss]))
    #             derivative = np.gradient(depth_window)
    #             derivative2 = np.gradient(derivative)
    #             x = np.linspace(0,inter,inter)
    #             plt.plot(x,depth_window,'blue',label = 'depthwindow')
    #             plt.plot(x,ldepth,'green',label = 'sddict')
    #             plt.plot(x,depth,'red', label = 'lddict')
    #             plt.plot(x,derivative,'orange', label = 'derivative')
    #             plt.plot(x,derivative2,'purple', label = 'derivative2')
    #             plt.legend()
                
    #             plt.savefig(para.FIG_DIR+id+'/'+gid+'.png')
    #             # plt.show()
    #             # plt.close('all')
    #     #         maxd = max(depth_plot)
    #     # maxi = np.where(depth_plot == maxd)[0][0] + ptss
    #             a = 0


























    # pltdict = {'p': 'o', 'i': 'x', 's':'*'}
    # fieldnames = ['gid','start', 'end', 'kind', 'algtss','papertss', 'diff']
    # writer = utils.csv_writer(para.TSS_DIR+id+'alg&paper.csv', fieldnames)
    # for gid, pallitem in pdict.items():
    #     if len(pallitem['p']) <= 0:
    #         continue
    #     iseq = gene_dict[gid]

    #     srand = iseq.srand
    #     start = iseq.start
    #     end = iseq.end
    #     depth = ddict[srand]
    #     sdepth = sddcit[srand]
    #     plots = max(0, start - 300)
    #     ptss = pallitem['p'][0].ctss
    #     depth_plot = depth[ptss : ptss + 200]
    #     maxd = max(depth_plot)
    #     maxi = np.where(depth_plot == maxd)[0][0] + ptss
    #     inter = maxd/10
    #     ploti = maxd
    #     textx = end - (end-start)/5
    #     # maxi = depth_plot.index(max(depth_plot)) + ptss
    #     plote = end
    #     x = np.linspace(plots,plote,plote-plots)
    #     l1 = plt.plot(x,depth[plots:plote])
    #     l2 = plt.plot(x,sdepth[plots:plote])
    #     # plt.legend([l1,l2],["long", "short"],loc='upper right')
    #     plt.axvline(start)
    #     plt.plot()
    #     for kind in para.plot_klist:
    #         for pitem in pallitem[kind]:
    #             ctss = pitem.ctss
    #             dtss = depth[ctss]
    #             plt.plot(ctss,depth[ctss],pltdict[kind])
    #             plt.text(textx, ploti, kind+'('+str(ctss)+','+str(int(dtss))+')', fontsize=10)
    #             ploti -= inter
    #             plt.plot(ctss,sdepth[ctss],pltdict[kind])
    #             plt.text(textx, ploti, 'short('+str(ctss)+','+str(int(sdepth[ctss]))+')', fontsize=10)
    #             ploti -= inter
    #             if gid in paper_dict and kind == 'p':
    #                 papertss =  paper_dict[gid]['p'][0] - 1
    #                 if srand == '-':
    #                     papertss = g_len - papertss
    #                 writer.writerow({'gid':gid, 'start':start, 'end':end, 'algtss':ctss, 'papertss':papertss, 'diff': papertss-ctss})
    #     plt.plot(maxi, maxd,'+')
    #     plt.text(maxi, maxd+15, 'max('+str(maxi)+','+str(int(maxd))+')', fontsize=10)
    #     if gid in paper_dict:
    #         pallitem = paper_dict[gid]
    #         for kind in para.plot_klist:
    #             for ctss in pallitem[kind]:
    #                 plt.plot(ctss, depth[ctss], pltdict[kind])
    #                 plt.text(textx, ploti, 'paper'+kind+'('+str(ctss)+','+str(int(depth[ctss]))+')', fontsize=10)
    #                 ploti -= inter
    #                 plt.plot(ctss,sdepth[ctss],pltdict[kind])
    #                 plt.text(textx, ploti, 'paper-short'+kind+'('+str(ctss)+','+str(int(sdepth[ctss]))+')', fontsize=10)
    #                 ploti -= inter
    #                 # writer.writerow({'gid':gid, 'start':start, 'end':end, 'kind':kind+'-paper', 'tss':ctss-1})
    #         # plt.plot(maxi,maxd,'+')
    #         # plt.text(maxi+15, maxd+15, 'max('+str(maxi)+','+str(int(maxd))+')', fontsize=10)
    #     plt.savefig(para.FIG_DIR+id+'/'+gid+'.png')
    #     # plt.show()
    #     plt.close("all")

    print('end')