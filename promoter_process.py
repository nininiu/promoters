# from ctypes import util
# from readline import parse_and_bind
# from unittest.util import strclass
import numpy as np
import csv
import os
import pandas as pd
from Bio import motifs
from Bio.Seq import Seq

import parameters as para
import utils

def search_tss(start, end, depth, g_len, iseq_tss, genome, iseq, kind, pallitem, fpmt, pmtlist):
    diff_last = 0
    for i in range(max(0, start), end):
        high = np.mean(depth[i:i+10])
        low = np.mean(depth[i-10:i])
        diff = high - low
        tss_list = []
        if utils.if_start(i, depth) and high >= 1.5*low:
            if len(tss_list) > 0 and i - tss_list[-1] <= 35:
                if diff < diff_last:
                    continue
                iseq_tss.pop()
            iseq_tss, iseq, pallitem, fpmt, pmtlist = tss_assign(iseq_tss, i, genome, iseq,kind, high, low, pallitem, fpmt, pmtlist)
            diff_last = diff
            # iseq_tss.append(i)
            # apmt = genome[max(i - para.promoter_len, 0) : i]
            # pitem = para.PMTITEM(iseq,apmt,i,kind)
            # pitem.strength = (high - low)/utils.copy_num(i)
            # pallitem[kind].append(pitem)
            # if len(apmt) == para.promoter_len:
            #     fpmt.write(apmt+'\n')
            #     pmtlist.append(apmt)
    return iseq_tss, iseq, pallitem, fpmt, pmtlist

def tss_assign(iseq_tss, tss, genome, iseq,kind, high, low, pallitem, fpmt, pmtlist):
    iseq_tss.append(tss)
    pmt = genome[max(tss - para.promoter_len, 0) : tss]
    pitem = para.PMTITEM(iseq,pmt,tss,kind)
    pitem.strength = (high - low)/utils.copy_num(tss)
    pallitem[kind].append(pitem)
                        
    if len(pmt) == para.promoter_len:
        fpmt.write(pmt+'\n')
        pmtlist.append(pmt)
    return iseq_tss, iseq, pallitem, fpmt, pmtlist

def potential_promoter_cdf(ddict, gdict, gene, id, g_len):
    fp = open(para.PROMOTER_DIR + id + 'ppromoter.txt','w')
    fs = open(para.PROMOTER_DIR + id + 'spromoter.txt','w')
    fi = open(para.PROMOTER_DIR + id + 'ipromoter.txt','w')
    fa = open(para.PROMOTER_DIR + id + 'apromoter.txt','w')
    plist = []
    slist = []
    ilist = []
    alist = []
    pdict = {}
    tss_dict = {'+':[], '-':[]}
    last_end = {'+': 0, '-':0}

    # sorted(gene.items(), key = lambda kv:(kv[1], kv[0]))
    
    for gid, iseq in gene.items():
        print('processing '+gid)
        srand = iseq.srand
        start = iseq.start
        end = iseq.end
        depth = ddict[srand]

        #ptss and stss
        max_diff = 0
        last_diff = 0
        pallitem = {'p':[], 's':[], 'i':[], 'a':[]}
        pdict[gid] = pallitem
        le = last_end[srand]
        last_end[srand] = end

        for i in range(max(start - para.search_area, 0), start):
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            diff = high - low
            if utils.if_start(i, depth) and high >= 1.5 * low:
                if iseq.ptss <=  0:
                    tss_dict[srand].append(i)
                    iseq.ptss = i
                    max_diff = diff
                    last_diff = diff
                    continue

                if i in tss_dict[srand]: # if i is already identified as a tss, continue
                    continue
                    # 判断当前diff是不是最大
                    #     前一个如果是stss 丢掉
                    #     直接ptss复制
                    # 判断当前diff是不是比上一个大
                    #     stss pop 
                    #     再赋值新的stss 

                if i - tss_dict[srand][-1] < 35: # if i - last tss < 35, we think they are the same tss, choose one
                    if diff > max_diff: # if current diff is larger, we choose the current one as the true tss, discard last one
                        if tss_dict[srand][-1] != iseq.ptss and len(iseq.stss) > 0: # if last tss is stss,
                            iseq.stss.pop()
                        tss_dict[srand][-1] = i #change the last tss as current one
                        iseq.ptss = i
                        max_diff = diff
                        last_diff = diff
                        continue
                    if diff > last_diff and len(iseq.stss) > 0:
                        tss_dict[srand].pop()
                        iseq.stss.pop()
                        last_diff = diff
                        iseq.stss, iseq, pallitem, fs, slist = tss_assign(iseq.stss, i, gdict[srand], iseq,'s',
                                                                             high, low, pallitem, fs, slist)
                    continue

                if iseq.ptss > 0 and diff > max_diff:
                    iseq.stss, iseq, pallitem, fs, slist = tss_assign(iseq.stss, iseq.ptss, gdict[srand], iseq,'s', 
                                                                        high, low, pallitem, fs, slist)
                    iseq.ptss = i
                    max_diff = diff
                else:
                    iseq.stss, iseq, pallitem, fs, slist = tss_assign(iseq.stss, i, gdict[srand], iseq,'s', 
                                                                        high, low, pallitem, fs, slist)
                last_diff = diff
        if iseq.ptss > 0:
            i = iseq.ptss
            high = np.mean(depth[i:i+10])
            low = np.mean(depth[i-10:i])
            tmp, iseq, pallitem, fp, plist = tss_assign([], i, gdict[srand], iseq,'p', high, low, pallitem, fp, plist)

        #itss
        iseq.itss, iseq, pallitem, fi, ilist = search_tss(start, end, depth, g_len, iseq.itss, gdict[srand], iseq, 'i', 
                                                            pallitem, fi, ilist)
                                                            
        #astss
        asrand = para.as_dict[srand]
        depth = ddict[asrand]
        iseq.atss, iseq, pallitem, fa, alist = search_tss(max(0, g_len - end - para.search_area), g_len - end, depth, g_len, 
                                                            iseq.itss, gdict[asrand], iseq, 'a', pallitem, fa, alist)

    fp.close()
    fs.close()
    fi.close()
    fa.close()
    np.save(para.GENE_DIR + id +'gene_after_tss.npy', gene)
    np.save(para.PROMOTER_DIR + id + '-raw-pdict.npy', pdict)
    utils.plot_logo(plist, id+'p')
    utils.plot_logo(slist, id+'s')
    utils.plot_logo(ilist, id+'i')
    utils.plot_logo(alist, id+'a')

    fieldnames = ['gid', 'srand','start', 'end', 'ptss', 'stss', 'itss', 'astss']
    writer = utils.csv_writer(para.TSS_DIR+id+'tss.csv', fieldnames)
    for gid, iseq in gene.items():
        writer.writerow({'gid':gid, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'ptss':iseq.ptss, 'stss':iseq.stss,
                         'itss':iseq.itss, 'astss':iseq.astss})
    return gene, pdict

def cal_str(tss, depth):
    a = 0
    depth_window = []
    inter = 500
    for i in range(tss, tss+inter):
        depth_window.append(np.mean(depth[i:i+20]))
    derivative = np.gradient(depth_window)
    for j in range(len(derivative)):
        if derivative[j] < 0:
            break
    return depth_window[j]

def find_tss_is(start, end, depth, tslist, iseqtss, genome, kind,iseq,pallitem,fpmt,pmtlist):
    for i in range(start, end):
        if utils.if_start(i, depth):
            if i in tslist:
                continue
            tslist.append(i)
            low = np.mean(depth[max(0,i-10):i])
            cur_str = (cal_str(i, depth) - low)/utils.copy_num(i)
            if len(pallitem[kind]) > 0 and i - pallitem[kind][-1].tss < 30:
                if cur_str > pallitem[kind][-1].strength:
                    pallitem[kind].pop()
                    # tslist.pop()
                    iseqtss.pop()
                else:
                    continue
            pmt = genome[max(i - para.promoter_len, 0) : i]
            pitem = para.PMTITEM(iseq,pmt,i,kind)
            pitem.strength = cur_str
            
            iseqtss.append(i)
            
            if len(pmt) < para.promoter_len:
                a = 0
                continue
            fpmt.write(pmt+'\n')
            pmtlist.append(pmt)
            # pitem.tss = i
            # pitem.promoter = pmt
            pallitem[kind].append(pitem)
                
    return tslist, iseqtss, iseq, pallitem, fpmt, pmtlist


def promoter_cdf_psi(ddict, gdict, gene, id):
    fp = open(para.PROMOTER_DIR + id + 'ppromoter.txt','w')
    fs = open(para.PROMOTER_DIR + id + 'spromoter.txt','w')
    fi = open(para.PROMOTER_DIR + id + 'ipromoter.txt','w')
    plist = []
    slist = []
    ilist = []
    pdict = {}
    tss_dict = {'+':[], '-':[]}
    
    for gid, iseq in gene.items():
        print('processing '+gid)
        srand = iseq.srand
        start = iseq.start
        end = iseq.end
        depth = ddict[srand]
        pallitem = {'p':[para.PMTITEM()], 's':[], 'i':[]}
        pdict[gid] = pallitem

        #stss
        tss_dict[srand], iseq.stss, iseq, pallitem, fs, slist = find_tss_is(max(0,start-para.search_area), start, depth, tss_dict[srand], 
                                                                    iseq.stss, gdict[srand],'s',iseq, pallitem, fs, slist)
        #ptss
        max_str = 0
        for pitem in pallitem['s']:
            if pitem.strength > max_str:
                iseq.ptss = pitem.tss
                max_str = pitem.strength
                pallitem['p'][0] = pitem
                pmt = pitem.promoter
                if len(pmt) == para.promoter_len:
                    fp.write(pmt+'\n')
                    plist.append(pmt)
        #itss
        tss_dict[srand], iseq.itss, iseq, pallitem, fi, ilist = find_tss_is(start, end, depth, tss_dict[srand], iseq.itss, gdict[srand],
                                                                     'i',iseq, pallitem, fi, ilist)


    fp.close()
    fs.close()
    fi.close()
    np.save(para.GENE_DIR + id +'gene_after_tss.npy', gene)
    np.save(para.PROMOTER_DIR + id + '-raw-pdict.npy', pdict)
    utils.plot_logo(plist, id+'p')
    utils.plot_logo(slist, id+'s')
    utils.plot_logo(ilist, id+'i')

    fieldnames = ['gid', 'srand','start', 'end', 'ptss', 'stss', 'itss']
    writer = utils.csv_writer(para.TSS_DIR+id+'tss.csv', fieldnames)
    for gid, iseq in gene.items():
        writer.writerow({'gid':gid, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'ptss':iseq.ptss, 'stss':iseq.stss,
                         'itss':iseq.itss})
    return gene, pdict

def meme(id, w, kind):
    infile = para.FASTA_DIR+id+kind+'.fasta'
    outfile = para.RESULT_DIR+id+kind
    os.system('meme '+infile+' -dna -mod oops -nmotifs '+str(para.nmotifs)+' -w '+str(w)+' -oc '+outfile)


def open_result10(gdict, pdict, id):
    f= open(para.RESULT_DIR+id+'10raw/meme.txt','r')
    f10f= open(para.FASTA_DIR+id+'10meme.fasta','w')
    fef= open(para.FASTA_DIR+id+'extmeme.fasta','w')
    f35f = open(para.FASTA_DIR+id+'35raw.fasta','w')

    fieldnames = ['geneid','srand','tss','ctss','m10','cm10','mext']
    writer = utils.csv_writer(para.TSS_DIR+id+'10correct.csv',fieldnames)

    m10list= []
    mextlist = []
    pmt85list = []
    mitrlist = []

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
            srand = keylist[1]
            if len(pdict[key]['p']) == 1:
                pitem = pdict[key]['p'][0]
                tss = pitem.tss
                m10 = item[4]
                inter = int(item[1]) - 1
                pitem.ctss = tss - (para.e10 - para.s10 - 6 - inter)
                ctss = pitem.ctss
                pitem.m10p = ctss-12
                pitem.m10pvalue = utils.get_pvalue(item[2])
                pitem.m10seq = m10
                pitem.mextseq = gdict[srand][pitem.ctss-17:pitem.ctss-13]
                mext = pitem.mextseq
                pitem.promoter = gdict[srand][max(pitem.ctss-100, 0):pitem.ctss]
                pmt85list.append(pitem.promoter[0:85])
                pitem.itrseq = gdict[srand][pitem.ctss:pitem.ctss+20]

                m10list.append(m10)
                mextlist.append(mext)
                mitrlist.append(pitem.itrseq)
                writer.writerow({'geneid':key,'srand':srand,'tss':tss,'ctss':ctss, 
                     'm10':m10,'cm10':gdict[srand][pitem.m10p:pitem.m10p+6], 'mext':mext})
                f10f.write('>'+str(item[0])+'\n'+item[4]+'\n')
                fef.write('>'+str(item[0])+'\n'+'AA'+ mext + '\n')
                p35 = ctss-12-18-7
                m35 = gdict[srand][p35:p35+9]
                f35f.write('>'+str(item[0])+'\n'+m35+'\n')
        if s > 0:
            s+=1
            continue
    f.close()
    f10f.close()
    fef.close()
    f35f.close()
    np.save(para.PROMOTER_DIR+id+'-10meme-pdict.npy',pdict)
    np.save(para.MOTIF_DIR+id+'m10list.npy',m10list)
    np.save(para.MOTIF_DIR+id+'mextlist.npy',mextlist)
    utils.get_pssm(m10list)
    utils.get_pssm(mextlist)
    utils.get_pssm(mitrlist)
    utils.plot_logo(m10list, id+'meme-10')
    utils.plot_logo(mextlist, id+'meme-ext')
    utils.plot_logo(mitrlist, id+'meme-itr')
    utils.plot_logo(pmt85list, id+'meme-85')
    return pdict, m10list, mextlist

def open_result10_psi(gdict, pdict, id):
    f= open(para.RESULT_DIR+id+'10raw/meme.txt','r')
    f10f= open(para.FASTA_DIR+id+'10meme.fasta','w')
    fef= open(para.FASTA_DIR+id+'extmeme.fasta','w')
    f35f = open(para.FASTA_DIR+id+'35raw.fasta','w')

    fieldnames = ['geneid','srand','tss','ctss','m10','cm10','mext']
    writer = utils.csv_writer(para.TSS_DIR+id+'10correct.csv',fieldnames)

    m10list= []
    mextlist = []
    pmt85list = []
    mitrlist = []

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
            i = int(keylist[1])
            pitem = pdict[key][i]
            tss = pitem.tss
            srand = pitem.srand
            m10 = item[4]
            inter = int(item[1]) - 1
            pitem.ctss = tss - (para.e10 - para.s10 - 6 - inter)
            ctss = pitem.ctss
            pitem.m10p = ctss-12
            pitem.m10pvalue = utils.get_pvalue(item[2])
            pitem.m10seq = m10
            pitem.mextseq = gdict[srand][pitem.ctss-17:pitem.ctss-13]
            mext = pitem.mextseq
            pitem.promoter = gdict[srand][max(pitem.ctss-100, 0):pitem.ctss]
            pmt85list.append(pitem.promoter[0:85])
            pitem.itrseq = gdict[srand][pitem.ctss:pitem.ctss+20]

            m10list.append(m10)
            mextlist.append(mext)
            mitrlist.append(pitem.itrseq)
            writer.writerow({'geneid':key,'srand':srand,'tss':tss,'ctss':ctss, 
                            'm10':m10,'cm10':gdict[srand][pitem.m10p:pitem.m10p+6], 'mext':mext})
            f10f.write('>'+str(item[0])+'\n'+item[4]+'\n')
            fef.write('>'+str(item[0])+'\n'+'AA'+ mext + '\n')
            p35 = ctss-12-18-7
            m35 = gdict[srand][p35:p35+9]
            f35f.write('>'+str(item[0])+'\n'+m35+'\n')
        if s > 0:
            s+=1
            continue
    f.close()
    f10f.close()
    fef.close()
    f35f.close()
    np.save(para.PROMOTER_DIR+id+'-10meme-pdict.npy',pdict)
    np.save(para.MOTIF_DIR+id+'m10list.npy',m10list)
    np.save(para.MOTIF_DIR+id+'mextlist.npy',mextlist)
    utils.get_pssm(m10list)
    utils.get_pssm(mextlist)
    utils.get_pssm(mitrlist)
    utils.plot_logo(m10list, id+'meme-10')
    utils.plot_logo(mextlist, id+'meme-ext')
    utils.plot_logo(mitrlist, id+'meme-itr')
    utils.plot_logo(pmt85list, id+'meme-85')
    return pdict, m10list, mextlist

def compare_paper(gdict, pdict, id, paper_dict):
    f= open(para.RESULT_DIR+id+'10raw/meme.txt','r')
    f10f= open(para.FASTA_DIR+id+'10meme.fasta','w')
    fef= open(para.FASTA_DIR+id+'extmeme.fasta','w')
    f35f = open(para.FASTA_DIR+id+'35raw.fasta','w')

    fieldnames = ['geneid','srand','tss','ctss','paper_tss','m10','cm10','paper_m10','mext']
    writer = utils.csv_writer(para.TSS_DIR+id+'10correct.csv',fieldnames)

    m10list= []
    mextlist = []
    pmt85list = []
    mitrlist = []

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
            srand = keylist[1]
            if len(pdict[key]['p']) == 1:
                pitem = pdict[key]['p'][0]
                tss = pitem.tss
                m10 = item[4]
                inter = int(item[1]) - 1
                pitem.ctss = tss - (para.e10 - para.s10 - 6 - inter)
                ctss = pitem.ctss
                pitem.m10p = ctss-12
                pitem.m10pvalue = utils.get_pvalue(item[2])
                pitem.m10seq = m10
                pitem.mextseq = gdict[srand][pitem.ctss-17:pitem.ctss-13]
                mext = pitem.mextseq
                pitem.promoter = gdict[srand][max(pitem.ctss-100, 0):pitem.ctss]
                pmt85list.append(pitem.promoter[0:85])
                pitem.itrseq = gdict[srand][pitem.ctss:pitem.ctss+20]

                m10list.append(m10)
                mextlist.append(mext)
                mitrlist.append(pitem.itrseq)
                if key in paper_dict:
                    paper_tss = paper_dict[key]['p'][0]
                    paper_m10 = gdict[srand][paper_tss-12:paper_tss-6]
                writer.writerow({'geneid':key,'srand':srand,'tss':tss,'ctss':ctss, 'paper_tss':paper_tss,
                     'm10':m10,'cm10':gdict[srand][pitem.m10p:pitem.m10p+6], 'paper_m10':paper_m10, 'mext':mext})
                f10f.write('>'+str(item[0])+'\n'+item[4]+'\n')
                fef.write('>'+str(item[0])+'\n'+'AA'+ mext + '\n')
                p35 = ctss-12-18-7
                m35 = gdict[srand][p35:p35+9]
                f35f.write('>'+str(item[0])+'\n'+m35+'\n')
        if s > 0:
            s+=1
            continue
    f.close()
    f10f.close()
    fef.close()
    f35f.close()
    np.save(para.PROMOTER_DIR+id+'-10meme-pdict.npy',pdict)
    np.save(para.MOTIF_DIR+id+'m10list.npy',m10list)
    np.save(para.MOTIF_DIR+id+'mextlist.npy',mextlist)
    utils.get_pssm(m10list)
    utils.get_pssm(mextlist)
    utils.get_pssm(mitrlist)
    utils.plot_logo(m10list, id+'meme-10')
    utils.plot_logo(mextlist, id+'meme-ext')
    utils.plot_logo(mitrlist, id+'meme-itr')
    utils.plot_logo(pmt85list, id+'meme-85')
    return pdict, m10list, mextlist

def open_result35(gdict, pdict, id):
    fext= open(para.RESULT_DIR+id+'extmeme/meme.txt','r')
    f35= open(para.RESULT_DIR+id+'35raw/meme.txt','r')
    f35f = open(para.FASTA_DIR+id+'35meme.fasta','w')

    fieldnames = ['geneid','srand','spacer','m35','cm35']
    writer = utils.csv_writer(para.TSS_DIR+id+'35correct.csv',fieldnames)

    m35list = []
    spacers = []
    muplist = []

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
            if len(pdict[key]['p']) == 1: 
                pdict[key]['p'][0].mextpvalue = utils.get_pvalue(item[2])
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
            # print(item)
            spacer = 19-int(item[1])
            m35 = item[4]
            if len(pdict[key]['p']) == 1:
                pitem = pdict[key]['p'][0]
                ctss = pitem.ctss
                upp=ctss-12-spacer-7-20
                m35p=ctss-12-spacer-7
                up = gdict[srand][upp:upp+20]

                pitem.m35p = m35p
                pitem.m35seq = m35
                pitem.spacer = spacer
                pitem.mupseq = up
                pitem.m35pvalue = utils.get_pvalue(item[2])

                m35list.append(m35)
                spacers.append(spacer)
                muplist.append(up)
                writer.writerow({'geneid':key,'srand':srand,'spacer':spacer,'m35':m35,
                    'cm35':gdict[srand][pitem.m35p:pitem.m35p+7]})
                f35f.write('>'+item[0]+'\n'+m35+'\n')
        if s > 0:
            s+=1
            continue
    np.save(para.PROMOTER_DIR+id+'pdict_35meme.npy',pdict)
    np.save(para.MOTIF_DIR+id+'m35list.npy',m35list)
    np.save(para.MOTIF_DIR+id+'spacers.npy',spacers)
    np.save(para.MOTIF_DIR+id+'muplist.npy',muplist)
    result = pd.value_counts(spacers)
    print(result)
    sm=np.mean(spacers)
    utils.get_pssm(m35list)
    utils.plot_logo(m35list, id+'meme35')
    utils.plot_logo(muplist, id+'memeup')
    fext.close()
    f35.close()
    f35f.close()
    return pdict, m35list, spacers

def open_result35_psi(gdict, pdict, id):
    fext= open(para.RESULT_DIR+id+'extmeme/meme.txt','r')
    f35= open(para.RESULT_DIR+id+'35raw/meme.txt','r')
    f35f = open(para.FASTA_DIR+id+'35meme.fasta','w')

    fieldnames = ['geneid','srand','spacer','m35','cm35']
    writer = utils.csv_writer(para.TSS_DIR+id+'35correct.csv',fieldnames)

    m35list = []
    spacers = []
    muplist = []

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
            i = int(keylist[1])
            pdict[key][i].mextpvalue = utils.get_pvalue(item[2])
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
            i = int(keylist[1])
            spacer = 19-int(item[1])
            m35 = item[4]
            pitem = pdict[key][i]
            srand = pitem.srand
            ctss = pitem.ctss
            upp=ctss-12-spacer-7-20
            m35p=ctss-12-spacer-7
            up = gdict[srand][upp:upp+20]

            pitem.m35p = m35p
            pitem.m35seq = m35
            pitem.spacer = spacer
            pitem.mupseq = up
            pitem.m35pvalue = utils.get_pvalue(item[2])

            m35list.append(m35)
            spacers.append(spacer)
            muplist.append(up)
            writer.writerow({'geneid':key,'srand':srand,'spacer':spacer,'m35':m35,
                    'cm35':gdict[srand][pitem.m35p:pitem.m35p+7]})
            f35f.write('>'+item[0]+'\n'+m35+'\n')
        if s > 0:
            s+=1
            continue
    np.save(para.PROMOTER_DIR+id+'pdict_35meme.npy',pdict)
    np.save(para.MOTIF_DIR+id+'m35list.npy',m35list)
    np.save(para.MOTIF_DIR+id+'spacers.npy',spacers)
    np.save(para.MOTIF_DIR+id+'muplist.npy',muplist)
    result = pd.value_counts(spacers)
    print(result)
    sm=np.mean(spacers)
    utils.get_pssm(m35list)
    utils.plot_logo(m35list, id+'meme35')
    utils.plot_logo(muplist, id+'memeup')
    fext.close()
    f35.close()
    f35f.close()
    return pdict, m35list, spacers

def iteration_base(gdict, pdict, id, kind, m10list, mextlist, m35list, spacers):
    id = id + kind
    pssm10 = utils.get_pssm(m10list)
    pssmext = utils.get_pssm(mextlist)
    pssm35 = utils.get_pssm(m35list)
    cntspacers = pd.value_counts(spacers)
        
    lspacers = utils.log_ods(cntspacers)
    print(cntspacers)
    print(lspacers)

    cnt = 0
    while 1:
        fieldnames = ['gid','tss','ctss','10','ext','spacer','35','up','itr','m10score','mextscore',
                        'm35score','spacerscore','promoter']
        writer = utils.csv_writer(para.PROMOTER_DIR+id+'-allpromoters.csv',fieldnames)

        m10ilist = []
        mextilist = []
        m35ilist = []
        spacersi = []
        muplist = []
        mitrlist = []
        p85list = []
        pcnt = 0
        for gid, pallitem in pdict.items():
            if gid == 'b3813':
                a = 0
            for pitem in pallitem[kind]:
                if pitem.ctss < 55:
                    continue
                srand = pitem.srand
                pcnt+=1
                mscore = para.INF_MIN
                for i in range(100-para.e10,100-para.s10):
                    ctss = pitem.tss-i+6
                    if ctss < 55:
                        continue
                    m10 = gdict[srand][ctss-12:ctss-6]
                    mext = gdict[srand][ctss-17:ctss-13]
                    m10score = pssm10.calculate(Seq(m10))
                    mextscore = pssmext.calculate(Seq(mext))
                    m35sscore = para.INF_MIN
                    tm35=''
                    tspacer = 0
                    tm35score=para.INF_MIN
                    tsscore=para.INF_MIN
                    for j in range(16,19):
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
                    if ctotal>mscore: #and (m10score>m35score or pitem.correcttss ==0):
                        ti = i
                        mscore = ctotal
                        up = gdict[srand][ctss-12-pitem.spacer-7-20:ctss-12-pitem.spacer-7]
                        itr = gdict[srand][ctss:ctss+20]
                        pmt = gdict[srand][ctss-100:ctss]

                        pitem.m10seq = m10
                        pitem.mextseq = mext
                        pitem.ctss = ctss
                        pitem.m10p = ctss-12
                        pitem.mextp = ctss-17
                        pitem.m10score = m10score
                        pitem.mextscore = mextscore
                        pitem.m35seq = tm35
                        pitem.spacer = tspacer
                        pitem.m35score = tm35score
                        pitem.mspacerscore = tsscore
                        pitem.mupseq = up
                        pitem.mitrseq = itr
                        pitem.promoter = pmt

                        m10 = pitem.m10seq
                        m35 = pitem.m35seq
                        mspacer = gdict[srand][ctss-12-pitem.spacer:ctss-12]
                        m35s10 = m35+mspacer+m10
                        target = gdict[srand][ctss-12-pitem.spacer-7:ctss-6]
                        if m35s10 != target:
                            a=0
                    
        f10f = open(para.FASTA_DIR+id+'10i.fasta','w')
        fef = open(para.FASTA_DIR+id+'exti.fasta','w')
        f35f = open(para.FASTA_DIR+id+'35i.fasta','w')

        pcnt = 0
        for gid, pallitem in pdict.items():
            for pitem in pallitem[kind]:
                ctss = pitem.ctss
                if ctss < 55 or pitem.spacer < 15:
                    continue
                # up = gdict[srand][ctss-12-pitem.spacer-7-20:ctss-12-pitem.spacer-7]
                f10f.write('>'+gid+'\n'+pitem.m10seq+'\n')
                fef.write('>'+gid+'\n'+'AA'+pitem.mextseq+'\n')
                f35f.write('>'+gid+'\n'+pitem.m35seq+'\n')

                m10ilist.append(pitem.m10seq)
                mextilist.append(pitem.mextseq)
                m35ilist.append(pitem.m35seq)
                spacersi.append(pitem.spacer)
                mitrlist.append(pitem.itrseq)
                muplist.append(pitem.mupseq)
                p85list.append(pitem.promoter[0:85])

                m10 = pitem.m10seq
                m35 = pitem.m35seq
                mspacer = gdict[srand][ctss-12-pitem.spacer:ctss-12]
                m35s10 = m35+mspacer+m10
                target = gdict[srand][ctss-12-pitem.spacer-7:ctss-6]
                if m35s10 != target:
                    a=0
                
                writer.writerow({'gid':gid,'tss':pitem.tss,'ctss':ctss,'10':pitem.m10seq,'ext':pitem.mextseq,'spacer':gdict[srand][ctss-12-pitem.spacer:ctss-12],
                    '35':pitem.m35seq,'up':pitem.mupseq,'m10score':pitem.m10score,'mextscore':pitem.mextscore,'m35score':pitem.m35score,
                    'spacerscore':pitem.spacerscore,'promoter':gdict[srand][ctss-100:ctss]})
                pcnt+=1
        utils.plot_logo(m10ilist,id+'-i10-'+str(cnt))
        utils.plot_logo(mextilist,id+'-iext-'+str(cnt))
        utils.plot_logo(m35ilist,id+'-i35-'+str(cnt))
        utils.plot_logo(muplist,id+'-iup-'+str(cnt))
        utils.plot_logo(mitrlist,id+'-iitr-'+str(cnt))
        utils.plot_logo(p85list,id+'-i85-'+str(cnt))

        pssm10 = utils.get_pssm(m10ilist)
        pssmext = utils.get_pssm(mextilist)
        pssm35 = utils.get_pssm(m35ilist)
        print(np.mean(spacers))
        f10f.close()
        fef.close()
        f35f.close()
        if utils.compare(m10ilist, m10list) and utils.compare(mextilist, mextlist) and utils.compare(m35ilist, m35list):
            np.save(para.MOTIF_DIR+id+'m10ilist.npy',m10ilist)
            np.save(para.MOTIF_DIR+id+'mextilist.npy',mextilist)
            np.save(para.MOTIF_DIR+id+'m35ilist.npy',m35ilist)
            np.save(para.PROMOTER_DIR+id+'pdict-iteration.npy',pdict)
            np.save(para.MOTIF_DIR+id+'spacersi.npy',spacersi)
            break
        m10list = m10ilist
        mextlist = mextilist
        m35list = m35ilist
        cnt+=1
        cntspacers = pd.value_counts(spacersi)
        lspacers = utils.log_ods(cntspacers)
        print(cntspacers)
        print(lspacers)
    return pdict, m10ilist, mextilist, m35ilist, spacersi

def iteration_base_psi(gdict, pdict, id, m10list, mextlist, m35list, spacers):
    pssm10 = utils.get_pssm(m10list)
    pssmext = utils.get_pssm(mextlist)
    pssm35 = utils.get_pssm(m35list)
    cntspacers = pd.value_counts(spacers)
        
    lspacers = utils.log_ods(cntspacers)
    print(cntspacers)
    print(lspacers)

    cnt = 0
    while 1:
        fieldnames = ['gid','tss','ctss','10','ext','spacer','35','up','itr','m10score','mextscore',
                        'm35score','spacerscore','promoter']
        writer = utils.csv_writer(para.PROMOTER_DIR+id+'-allpromoters.csv',fieldnames)

        m10ilist = []
        mextilist = []
        m35ilist = []
        spacersi = []
        muplist = []
        mitrlist = []
        p85list = []
        pcnt = 0
        for gid, pitemlist in pdict.items():
            for pitem in pitemlist:
                if pitem.ctss <= 55:
                    continue
                srand = pitem.srand
                pcnt+=1
                mscore = para.INF_MIN
                for i in range(100-para.e10,100-para.s10):
                    ctss = pitem.tss-i+6
                    if ctss == 0:
                        a = 0
                    m10 = gdict[srand][ctss-12:ctss-6]
                    mext = gdict[srand][ctss-17:ctss-13]
                    m10score = pssm10.calculate(Seq(m10))
                    mextscore = pssmext.calculate(Seq(mext))
                    m35sscore = para.INF_MIN
                    tm35=''
                    tspacer = 0
                    tm35score=para.INF_MIN
                    tsscore=para.INF_MIN
                    for j in range(16,19):
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
                    if ctotal>mscore: #and (m10score>m35score or pitem.correcttss ==0):
                        ti = i
                        mscore = ctotal
                        up = gdict[srand][ctss-12-pitem.spacer-7-20:ctss-12-pitem.spacer-7]
                        itr = gdict[srand][ctss:ctss+20]
                        pmt = gdict[srand][ctss-100:ctss]

                        pitem.m10seq = m10
                        pitem.mextseq = mext
                        pitem.ctss = ctss
                        pitem.m10p = ctss-12
                        pitem.mextp = ctss-17
                        pitem.m10score = m10score
                        pitem.mextscore = mextscore
                        pitem.m35seq = tm35
                        pitem.spacer = tspacer
                        pitem.m35score = tm35score
                        pitem.mspacerscore = tsscore
                        pitem.mupseq = up
                        pitem.itrseq = itr
                        pitem.promoter = pmt

                        m10 = pitem.m10seq
                        m35 = pitem.m35seq
                        mspacer = gdict[srand][ctss-12-pitem.spacer:ctss-12]
                        m35s10 = m35+mspacer+m10
                        target = gdict[srand][ctss-12-pitem.spacer-7:ctss-6]
                        if m35s10 != target:
                            a=0
                    
        f10f = open(para.FASTA_DIR+id+'10i.fasta','w')
        fef = open(para.FASTA_DIR+id+'exti.fasta','w')
        f35f = open(para.FASTA_DIR+id+'35i.fasta','w')

        pcnt = 0
        for gid, pitemlist in pdict.items():
            for pitem in pitemlist:
                if pitem.tss <=0:
                    continue
                ctss = pitem.ctss
                # up = gdict[srand][ctss-12-pitem.spacer-7-20:ctss-12-pitem.spacer-7]
                f10f.write('>'+gid+'\n'+pitem.m10seq+'\n')
                fef.write('>'+gid+'\n'+'AA'+pitem.mextseq+'\n')
                f35f.write('>'+gid+'\n'+pitem.m35seq+'\n')

                m10ilist.append(pitem.m10seq)
                mextilist.append(pitem.mextseq)
                m35ilist.append(pitem.m35seq)
                spacersi.append(pitem.spacer)
                mitrlist.append(pitem.itrseq)
                muplist.append(pitem.mupseq)
                p85list.append(pitem.promoter[0:85])

                m10 = pitem.m10seq
                m35 = pitem.m35seq
                mspacer = gdict[srand][ctss-12-pitem.spacer:ctss-12]
                m35s10 = m35+mspacer+m10
                target = gdict[srand][ctss-12-pitem.spacer-7:ctss-6]
                if m35s10 != target:
                    a=0
                
                writer.writerow({'gid':gid,'tss':pitem.tss,'ctss':ctss,'10':pitem.m10seq,'ext':pitem.mextseq,'spacer':gdict[srand][ctss-12-pitem.spacer:ctss-12],
                    '35':pitem.m35seq,'up':pitem.mupseq,'m10score':pitem.m10score,'mextscore':pitem.mextscore,'m35score':pitem.m35score,
                    'spacerscore':pitem.mspacerscore,'promoter':gdict[srand][ctss-100:ctss]})
                pcnt+=1
        utils.plot_logo(m10ilist,id+'-i10-'+str(cnt))
        utils.plot_logo(mextilist,id+'-iext-'+str(cnt))
        utils.plot_logo(m35ilist,id+'-i35-'+str(cnt))
        utils.plot_logo(muplist,id+'-iup-'+str(cnt))
        utils.plot_logo(mitrlist,id+'-iitr-'+str(cnt))
        utils.plot_logo(p85list,id+'-i85-'+str(cnt))

        pssm10 = utils.get_pssm(m10ilist)
        pssmext = utils.get_pssm(mextilist)
        pssm35 = utils.get_pssm(m35ilist)
        print(np.mean(spacers))
        f10f.close()
        fef.close()
        f35f.close()
        if utils.compare(m10ilist, m10list) and utils.compare(mextilist, mextlist) and utils.compare(m35ilist, m35list):
            np.save(para.MOTIF_DIR+id+'m10ilist.npy',m10ilist)
            np.save(para.MOTIF_DIR+id+'mextilist.npy',mextilist)
            np.save(para.MOTIF_DIR+id+'m35ilist.npy',m35ilist)
            np.save(para.PROMOTER_DIR+id+'pdict-iteration.npy',pdict)
            np.save(para.MOTIF_DIR+id+'spacersi.npy',spacersi)
            break
        m10list = m10ilist
        mextlist = mextilist
        m35list = m35ilist
        cnt+=1
        cntspacers = pd.value_counts(spacersi)
        lspacers = utils.log_ods(cntspacers)
        print(cntspacers)
        print(lspacers)
    return pdict, m10ilist, mextilist, m35ilist, spacersi

def pwm_tss(gdict, pdict, id, kind, m10list, mextlist, m35list, spacers):
    id = id + kind
    pssm10 = utils.get_pssm(m10list)
    pssmext = utils.get_pssm(mextlist)
    pssm35 = utils.get_pssm(m35list)
    cntspacers = pd.value_counts(spacers)
        
    lspacers = utils.log_ods(cntspacers)
    print(cntspacers)
    print(lspacers)

    fieldnames = ['gid','tss','ctss','10','ext','spacer','35','up','itr','m10score','mextscore',
                        'm35score','spacerscore','promoter']
    writer = utils.csv_writer(para.PROMOTER_DIR+id+'-allpromoters.csv',fieldnames)

    m10ilist = []
    mextilist = []
    m35ilist = []
    spacersi = []
    muplist = []
    mitrlist = []
    p85list = []
    pcnt = 0
    for gid, pallitem in pdict.items():
        if len( pallitem[kind]) == 0:
            continue
        for pitem in pallitem[kind]:
            if pitem.tss <= 55 or pitem.ctss <= 55:
                continue
            srand = pitem.srand
            pcnt+=1
            mscore = para.INF_MIN
            for i in range(5,20):
                ctss = pitem.tss-i+6
                m10 = gdict[srand][ctss-12:ctss-6]
                mext = gdict[srand][ctss-17:ctss-13]
                m10score = pssm10.calculate(Seq(m10))
                mextscore = pssmext.calculate(Seq(mext))
                m35sscore = para.INF_MIN
                tm35=''
                tspacer = 0
                tm35score=para.INF_MIN
                tsscore=para.INF_MIN
                for j in range(16,19):
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
                if ctotal>mscore: #and (m10score>m35score or pitem.correcttss ==0):
                    mscore = ctotal
                    up = gdict[srand][ctss-12-pitem.spacer-7-20:ctss-12-pitem.spacer-7]
                    itr = gdict[srand][ctss:ctss+20]
                    pmt = gdict[srand][ctss-100:ctss]

                    pitem.m10seq = m10
                    pitem.mextseq = mext
                    pitem.ctss = ctss
                    pitem.m10p = ctss-12
                    pitem.mextp = ctss-17
                    pitem.m10score = m10score
                    pitem.mextscore = mextscore
                    pitem.m35seq = tm35
                    pitem.spacer = tspacer
                    pitem.m35score = tm35score
                    pitem.mspacerscore = tsscore
                    pitem.mupseq = up
                    pitem.itrseq = itr
                    pitem.promoter = pmt

                    m10 = pitem.m10seq
                    m35 = pitem.m35seq
                    mspacer = gdict[srand][ctss-12-pitem.spacer:ctss-12]
                
    f10f = open(para.FASTA_DIR+id+'10i.fasta','w')
    fef = open(para.FASTA_DIR+id+'exti.fasta','w')
    f35f = open(para.FASTA_DIR+id+'35i.fasta','w')

    pcnt = 0
    for gid, pallitem in pdict.items():
        if len( pallitem[kind]) == 0:
            continue
        for pitem in pallitem[kind]:
            ctss = pitem.ctss
            srand = pitem.srand
            f10f.write('>'+gid+'\n'+pitem.m10seq+'\n')
            fef.write('>'+gid+'\n'+pitem.mextseq+'\n')
            f35f.write('>'+gid+'\n'+pitem.m35seq+'\n')

            m10ilist.append(pitem.m10seq)
            mextilist.append(pitem.mextseq)
            m35ilist.append(pitem.m35seq)
            spacersi.append(pitem.spacer)
            mitrlist.append(pitem.itrseq)
            muplist.append(pitem.mupseq)
            p85list.append(pitem.promoter[0:85])

            m10 = pitem.m10seq
            m35 = pitem.m35seq

            writer.writerow({'gid':gid,'tss':pitem.tss,'ctss':ctss,'10':pitem.m10seq,'ext':pitem.mextseq,'spacer':gdict[srand][ctss-12-pitem.spacer:ctss-12],
                '35':pitem.m35seq,'up':pitem.mupseq,'m10score':pitem.m10score,'mextscore':pitem.mextscore,'m35score':pitem.m35score,
                'spacerscore':pitem.spacerscore,'promoter':gdict[srand][ctss-100:ctss]})
            pcnt+=1
    utils.plot_logo(m10ilist,id+'-i10')
    utils.plot_logo(mextilist,id+'-iext')
    utils.plot_logo(m35ilist,id+'-i35')
    utils.plot_logo(muplist,id+'-iup')
    utils.plot_logo(mitrlist,id+'-iitr')
    utils.plot_logo(p85list,id+'-i85')

    pssm10 = utils.get_pssm(m10ilist)
    pssmext = utils.get_pssm(mextilist)
    pssm35 = utils.get_pssm(m35ilist)
    print(np.mean(spacers))
    f10f.close()
    fef.close()
    f35f.close()
    np.save(para.PROMOTER_DIR+id+'pdict-pwmtss.npy',pdict)
    return pdict