import numpy as np
import csv
import pandas as pd
from Bio import motifs
from Bio.Seq import Seq
import parameters as para
import utils

def get_motif_list(pdict, gdict, id, ddict,sddcit, g_len):
    pdict2 = {}
    all_tss = {'+':[], '-':[]}
    allmotifs = {'p':[], 's':[], 'i':[]}
    # {'m10':[], 'm35':[], 'ext':[], 'spacer':[], 'itr':[], 'disc':[], 'up':[]}
    fieldnames = ['gid', 'kind', 'tss','lstr', 'sstr']
    writer = utils.csv_writer(para.TSS_DIR+id+'lsstrength.csv', fieldnames)
    for kind in para.plot_klist:
        for gid, pallitem in pdict.items():
            pitemlist = pallitem[kind]
            if gid not in pdict2:
                pdict2[gid] = {'p':[], 's':[], 'i':[]}
            for pitem in pitemlist:
                tss = pitem.ctss
                srand = pitem.srand
                if tss in all_tss[srand]:
                    continue
                all_tss[srand].append(tss)
                pitem.trainseq = gdict[srand][tss-100:tss+20]
                pitem.discseq = gdict[srand][tss-6:tss]
                pitem.spacerseq = gdict[srand][tss-12-pitem.spacer:tss-12]
                pitem.sstrength = utils.localmax(tss, sddcit[srand], g_len)
                pitem.lstrength = utils.localmax(tss, ddict[srand], g_len)
                writer.writerow({'gid':gid, 'kind':kind, 'tss':tss, 'lstr':pitem.lstrength, 'sstr':pitem.sstrength})
                allmotifs[kind].append([pitem.itrseq, pitem.discseq, pitem.m10seq, pitem.mextseq, pitem.spacerseq, pitem.m35seq, pitem.mupseq])
                pdict2[gid][kind].append(pitem)
    np.save(para.PROMOTER_DIR + id + 'pdict_after_del_repeat.npy', pdict2)
    np.save(para.MOTIF_DIR+id+'allmotifs.npy', allmotifs)
    return pdict2, allmotifs

def get_score(pdict, allmotifs, id):
    motif_kinds = np.array([])
    
    for kind in para.plot_klist:
        motif_kinds = list(motif_kinds) + allmotifs[kind]
        motif_kinds = np.array(motif_kinds)
        pwm = []
        i= 0
        for col in motif_kinds.T:
            if i == 4:
                spacers = list(map(lambda x:len(x), col))
                # tmp = spacers(col)
                pwm.append(utils.log_ods(pd.value_counts(spacers)))
            else:
                pwm.append(utils.get_pwm(col))
            i+=1

        # pwm_disc = utils.get_pwm(allmotifs)
        for gid, pallitem in pdict.items():
            pitemlist = pallitem[kind]
            for pitem in pitemlist:
                tss = pitem.ctss
                srand = pitem.srand
                pitem.itrscore = utils.get_pwm_score(pwm[0],pitem.itrseq)
                pitem.discscore = utils.get_pwm_score(pwm[1],pitem.discseq)
                pitem.m10score = utils.get_pwm_score(pwm[2],pitem.m10seq)
                pitem.mextscore = utils.get_pwm_score(pwm[3],pitem.mextseq)
                pitem.spacerscore = pwm[4][pitem.spacer]
                pitem.m35score = utils.get_pwm_score(pwm[0],pitem.m35seq)
                pitem.upscore = utils.get_pwm_score(pwm[0],pitem.mupseq)
    
    anderson_list = utils.open_npy_list(para.TRAIN_DIR+'anderson_list.npy')
    for pitem in anderson_list:
        pitem.itrscore = utils.get_pwm_score(pwm[0],pitem.itrseq)
        pitem.discscore = utils.get_pwm_score(pwm[1],pitem.discseq)
        pitem.m10score = utils.get_pwm_score(pwm[2],pitem.m10seq)
        pitem.mextscore = utils.get_pwm_score(pwm[3],pitem.mextseq)
        pitem.spacerscore = pwm[4][pitem.spacer]
        pitem.m35score = utils.get_pwm_score(pwm[0],pitem.m35seq)
        pitem.upscore = utils.get_pwm_score(pwm[0],pitem.mupseq)
    np.save(para.PROMOTER_DIR + id + 'pdict_after_score.npy', pdict)
    np.save(para.TRAIN_DIR + 'anderson_list_after_score.npy', anderson_list)
    return pdict

def gen_train_data(pdict, gdict, id):
    train_pmt = [] # all kinds training data, including n promoter seq * 3
    train_motifs = []  # all kinds trainning data, including n promoter seq and motif * 3
    trains_scores = [] # all kinds trainning data, including n promoter motif score[] * 3
    trains_scores_1035 = [] # all kinds trainning data, including n promoter motif score[] * 3
    ly = [] # all kinds long read y, including n y * 3

    anderson_scores = [] # all kinds anderson data, including n promoter motif score[] * 3
    anderson_scores_1035 = [] # all kinds anderson data, including n promoter motif score[] * 3
    
    
    sy = [] # all kinds short read y, including n y * 3
    train_pmt_s = [] # all kinds training data, including n promoter seq * 3
    train_motifs_s = []  # all kinds trainning data, including n promoter seq and motif * 3
    trains_scores_s = [] # all kinds trainning data, including n promoter motif score[] * 3
    trains_scores_1035 = [] # all kinds trainning data, including n promoter motif score[] * 3
    
    # score detail: pitem.itrseq, pitem.discseq, pitem.m10seq, pitem.mextseq, pitem.spacerseq, pitem.m35seq, pitem.mupseq
    
    i = 0
    anderson_list = utils.open_npy_list(para.TRAIN_DIR+'anderson_list_after_score.npy')
    for kind in para.plot_klist:
        kind_pmts = [] # current kind data, including n promoter seq
        kind_motifs = [] # current kind data, including n promoter and motif
        kind_ly = []
        kind_sy = []
        kind_scores = []
        kind_anderson_scores = []
        for gid, pallitem in pdict.items():
            pitemlist = pallitem[kind]
            for pitem in pitemlist: 
                tss = pitem.ctss
                srand = pitem.srand
                pitem.itrseq = gdict[srand][tss:tss+20]

                cur_pmt = utils.seq2onehot(gdict[srand][tss-100:tss+20])
                # cur_motifs = []
                # cur_motifs.append(cur_pmt)
                # cur_motifs.append(utils.seq2onehot(pitem.itrseq))
                # cur_motifs.append(utils.seq2onehot(pitem.discseq))
                # cur_motifs.append(utils.seq2onehot(pitem.m10seq))
                # cur_motifs.append(utils.seq2onehot(pitem.mextseq))
                # cur_motifs.append(utils.seq2onehot(pitem.spacerseq))
                # cur_motifs.append(utils.seq2onehot(pitem.m35seq))
                # cur_motifs.append(utils.seq2onehot(pitem.mupseq))
                cur_motifs = utils.seq2onehot(gdict[srand][tss-100:tss+20]) + para.onehot0
                cur_motifs = cur_motifs + utils.seq2onehot(pitem.itrseq) + para.onehot0
                cur_motifs = cur_motifs + utils.seq2onehot(pitem.discseq) + para.onehot0
                cur_motifs = cur_motifs + utils.seq2onehot(pitem.m10seq) + para.onehot0
                cur_motifs = cur_motifs + utils.seq2onehot(pitem.mextseq) + para.onehot0
                cur_motifs = cur_motifs + (18-pitem.spacer)*para.onehot0 + utils.seq2onehot(pitem.spacerseq) + para.onehot0
                cur_motifs = cur_motifs + utils.seq2onehot(pitem.m35seq) + para.onehot0
                cur_motifs = cur_motifs + utils.seq2onehot(pitem.mupseq)

                cur_scores = [pitem.itrscore,pitem.discscore,pitem.m10score,pitem.mextscore,pitem.spacerscore,pitem.m35score,pitem.upscore]
                
                kind_ly.append(pitem.lstrength)
                kind_sy.append(pitem.sstrength)
                kind_pmts.append(cur_pmt)
                kind_motifs.append(cur_motifs)
                kind_scores.append(cur_scores)
        for pitem in anderson_list:
            cur_scores = [pitem.itrscore,pitem.discscore,pitem.m10score,pitem.mextscore,pitem.spacerscore,pitem.m35score,pitem.upscore]
            kind_anderson_scores.append(cur_scores)
        if i > 0:
            kind_pmts = train_pmt[-1] + kind_pmts
            kind_motifs = train_motifs[-1] + kind_motifs
            kind_scores = trains_scores[-1] + kind_scores
            # kind_anderson_scores = anderson_scores[-1] + kind_anderson_scores
            kind_ly = ly[-1] + kind_ly
            kind_sy = sy[-1] + kind_sy
        train_pmt.append(kind_pmts)
        train_motifs.append(kind_motifs)
        trains_scores.append(kind_scores)
        anderson_scores.append(kind_anderson_scores)
        ly.append(kind_ly)
        sy.append(kind_sy)
        i+=1
    np.save(para.TRAIN_DIR+id+'train_pmt.npy',train_pmt)
    np.save(para.TRAIN_DIR+id+'train_motif.npy',train_motifs)
    np.save(para.TRAIN_DIR+id+'trains_scores.npy',trains_scores)
    np.save(para.TRAIN_DIR+id+'anderson_scores.npy',anderson_scores)
    np.save(para.TRAIN_DIR+id+'ly.npy',ly,allow_pickle=True)
    np.save(para.TRAIN_DIR+id+'sy.txt',sy,allow_pickle=True)
    return train_pmt,train_motifs,trains_scores,ly,sy,anderson_scores
