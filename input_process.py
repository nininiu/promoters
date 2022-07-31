import numpy as np
import csv
from Bio import motifs
from Bio.Seq import Seq
import parameters as para
import utils

def gen_gene_dict(id, gdict):
    fpkm_dict, tpm_dict, cov_dict = utils.parse_fpkm(id)
    g_len = len(gdict['+'])
    gene_dict = {}
    locus_dict = {}

    fieldnames = ['locus_tag', 'srand','start', 'end', 'len', 'start_codon', 'stop_codon', 'fpkm', 'tpm', 'cov',
                     'geneid', 'locus', 'protein_product', 'seq']
    writer = utils.csv_writer(para.GENE_DIR + id+'.csv', fieldnames)

    csv_reader = csv.reader(open(para.GTF_DIR+para.protein_csv_file))
    csv_reader.__next__()
    for line in csv_reader:
        iseq = para.GENEITEM()
        start = int(line[2])
        end = int(line[3])
        iseq.srand = line[4]
        iseq.geneid = line[5]
        iseq.locus = line[6]
        iseq.locus_tag = line[7]
        gid = iseq.locus_tag
        iseq.protein_product = line[8]
        istart = iseq.get_start(start, end, g_len)
        iend = iseq.get_end(start, end, g_len)
        iseq.set_seq(gdict, istart, iend)
        iseq.start_codon = iseq.seq[0:3]
        iseq.stop_codon = iseq.seq[-3:]
        iseq.set_len(start, end)
        if iseq.start_codon not in para.protein_start or iseq.stop_codon not in para.protein_end:
            a = 0
        if gid in fpkm_dict:
            iseq.fpkm = fpkm_dict[gid]
            iseq.cov = cov_dict[gid]
            iseq.tpm = tpm_dict[gid]
        writer.writerow({'locus_tag':iseq.locus_tag, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'len':iseq.len,
                     'start_codon':iseq.start_codon, 'stop_codon':iseq.stop_codon, 'fpkm':iseq.fpkm, 'tpm':iseq.tpm, 
                     'cov':iseq.cov, 'geneid':iseq.geneid, 'locus':iseq.locus, 'protein_product':iseq.protein_product, 'seq':iseq.seq})
        gene_dict[gid] = iseq
        locus_dict[iseq.locus] = gid
    np.save(para.GENE_DIR+para.bid+id+'csv_gene_dict.npy', gene_dict)
    np.save(para.GENE_DIR+para.bid+'locus_dict.npy', locus_dict)
    return gene_dict

def gen_ddict(id, g_len, suffix = para.suffix, inter=para.dp_para):
    depthf = open_depth(para.DEPTH_DIR + id +'_f.'+suffix, g_len, inter=inter)
    depthr = open_depth(para.DEPTH_DIR + id +'_r.'+suffix, g_len, inter=inter)
    depthr = depthr[::-1]
    ddict = {'+':depthf,'-':depthr}
    np.save(para.DEPTH_DIR + id +'ddict.npy',ddict)
    return ddict

def open_depth(input_file, g_len, inter = para.dp_para):
    dth = np.zeros((g_len))
    f= open(input_file,'r')
    for item in f:
        item = item.split()
        dth[int(item[1-inter])-1]=float(item[2-inter])
    f.close()
    d_mean = np.mean(dth)
    dth = dth/d_mean * 20
    return dth

def open_depth_wig(input_file, g_len):
    dth = np.zeros((g_len))
    f= open(input_file,'r')
    for item in f:
        item = item.split()
        dth[int(item[0])-1]=float(item[1])
    f.close()
    return dth

def gen_paper_dict(gene):
    paper_dict = {}
    csv_reader = csv.reader(open(para.paper_tss_csv))
    csv_reader.__next__()
    for line in csv_reader:
        gid = line[0]
        srand = line[1]
        if gid in gene and srand == gene[gid].srand:
            pallitem = {'p':[], 's':[], 'i':[], 'a':[]}
            pallitem['p'] = [int(line[2])]
            pallitem['s'] = utils.str2list(line[3])
            pallitem['i'] = utils.str2list(line[4])
            pallitem['a'] = utils.str2list(line[5])
            paper_dict[gid] = pallitem
    np.save(para.GENE_DIR+'paper_dict.npy', paper_dict)
    return paper_dict

def open_genome():
    genome = ''
    fg = para.GTF_DIR+para.genome_file
    # './00_genome/Ecoli_K12_MG1655_NC_000913.2.fasta'
    f= open(fg,'r')
    i = 0
    for item in f:
        if i > 0:
            genome = genome + item[0:-1]
        i = i + 1
    f.close()
    print('genome len:', len(genome), end='\n')
    genome_seq = Seq(genome)
    return genome_seq

def gen_gdict():
    genome = open_genome() #Positive strand
    genomer = genome.reverse_complement()  #negative strand
    gdict = {'+':str(genome),'-':str(genomer)} #whole genome dict
    np.save(para.GTF_DIR+para.bid+'_gdict.npy',gdict)
    return gdict

def get_gdict():
    return utils.open_npy(para.GTF_DIR+para.bid+'_gdict.npy')

def parse_gff(id, gdict):
    fpkm_dict, tpm_dict, cov_dict = utils.parse_fpkm(id)
    g_len = len(gdict['+'])
    gene_dict = {}

    fieldnames = ['locus_tag', 'srand','start', 'end', 'len', 'start_codon', 'stop_codon', 'fpkm', 'tpm', 'cov',
                     'geneid', 'protein_product', 'seq']
    writer = utils.csv_writer(para.GENE_DIR + id +'.csv', fieldnames)

    data = np.loadtxt(para.MG1655_9132_gtf_file, dtype=np.str, delimiter='\t').tolist()
    for item in data:
        if len(item) > 8 and item[2] == 'CDS':
            # feature= item[2]
            iseq = para.GENEITEM()
            start = int(item[3])
            end = int(item[4])
            iseq.srand = item[6]
            att_list = item[8].split(';')
            for att in att_list:
                if 'locus_tag' in att:
                    iseq.locus_tag = att.split('=')[1]
                elif 'Name=' in att:
                    iseq.protein_product = att.split('=')[1]
                elif 'gene=' in att:
                    iseq.geneid = att.split('=')[1]
            
            # start = int(line[2])
            # end = int(line[3])
            # iseq.srand = line[4]
            # iseq.geneid = line[5]
            # iseq.locus = line[6]
            # iseq.locus_tag = line[7]
            gid = iseq.locus_tag
            # iseq.protein_product = line[8]
            istart = iseq.get_start(start, end, g_len)
            iend = iseq.get_end(start, end, g_len)
            iseq.set_seq(gdict, istart, iend)
            iseq.start_codon = iseq.seq[0:3]
            iseq.stop_codon = iseq.seq[-3:]
            iseq.set_len(start, end)
            if iseq.start_codon not in para.protein_start or iseq.stop_codon not in para.protein_end:
                continue
            if gid in fpkm_dict:
                iseq.fpkm = fpkm_dict[gid]
                iseq.cov = cov_dict[gid]
                iseq.tpm = tpm_dict[gid]
            writer.writerow({'locus_tag':iseq.locus_tag, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'len':iseq.len,
                        'start_codon':iseq.start_codon, 'stop_codon':iseq.stop_codon, 'fpkm':iseq.fpkm, 'tpm':iseq.tpm, 
                        'cov':iseq.cov, 'geneid':iseq.geneid, 'protein_product':iseq.protein_product, 'seq':iseq.seq})
            gene_dict[gid] = iseq
    np.save(para.GENE_DIR+para.bid+id+'csv_gene_dict.npy', gene_dict)
    return gene_dict

def srand_gff(id, gdict, gene_dict):
    g_len = len(gdict['+'])
    fieldnames = ['locus_tag', 'srand','start', 'end', 'len', 'start_codon', 'stop_codon', 'fpkm', 'tpm', 'cov',
                     'geneid', 'protein_product', 'seq']
    writer = utils.csv_writer(para.GENE_DIR + id +'.csv', fieldnames)

    data = np.loadtxt(para.GTF_DIR+para.gff_file, dtype=np.str, delimiter='\t').tolist()
    for item in data:
        if len(item) > 8 and item[2] == 'gene':
            att_list = item[8].split(';')
            locus_tag = ''
            for att in att_list:
                if 'locus_tag' in att:
                    locus_tag = att.split('=')[1]
                    break
            if locus_tag not in gene_dict:
                continue
            iseq = gene_dict[locus_tag]
            start = int(item[3])
            end = int(item[4])
            iseq.srand = item[6]
            istart = iseq.get_start(start, end, g_len)
            iend = iseq.get_end(start, end, g_len)
            iseq.set_seq(gdict, istart, iend)
            iseq.start_codon = iseq.seq[0:3]
            iseq.stop_codon = iseq.seq[-3:]
            iseq.set_len(start, end)
            writer.writerow({'locus_tag':iseq.locus_tag, 'srand':iseq.srand,'start':iseq.start, 'end':iseq.end, 'len':iseq.len,
                        'start_codon':iseq.start_codon, 'stop_codon':iseq.stop_codon, 'fpkm':iseq.fpkm, 'tpm':iseq.tpm, 
                        'cov':iseq.cov, 'geneid':iseq.geneid, 'protein_product':iseq.protein_product, 'seq':iseq.seq})
    np.save(para.GENE_DIR+para.bid+id+'csv_gene_dict.npy', gene_dict)
    return gene_dict


def anderson_process():
    anderson_list = []
    itr = 'CTGAAGCTGTCACCGGATGT'
    m35forwad = 'CCAATTATTGAAGGCCTCCCTAACGGGGGGCCTTTTTTTGTTTCTGGTCTCCCgcttGATAAGTCCCTAACT'.upper()[-65:]
    up = m35forwad[-21:-1]
    itr_score= utils.seq2onehot(itr)
    up_score = utils.seq2onehot(up)
    pmt = [] # all kinds training data, including n promoter seq * 4
    pmt_motif = []  # all kinds trainning data, including n promoter seq and motif * 4
    y = []
    f= open(para.TRAIN_DIR+'anderson.txt','r')
    for item in f:
        pitem = para.PMTITEM()
        item = item.split()
        pmt_main = item[1].upper()
        y.append(float(item[2]))
        pitem.promoter = m35forwad+pmt_main+itr
        pitem.itrseq = itr
        pitem.discseq = pmt_main[-6:]
        pitem.m10seq = pmt_main[-12:-6]
        pitem.mextseq = pmt_main[-17:-13]
        pitem.spacerseq = pmt_main[-12-17:-12]
        pitem.m35seq = m35forwad[-1] + pmt_main[:6]
        pitem.mupseq = up
        pitem.spacer = 18

        cur_pmt = utils.seq2onehot(m35forwad+pmt_main+itr) 
        
        pmt.append(cur_pmt)
        disc = utils.seq2onehot(pmt_main[-6:])
        m10 = utils.seq2onehot(pmt_main[-12:-6])
        ext = utils.seq2onehot(pmt_main[-17:-13])
        spacer = utils.seq2onehot(pmt_main[-12-17:-12])
        m35 = utils.seq2onehot(m35forwad[-1] + pmt_main[:6])
        pmt_motif.append(cur_pmt+para.onehot0
                        +itr_score+para.onehot0
                        +disc+para.onehot0
                        +m10+para.onehot0
                        +ext+para.onehot0*2
                        +spacer+para.onehot0
                        +m35+para.onehot0+up_score)
        anderson_list.append(pitem)
    f.close()
    np.save(para.TRAIN_DIR+'anderson_pmt.npy',pmt)
    np.save(para.TRAIN_DIR+'anderson_pmt_motif.npy',pmt_motif)
    np.save(para.TRAIN_DIR+'anderson_y.npy',y)
    np.save(para.TRAIN_DIR+'anderson_list.npy',anderson_list)
    return pmt, pmt_motif, y