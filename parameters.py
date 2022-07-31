GTF_DIR = './00_genome/'
DEPTH_DIR = './01_depth/'
TRAIN_DIR = './02_train/'
TSS_DIR = './04_tss/'
PROMOTER_DIR = './05_promoterseq/'
FASTA_DIR = './06_fasta/'
RESULT_DIR = './07_result/'
GENE_DIR = './08_gene_data/'
MOTIF_DIR = './09_motif/'
FIG_DIR = './10_fig/'
LOGO_DIR = './11_logo/'
STIE_DIR = './13_stringtie/'
bid = 'MG1655_9132'
sid = 'SRR1173983'


bacteria = {'MG1655_9132':['Ecoli_K12_MG1655_NC_000913.2.fasta','Ecoli_K12_MG1655_NC_000913.2.gff3','proteins_167_161521.csv'],
            'MG1655_9133':['','','proteins_167_161521.csv'],
            'Bacillus168':['Bacillus_subtilis_str_168.fna','Bacillus_subtilis_str_168.gff','Bacillus_subtilis_str_168.csv'],
            'Cor_glutamicum':['Corynebacterium_glutamicum.fna','','Corynebacterium_glutamicum.csv']} 

genome_file = bacteria[bid][0]
gff_file = bacteria[bid][1]
protein_csv_file = bacteria[bid][2]




shortlist = ['GSE53767', 'wig']
S986list = ['SRR1173970', 'depth.txt']
PRJNA827912list = ['SRR1232437', 'depth.txt'] # bacillus
PRJNA798895list = ['SRR19091943', 'depth.txt'] # glu
BS168list = ['SRR6435132',
'SRR6435133',
'SRR6435134',
'SRR6435135',
'SRR6435136',
'SRR6435137',
'SRR6435138',
'SRR6435139',
'SRR6435140',
'SRR6435141',
'SRR6435142',
'SRR6435143',
'SRR6435144',
'SRR6435145']

dp_dict = {'wig':1, 'depth.txt':0}

input_list = PRJNA798895list
# input_list[0]
suffix = input_list[1]
dp_para = dp_dict[suffix]

short_id = shortlist[0]
short_suffix = shortlist[1]
short_dp = dp_dict[short_suffix]

BS168_LIST = []
CGLU_LIST = []

CORR_LIST = ['SRR10907652', 'SRR10907653', 'SRR10907655', 'SRR10907660', 'SRR10907661', 'SRR10907665', 'SRR10907668',
                'SRR10907669', 'SRR10920132']
EXP_LIST = ['SRR1173965', 'SRR1173966', 'SRR1173967', 'SRR1173968', 'SRR1173969', 'SRR1173970', 'SRR1173981', 'SRR1173982',
            'SRR1173983', 'SRR1173984', 'SRR1173985', 'SRR1173986']
STA_LIST = ['SRR1173971', 'SRR1173972', 'SRR1173973', 'SRR1173974', 'SRR1173975', 'SRR1173976', 'SRR1173977', 'SRR1173978',
            'SRR1173979', 'SRR1173980']
EXP_LIST_NEG = ['SRR1173965', 'SRR1173967', 'SRR1173968', 'SRR1173981', 'SRR1173983', 'SRR1173984']
STA_LIST_NEG = ['SRR1173971', 'SRR1173972', 'SRR1173975', 'SRR1173976', 'SRR1173977']
EXP_LIST_POS = ['SRR1173966', 'SRR1173969', 'SRR1173970', 'SRR1173982', 'SRR1173985', 'SRR1173986']
STA_LIST_POS = ['SRR1173973', 'SRR1173974', 'SRR1173978', 'SRR1173979', 'SRR1173980']

stage_dict = {'exp':{'neg':EXP_LIST_NEG, 'pos':EXP_LIST_POS}, 'sta':{'neg':STA_LIST_NEG,'pos':STA_LIST_POS}}

# sid = 'SRR18803111'

MG1655_protein_csv = './00_genome/proteins_167_161521.csv'
paper_tss_csv = './00_genome/MG1655_tss_all_paper.csv'
MG1655_9132_gtf_file = './00_genome/Ecoli_K12_MG1655_NC_000913.2.gff3'
protein_start = ['ATG','GTG','TTG']
protein_end = ['TAA','TAG','TGA']

as_dict = {'+':'-', '-':'+'}
klist = ['s', 'i', 'a']
plot_klist = ['p', 's', 'i']
kalllist = ['p', 's', 'i', 'a']

terminus = 1600000
origin = 3900000
search_area = 300
promoter_len = 100
nmotifs = 1
str_search_len = 200

INF_MIN = float('-inf')




s10 = 75
e10 = 99
onehot0 = [0,0,0,0]
data_kind = {0:'p',1:'p+s',2:'p+s+i'}
class GENEITEM():
    def __init__(self):
        #from gtf
        self.start = -1
        self.end = -1
        self.len = 0
        self.srand = '+'
        self.geneid = ''
        self.locus = ''
        self.locus_tag = ''
        self.protein_product = ''
        self.start_codon = ''
        self.stop_codon = ''
        self.fpkm = 0
        self.tpm = 0
        self.cov = 0
        #genome
        self.seq = []
        #compute
        self.strength = 0
        #algorithm
        self.tss = -1
        self.ptss = -1
        self.itss = []
        self.stss = []
        self.astss = []

    def get_start(self, start, end, g_len):
        if self.srand == '+':
            self.start = start - 1
        else:
            self.start = g_len - end
        return self.start

    def get_end(self, start, end, g_len):
        if self.srand == '+':
            self.end = end - 1
        else:
            self.end = g_len - start
        return self.end

    def set_len(self, start, end):
        self.len += abs(end - start) + 1

    def set_srand(self, srand):
        self.srand = srand

    def set_seq(self,gdict,start, end):
        self.seq = gdict[self.srand][start:end+1]

class PMTALL():
    def __init__(self):
        self.allp = {'p':[], 's':[], 'i':[], 'a':[]}

class PMTITEM():
    def __init__(self,iseq = GENEITEM(),pmt = '',tss = 0,kind = ''):
        self.srand = iseq.srand
        self.kind = kind
        self.tss = tss
        self.ctss = 0
        self.lstrength = 0
        self.sstrength = 0
        self.m10p = 0
        self.m35p = 0
        self.spacer = 0
        #p-value
        self.m10pvalue = 0.0
        self.m35pvalue = 0.0
        self.mextpvalue = 0.0
        #score
        self.m10score = 0.0
        self.m35score = 0.0
        self.mextscore = 0.0
        self.spacerscore = 0.0
        self.upscore = 0.0
        self.itrscore = 0.0
        self.discscore = 0.0
        self.mscore = 0
        #motif seq
        self.discseq = ''
        self.m10seq = ''
        self.m35seq = ''
        self.mextseq = ''
        self.mupseq = ''
        self.spacerseq = ''
        self.itrseq = ''
        self.promoter = pmt
        self.trainseq = ''