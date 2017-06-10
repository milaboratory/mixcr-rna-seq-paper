import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tre
import os
import re
import matplotlib.gridspec as gridspec

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif']=["Arial"] 


# Global options
ROOT_DIRECTORY="./"
STAR_PATH="%s/star"%ROOT_DIRECTORY
MIXCR_PATH="%s/mixcr"%ROOT_DIRECTORY
TRUST_PATH="%s/trust"%ROOT_DIRECTORY
FIGURES="%sfigures"%ROOT_DIRECTORY


def parse_true_reads(chain):
    """
    Parse initial in-silico generated V-D-J reads into DataFrame
    
    Parameters:
        chain -- the immunological chain (TRA/IGH/...)
    """
    
    file_name = '%s/in_silico_%s.fasta'%(ROOT_DIRECTORY,chain)
    
    data = {
        "nSeqCDR3": [], # nucleotide CDR3 sequence
        "aaSeqCDR3": [] # amino acid CDR3 sequence
    }
        
    with open(file_name, 'r', 4096) as f:
        content = f.readlines()
        
        for line in content:
            if not line.startswith('>'):
                continue
            split = line.split("|")
            data['nSeqCDR3'].append(split[1])
            data['aaSeqCDR3'].append(split[2])
    data = pd.DataFrame.from_dict(data)
    data = data.groupby("nSeqCDR3").first()
    data = data.reset_index()

    return data

def parse_bam_file_name(fname):
    """
    Auxiliary method for parsing a file name from star/* directory

    Arguments:
        fname -- the BAM file name
        
    Returns:
        {
            vdj    : true/false
            len    : read length (bp)
            ref    : hg38/hg37
            paired : paired/single-end input
        }
    """
        
    if not fname.endswith('.bam'):
        raise Exception()
    
    fname = fname.replace('.bam', '')
    fname = fname.replace('.sorted', '')
    
    sample = {'sample_name' : fname.replace('in_silico_RNA_Seq_','')}
    patt = re.search('in_silico_RNA_Seq_*([a-zA-Z_]*)_([0-9]+)bp\.(hg3[78])\.([a-z]+)', fname)
    sample['vdj'] = patt.group(1) == ''
    sample['len'] = int(patt.group(2))
    sample['ref'] = patt.group(3)
    sample['paired'] = patt.group(4) == 'paired'
    return sample

def parse_mixcr(sample):
    """
    Parse MiXCR results into DataFrame for specified sample 

    Arguments:
        sample -- the sample info of the form:
                    {
                        vdj    : true/false
                        len    : read length (bp)
                        ref    : hg38/hg37
                        paired : paired/single-end input
                    }
        
    Returns:
        DataFrame with MiXCR results obtained for the sample
    """
    
    return pd.read_table('%s/in_silico_RNA_Seq_%s%sbp.%s.txt'%(
        MIXCR_PATH, '' if sample['vdj'] else 'no_VDJ_', 
        sample['len'], 'paired' if sample['paired'] else 'single'))

def safe_float(s):
    try:
        return float(s)
    except ValueError:
        return float('nan')
    
def parse_trust(sample):
    """
    Parse TRUST results into DataFrame for specified sample 

    Arguments:
        sample -- the sample info of the form:
                    {
                        vdj    : true/false
                        len    : read length (bp)
                        ref    : hg38/hg37
                        paired : paired/single-end input
                    }
        
    Returns:
        DataFrame with TRUST results obtained for the sample
    """
    
    fileName = '%s/in_silico_RNA_Seq_%s%sbp.%s.%s.sorted.bam.fa'%(
        TRUST_PATH, '' if sample['vdj'] else 'no_VDJ_', 
        sample['len'], sample['ref'], 'paired' if sample['paired'] else 'single')
    data = {
        "contig": [],
        "V": [],
        "J": [],
        "nSeqCDR3": [],
        "aaSeqCDR3": [],
        "freq": []
    }
    with open(fileName, 'r', 4096) as f:
        content = f.readlines()
        
        prev = None
        for line in content:

            if line.startswith("#"):
                continue

            if prev and prev.startswith(">"):
                split = prev[1:-1].split("+")
                data["contig"].append(line[:-1])
                data["V"].append(split[4])
                data["J"].append(split[5])
                data["nSeqCDR3"].append(split[9])
                data["aaSeqCDR3"].append(split[7])
                data["freq"].append(safe_float(split[1][16:]))

            prev = line

    trust = pd.DataFrame.from_dict(data)
    del trust["freq"]
    return trust


# will use only TRB clones from MiXCR
def parse_mixcr_tcr(sample):
    mixcr = parse_mixcr(sample)
    return mixcr[mixcr['allVHitsWithScore'].str.contains('^TR[ABDG]V')]

def parse_mixcr_ig(sample):
    mixcr = parse_mixcr(sample)
    return mixcr[mixcr['allVHitsWithScore'].str.contains('^IG[HKL]V')]


class TrueClonesDb:
    """
    A database of the initial in-silico generated clones which 
    allows fast search 
    
    """
    
    def __init__(self, true_clones_df, nt=True):
        self.true_clones_df = true_clones_df
        
        column = 'nSeqCDR3' if nt else 'aaSeqCDR3'
        patterns_end_to_end = []
        patterns_any = []
        for row in true_clones.iterrows():
            patterns_end_to_end += [(row[1][column], tre.compile("^" + row[1][column] + "$"))]
            #patternsAny += [tre.compile(row[1].cdr3)]
            patterns_any += [(row[1][column], tre.compile(row[1][column][3:-3]))]

        self.patterns_end_to_end = patterns_end_to_end
        self.patterns_any = patterns_any
        
    def search_cdr3_seq(self, seq, mmaxerr, end_to_end=True):
        for maxerr in range(0, 1 + mmaxerr):
            fuzzyness = tre.Fuzzyness(maxerr = maxerr)
            for p in (self.patterns_end_to_end if end_to_end else self.patterns_any):
                if p[1].search(seq, fuzzyness):
                    return p[0]    
        return np.nan


"""
Parse all available BAM files from star/ directory
"""
AllSamples = [parse_bam_file_name(f) for f in os.listdir(STAR_PATH) if f.endswith('bam')]
# todo: inferr chains automatically
# true_clones = parse_true_reads('TRB').append(parse_true_reads('IGH'), verify_integrity=True, ignore_index=True)
true_clones = parse_true_reads('TRB')
true_clones_nt_db = TrueClonesDb(true_clones, nt=True)
true_clones_aa_db = TrueClonesDb(true_clones, nt=False)



MAX_MAXERR=5

def get_stats(df_table, nt=True, prefix=''):
    result = { '%stotal'%prefix : len(df_table) }
    column = 'nSeqCDR3' if nt else 'aaSeqCDR3'
    db = true_clones_nt_db if nt else true_clones_aa_db
    for maxerr in range(0, MAX_MAXERR):
        tmp = df_table[column].apply(lambda seq: db.search_cdr3_seq(seq, maxerr))        
        result['%smatched_clones_%s'%(prefix,maxerr)] = pd.notnull(tmp.unique()).sum()
        result['%smatched_records_%s'%(prefix,maxerr)] = pd.notnull(tmp).sum()
        result['%sunmatched_clones_%s'%(prefix,maxerr)] = len(df_table.loc[tmp.isnull(), column].unique())
        result['%sunmatched_records_%s'%(prefix,maxerr)] = tmp.isnull().sum()
    return result


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

AllStats = []
for sample in AllSamples:
    mixcr_stats = get_stats(parse_mixcr_tcr(sample), prefix='mixcr_')
    trust_stats = get_stats(parse_trust(sample), prefix='trust_')
    AllStats.append(merge_dicts(sample, mixcr_stats, trust_stats))
    
AllStats = pd.DataFrame.from_records(AllStats)




matplotlib.style.use('default')
def plot_false_positive_in_silico_left_only(paired=True, star_reference = 'hg37', ax = None):
    if ax is None:
        ax = plt.subplot()

    ax_true_vdj  = ax
    
    lengths   = [50, 75, 100]
    softwares = ['MiXCR', 'TRUST']
    software_offset = 2.5
    width           = 1.5 
    length_offset   = software_offset * len(softwares) + width
    
    all_data = AllStats[(AllStats.paired == paired) & (AllStats.ref == star_reference)]
    re_total="(" + "|".join([f.lower() for f in softwares]) + ")_total"

    max_total_true_vdj = all_data[all_data.vdj == True].filter(regex=re_total).max().max()
    max_total_false_vdj = all_data[all_data.vdj == False].filter(regex=re_total).max().max()
    max_total_all = max(max_total_true_vdj, max_total_false_vdj)
    labels = []
    label_positions = []
    for true_vdj in [True]:
        ax   = ax_true_vdj if true_vdj else ax_false_vdj
        delta = 0 if true_vdj else length_offset * 3

        for i in range(len(lengths)):
            length = lengths[i]
            
            maxy = max_total_true_vdj if true_vdj else max_total_false_vdj
            lendata =  all_data[(all_data['len'] == length) & (all_data['vdj'] == true_vdj)]
            
            gmax = 0
            for j in range(len(softwares)):

                software = softwares[j].lower().replace('\'', '')
                begin = delta + i * length_offset + j * software_offset

                labels.append(softwares[j])
                label_positions.append(begin)

                data = lendata.filter(regex=software)
                matched_0 = data['%s_matched_clones_0'%software].values[0]
                matched_1 = data['%s_matched_clones_1'%software].values[0]
                matched_2 = data['%s_matched_clones_2'%software].values[0]
                matched_3 = data['%s_matched_clones_3'%software].values[0]
                matched_4 = data['%s_matched_clones_4'%software].values[0]
                unmatched = data['%s_unmatched_clones_4'%software].values[0]
                gmax = max(gmax, unmatched + matched_4)

                gr_col = (0.24000000000000016, 0.5819607843137254, 0.2180392156862746) 
                r_col  = (0.8113725490196079, 0.1807843137254901, 0.18705882352941164) 
                ax.bar([begin], [matched_0 -    0     ], bottom = [    0    ], color=gr_col, alpha=1.0, width=width)
                ax.bar([begin], [matched_1 - matched_0], bottom = [matched_0], color=gr_col, alpha=0.8, width=width)
                ax.bar([begin], [matched_2 - matched_1], bottom = [matched_1], color=gr_col, alpha=0.6, width=width)
                ax.bar([begin], [matched_3 - matched_2], bottom = [matched_2], color=gr_col, alpha=0.4, width=width)
                ax.bar([begin], [matched_4 - matched_3], bottom = [matched_3], color=gr_col, alpha=0.2, width=width)
                if unmatched != 0:
                    ax.bar([begin], [       unmatched     ], bottom = [matched_4], color=r_col, alpha=1.0, width=width)
        
            ax.text(delta + i * length_offset + software_offset * len(softwares) / 2 - width / 2, 
                    gmax + maxy / 30.0 if maxy != 0 else max_total_all / 10.0, str(length) + 'bp', 
                    ha='center', va='bottom', fontsize=13, color='k')

    
    ax_true_vdj.set_xticks(label_positions)
    ax_true_vdj.set_xticklabels(labels, rotation=75, ha='center', fontsize=11)
    ax_true_vdj.get_xaxis().set_tick_params(length=0)
    
    ax_true_vdj.set_ylim(0, 1.3 * max_total_true_vdj)
    ax_true_vdj.set_ylabel("Number of CDR3 clones", fontsize=13)


def panel_letter(ax, letter):
        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=19, fontweight='bold', va='top', ha='right')



plt.figure(figsize=(14,5))

left = plt.subplot(121)
plot_false_positive_in_silico_left_only(paired=True, star_reference='hg37',ax=left)
plt.title('Paired-end')
panel_letter(left, 'a')

right=plt.subplot(122)
plot_false_positive_in_silico_left_only(paired=False, star_reference='hg37',ax=right)
plt.title('Single-end')
panel_letter(right, 'b')

plt.subplots_adjust(wspace=0.45)
plt.savefig("SupplementaryFigure_in_silico.pdf")

