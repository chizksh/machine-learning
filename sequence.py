import numpy as np
import pandas as pd
import Levenshtein as lvs

from string import maketrans
def FlipAsciiSeq(ascii_seq):
   trans = maketrans('CGAT', 'GCTA')
   return ascii_seq.translate(trans)[::-1]
def FlipUnicodeSeq(unicode_seq):
   return unicode_seq.translate(dict(zip(u'CGAT', u'GCTA')))[::-1]
def sgRNACAS9distance(x, y):
   arr_x = np.fromstring(x, dtype=np.int8)
   arr_y = np.fromstring(y, dtype=np.int8)
   return np.sum(arr_x!=arr_y)
def CAS9distance(x, y):
   arr_x = np.fromstring(x[20:], dtype=np.int8)
   arr_y = np.fromstring(y[20:], dtype=np.int8)
   return np.sum(arr_x!=arr_y)
def sgRNAdistance(x, y):
   arr_x = np.fromstring(x[:20], dtype=np.int8)
   arr_y = np.fromstring(y[:20], dtype=np.int8)
   return np.sum(arr_x!=arr_y)


def sgRNAlvsdistance(x, y):
    return lvs.distance(x[:20], y[:20])




def CountCleavageSites(cleavage_sites):
    cleavage_count=cleavage_sites.groupby(['ontarget_sgRNAseq'])['mismatch'].value_counts()
    cleavage_count.index.names=['ontarget_sgRNAseq','mismatch']
    cleavage_count=cleavage_count.unstack('mismatch')
    cleavage_count.fillna(0, inplace=True)
    return cleavage_count

def StringEntropy1(seq):
    a=pd.DataFrame(data=list(seq))
    a.columns=['bp']
    ps = a.bp.value_counts() * 0.05
    return ps.map(lambda x : -x * np.log(x)).sum()

def StringEntropy2(seq):
    two_grams = [''.join(x) for x in zip(seq[1:], seq[:-1])]
    b=pd.DataFrame(two_grams)
    b.columns=['two_grams']
    ps = b.two_grams.value_counts() * (1.0/19)#0.0625
    return ps.map(lambda x : -x * np.log(x)).sum()

def GCcontents(seq):
    return (seq.count('C')+seq.count('G'))*1.0/len(seq)
def Tcontents(seq):
    return (seq.count('T'))*1.0/len(seq)

def FillOntargetSeqFeatures(df):
    df['entropy1']=df.index.map(StringEntropy1)
    df['entropy2']=df.index.map(StringEntropy2)
    df['GCcontents']=df.index.map(GCcontents)
    cols_low_mismatch=df.columns[df.columns<4]
    cols_mismatch    =df.columns[df.columns<7]
    df['lowMismatchCleavageRatio']=df[cols_low_mismatch].sum(axis=1)/df[cols_mismatch].sum(axis=1)