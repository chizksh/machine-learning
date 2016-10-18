import pandas as pd
import numpy as np
import sys
#import FindNearsestOffsites as *


class FindNearsestOffsites:
    def __init__(self, offsites):
        self.offsites=offsites.sort_values(by=['clv_start'])
        self.clv_sites={}
        for key, sites in self.offsites.groupby(['chroms'])['clv_start']:
            self.clv_sites[key]=sites
    def Find(self, chroms, sites):
        clv_sites=self.clv_sites[chroms]
        idxs=clv_sites.searchsorted(sites)
        idxs[idxs==len(clv_sites)]=0
        #mask= (idxs<len(clv_sites))&(idxs>0)
        #idxs=idxs[mask]
        #clv_sites=clv_sites[mask]
        lefts=clv_sites.iloc[idxs-1]
        rights=clv_sites.iloc[idxs]
        idx_nearby=np.where(np.abs(rights.values-sites.values)<np.abs(lefts.values-sites.values), 
                            idxs, idxs-1)
        return clv_sites.iloc[idx_nearby]
    def FindAll(self, cleavege_sites):
        def helper():
            for chroms, df in cleavege_sites.groupby('chroms'):
                nearby=self.Find(chroms, df.clv_start)
                df=df.copy()
                df['idx']=nearby.index
                df['site_nearby']=nearby.values
                yield df
        df=pd.concat(helper(), axis=0)
        df['offset']=np.abs(df.clv_start-df.site_nearby)
        return df
    def __call__(self, cleavege_sites, ordering, offset_cut):
        df=self.FindAll(cleavege_sites)
        
        df=df[df.offset<offset_cut]
        df=df[['score','offset','idx']].merge(self.offsites, left_on='idx', right_index=True)
        print 'before:', df.shape
        df=df.sort_values(ordering).drop_duplicates(['chroms','clv_start'])
        print 'after:', df.shape
        return df

		
#import load_ontargets as *


#ontargets_rgen101 = pd.read_csv('data/ontargets.101rgen', header=None, names=['seq'])

ontarget_seq = sys.argv[1]
ontargets_user = pd.read_csv("./"+ontarget_seq, header=None, names=['seq'])
#assert(ontargets_val1test.seq.isin(ontargets_user.seq).sum()+ontargets_val1train.seq.isin(ontargets_user.seq).sum()==9)

#ibs_seqs=pd.read_csv('data/seq_name.csv')
#seq_name=dict(zip(ibs_seqs.seq, ibs_seqs.name))
#seq_name.setdefault(np.NaN)
#name_seq=dict(zip(ibs_seqs.name, ibs_seqs.seq))

#rank_val = pd.read_csv('data/rank_val.csv', sep='\t', header=None, names=['gene','seq'])
		
		
#part1 - Load

#load raw_data(digenome, guideseq)
score_ibs=pd.read_csv('./101rgen.score')
score_guideseq = pd.read_csv('./guide-seq.score')

#add mt11 score to guide-seq
score_guideseq['score']=score_guideseq.score_GuideSeq
#score_guideseq_mt11 = score_guideseq[score_guideseq.ontarget_sgRNAseq.isin(ontargets_user.seq)]

#load data_preparing output(from previous step)
offsites = pd.read_csv("./4."+ontarget_seq+".bulge.p6.v3.csv")
offsites_nobulge = pd.read_csv("./5."+ontarget_seq+".nobulge.p6.v3.csv")

"""
#load raw_data null exp / digenome / guide-seq
mt0=pd.read_csv('data/mt11rgen.score')
rgen0=pd.read_csv('data/101rgen.score')
null0=pd.read_csv('data/101rgen.nullexp.score')
guideseq=pd.read_csv('data/guide-seq.score')
"""

#part2 - Calculate dNOC

#1. FindNearsestOffsites class initialize
#initialize: load clv_sites by chromosome -> dict
bulgeOffinder = FindNearsestOffsites(offsites)
nobulgeOffinder = FindNearsestOffsites(offsites_nobulge)

#2. Getting d_nearest 
"""
#calculate offset value by offtarget site -> export csv file
#100rgen null
nobulgeOffinder.FindAll(null0).to_csv('data/100rgen.nobulge.score.null.offset.csv', index=False)
bulgeOffinder.FindAll(null0).to_csv('data/100rgen.score.null.offset.csv', index=False)
#100rgen digenome
nobulgeOffinder.FindAll(rgen0).to_csv('data/100rgen.nobulge.score.offset.csv', index=False)
bulgeOffinder.FindAll(rgen0).to_csv('data/100rgen.score.offset.csv', index=False)
#11rgen digenome
nobulgeOffinder.FindAll(mt0).to_csv('data/11rgen.nobulge.score.offset.csv', index=False)
bulgeOffinder.FindAll(mt0).to_csv('data/11rgen.score.offset.csv', index=False)
#guide-seq digenome
nobulgeOffinder.FindAll(guideseq).to_csv('data/guideseq.nobulge.score.offset.csv', index=False)
bulgeOffinder.FindAll(guideseq).to_csv('data/guideseq.score.offset.csv', index=False)
"""

#find nearest offtarget site -> drop_duplicate

"""
offset_cut=2
a =nobulgeOffinder(rgen0, 'mismatch',offset_cut)
a0=nobulgeOffinder(null0, 'mismatch',offset_cut)
b =bulgeOffinder(rgen0, 'penalty',offset_cut)
b0=bulgeOffinder(null0, 'penalty',offset_cut)
a0.to_csv('data/100rgen.nobulge.offset_cut.null.offset.csv', index=False)
b0.to_csv('data/100rgen.offset_cut.null.offset.csv', index=False)
a.to_csv('data/100rgen.nobulge.offset_cut.offset.csv', index=False)
b.to_csv('data/100rgen.offset_cut.offset.csv', index=False)
"""
#101rgen/guide-seq
offtargets_ibs=bulgeOffinder(score_ibs, 'penalty',2)
offtargets_guideseq=bulgeOffinder(score_guideseq, 'penalty',2)
#offtargets_nobulge_ibs=nobulgeOffinder(score_ibs, 'mismatch',2)
#offtargets_nobulge_guideseq=nobulgeOffinder(score_guideseq, 'mismatch',2)
#offtargets_guideseq_mt11=bulgeOffinder(score_guideseq_mt11, 'penalty',2)
#merge(ibs,guideseq)
offtargets_with_score=offsites.merge(offtargets_ibs[['score','offset','idx']].set_index('idx'), left_index=True, right_index=True, how='left')
offtargets_with_scores=offtargets_with_score.merge(offtargets_guideseq[['score','offset','idx']].add_suffix('_guideseq').set_index('idx_guideseq'), left_index=True, right_index=True, how='left')
#export to file
offtargets_with_scores.to_csv('./101rgen.NRG.bulge.p6.csv', index=False)


#part3 - Prepare ontargets and split data into test&train sets.
rawdata_test2=offtargets_with_scores[offtargets_with_scores.ontarget_sgRNAseq.isin(ontargets_user.seq)]
rawdata_train2=offtargets_with_scores[-offtargets_with_scores.ontarget_sgRNAseq.isin(ontargets_user.seq)]
rawdata_test2.to_csv("./6."+ontarget_seq+".test.csv", index=False)
rawdata_train2.to_csv("./7."+ontarget_seq+".train.csv", index=False)