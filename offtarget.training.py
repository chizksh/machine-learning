import numpy as np
import pandas as pd
import sys

from rgen_offtarget_learning import *

trainset_name = sys.argv[1] #train.msparse_rep
testset_name = sys.argv[2]  #test.msparse_rep
train_info_name = sys.argv[3] #.train.csv
test_info_name = sys.argv[4] #.test.csv

id, score_cut, n_epoch = 0, 3000,60

offset_cut=2
dataset_name = sys.argv[5]

trainset= pd.read_csv(trainset_name).drop('idx', axis=1)
testset = pd.read_csv(testset_name).drop('idx', axis=1)
train_info= pd.read_csv(train_info_name)
test_info = pd.read_csv(test_info_name)
assert(np.all(trainset[trainset.y.notnull()].y == train_info.score[train_info.score.notnull()]))

score_cutoff=lambda x: (x.score>score_cut)&(x.offset<offset_cut)
trainset.y=np.where(score_cutoff(train_info), 1, 0)
testset.y=np.where(score_cutoff(test_info), 1, 0)
    
n_neurons0=2048
p_dropout0=0.5
learning_rate=0.002

learner=RGENPrediction(LasagneNN1(2, trainset.shape[1]-1, n_neurons0,p_dropout0,learning_rate,n_epoch))
#Note that name of a label column of dataset must be 'y'.
learner.Train(trainset)

ps_test = learner.model.predict_proba(testset.drop('y',axis=1)).T[1]
test_info['ibs'] = np.log(ps_test)
test_info['PAM']=test_info.sgRNAcas9seq.apply(lambda x : x[-3:])
col_out = [u'chroms', u'strand', u'mismatch', 'PAM',u'bulge', u'ed', u'penalty', u'clv_start',\
           u'ontarget_sgRNAseq', u'score', u'offset', u'score_guideseq', u'offset_guideseq', u'ibs']
test_info[col_out].to_csv(dataset_name+'.test.pred.%d'%id, index=False)