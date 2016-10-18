#ipython casoffinder_gen_inputfile.py $INPUT $OUTPUT
import pandas as pd
import sys

#sys.argv[1] : input file
print sys.argv[1]
#sys.argv[2] : outputfile (casoffinder input)
print sys.argv[2]

df_default= pd.read_csv("./ontargets.default", sep="\t", header=None, names=['sgRNAseq'])
df_user= pd.read_csv(sys.argv[1], sep="\t", header=None, names=['sgRNAseq'])

df=pd.concat([df_default, df_user])
df['input_seq']=df.sgRNAseq.apply(lambda x : x+'NNN')
df['mismatch']=6
df[['input_seq','mismatch']].to_csv(sys.argv[2], sep=' ', header=False, index=False)
print "seqeunce writed"

#sys.argv[3] : reference_genome path
with file(sys.argv[2], 'r') as original: data = original.read()
with file(sys.argv[2], 'w') as modified: modified.write(sys.argv[3]+"\nNNNNNNNNNNNNNNNNNNNNNRG 2 1\n" + data)
 
print "header writed"
