import numpy
import pandas
import sys
import re
import itertools
import thread
import pp
ncpus = 20
ppservers = ()
#job_server = pp.Server(ppservers=ppservers)
job_server = pp.Server()

def SeqToSparseVector(seq):
    arr=numpy.array(list(seq))
    trans = {'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],'T':[0,0,0,1],'N':[1,1,1,1],'-':[0,0,0,0]}
    
    return numpy.array(list(itertools.chain.from_iterable([trans[x] for x in seq])))


def SeqToMSparseVector(ontarget_sgRNA, offsite_sgRNAcas9):
    seq1=ontarget_sgRNA
    seq2=offsite_sgRNAcas9[:20]
    arr1=numpy.array(list(seq1))
    arr2=numpy.array(list(seq2))
    trans = {'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],'T':[0,0,0,1],'N':[1,1,1,1],'-':[0,0,0,0]}
    rep_sgRNA=numpy.concatenate([numpy.outer(numpy.array(trans[x]),numpy.array(trans[y])).flatten() for x,y in zip(arr1,arr2)])
    rep_cas9 = SeqToSparseVector(offsite_sgRNAcas9[-3:])
    return numpy.append(rep_sgRNA,rep_cas9)

    
def line_parallel(line_chunk, col_ontarget_sgRNAseq, col_sgRNAcas9seq, col_bulge, col_score):
    output = ''
    for each_line in line_chunk:
        words=each_line.strip().split(',')
        idx=words[0]
        ontarget_sgRNAseq = words[col_ontarget_sgRNAseq]
        sgRNAcas9seq = words[col_sgRNAcas9seq]
        bulge = words[col_bulge]
        score = words[col_score]
        features = SeqToMSparseVector(ontarget_sgRNAseq, sgRNAcas9seq)
        x_arrstr = numpy.char.mod('%d', features)
        output += ','.join(x_arrstr) + ',%s,%s,%s\n'%(score,bulge,idx)
    #print output
    return output
    
#.csv file
def TransformUnlabeledDataset(f_in_name,f_out_name):
    print f_in_name, f_out_name
    f_in = open(f_in_name, 'r')
    f_out = open(f_out_name, 'w')
    #loading total lines
    print "loading total lines"
    lines = f_in.readlines()
    
    #header writing
    header = ','.join([str(x) for x in range(332)]) +',y,bulge,idx\n'
    f_out.write(header)
    print "header writed"
    
    #column (line1)
    csv_columns=lines[0].strip().split(',')
    col_ontarget_sgRNAseq = csv_columns.index('ontarget_sgRNAseq')
    col_sgRNAcas9seq = csv_columns.index('sgRNAcas9seq')
    col_bulge = csv_columns.index('bulge')
    col_score = csv_columns.index('score')
    #init
    start = 0
    end = len(lines)-1
    print start, end
    job_server.set_ncpus(min(ncpus,(end-start+1)))
    parts = job_server.get_ncpus()
    line_split = numpy.array_split(lines[1:], parts) #line2~
    
    print "Starting with", job_server.get_ncpus(), "workers"
    
    #scatter by file length
    translate_jobs = []
    
    for line_chunk in line_split:
        print len(line_chunk),
        translate_jobs.append(job_server.submit(line_parallel, (line_chunk, col_ontarget_sgRNAseq, col_sgRNAcas9seq, col_bulge, col_score,), (SeqToSparseVector, SeqToMSparseVector,), ("pandas","numpy","itertools",)))
    
    print "concat"
    outputs = '\n'.join([test_job() for test_job in translate_jobs])
    job_server.print_stats()
    
    f_out.write(outputs)
    print "detail writed"
    
    f_in.close()
    f_out.close()

TransformUnlabeledDataset(sys.argv[1],sys.argv[2])