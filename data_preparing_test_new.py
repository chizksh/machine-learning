import pandas
import numpy
import thread
import pp
import glob
import sys
ncpus = 20

def FormattingCasoffinderOutput(df):
    df.columns = [u'bulge_type','crRNA','DNA','chroms','Position', u'strand', u'mismatch', u'bulge_len']
    df['bulge']=numpy.where(df['bulge_len']>0,1,0)
    df['ed'] =df.mismatch+  df.bulge
    df['penalty']=df.mismatch+2*df.bulge
    return df

def CleavageOffset(elm):
    seq=elm.crRNA
    tmp=numpy.where(numpy.array(list(seq))=='-',0,1).cumsum()
    clv_offset=numpy.where(tmp==17)[0][0]+1
    dna_offset_factor=elm.DNA[:clv_offset].count('-')
    if(elm.strand=='+'):
        clv_start = elm.bed_start+(clv_offset-dna_offset_factor)
    else:
        clv_start = elm.bed_end-(clv_offset-dna_offset_factor)
    return clv_start

def FillColumns(df):
    df['bed_start']=numpy.where(df.strand=='+', df.Position,\
                             numpy.where(df.bulge_type=='RNA',\
                                      df.Position-4+df.bulge_len,df.Position-2+df.bulge_len))
    df['bed_end']  =numpy.where(df.bulge_type=='DNA',\
                             df.bed_start+23+df.bulge_len,\
                             df.bed_start+23-df.bulge_len)
    df['clv_start']=df.apply(CleavageOffset, axis=1)
    df['clv_end'  ]=df.clv_start+1
    df['ontarget_sgRNAseq']=df.crRNA.str.upper().str.replace('-','').apply(lambda x: x[:-3])  
    
def ProcessBulge(df_dedup):
    df_dedup=df_dedup[df_dedup.DNA.apply(lambda x : len(x)) >22]
    df_dedup['sgRNAcas9seq']=df_dedup.apply(lambda x: numpy.array(list(x.DNA))[numpy.array(list(x.crRNA))!='-'].tostring().upper(), axis=1)
    return df_dedup

    
    
    
ppservers = ()
#job_server = pp.Server(ppservers=ppservers)
job_server = pp.Server()

def offtarget6(file_name, index_lst): ##1~4
    #0.file_name
    infiles = []
    for index in index_lst:
        infiles += [file_name+".%02d"%index]
    #1. read data from file
    df=pandas.concat([pandas.read_csv(infile,sep='\t') for infile in infiles])
    #2. FormattingCasoffinder
    df=FormattingCasoffinderOutput(df)
    #3. copy of penalty < 7 --> no_bulge6
    offtarget6=df[df.penalty<7].copy()
    #4. FillColumns
    FillColumns(offtarget6)
    
    return offtarget6
    
def processbulge(df, start, end): ##6. process_bulge
    df_dedup=ProcessBulge(df[start:end])
    return df_dedup
    

#main
if __name__ == '__main__':
    ontarget_seq = sys.argv[1]
    no_bulge_file = "./3."+ontarget_seq+".nobulge"
    bulge_file = "./2."+ontarget_seq+".bulge.[0-9][0-9]"
    
    ##no_bulge part
    print "no_bulge"
    #1. read data from file
    nobulge= pandas.read_csv(no_bulge_file, sep='\t')
    #2. FormattingCasoffinder
    nobulge=FormattingCasoffinderOutput(nobulge)
    #3. copy of penalty < 7 --> no_bulge6
    nobulge6=nobulge[nobulge.penalty<7].copy()
    #4. FillColumns
    FillColumns(nobulge6)
    #5. sort_by_value, and get top of them --> drop_duplicate
    nobulge6_dedup=nobulge6.sort_values(['penalty','bulge_len']).drop_duplicates(['ontarget_sgRNAseq','chroms','clv_start']).copy()
    #6. Process_bulge
    nobulge6_dedup=ProcessBulge(nobulge6_dedup)
    #7. write data to file
    nobulge6_dedup.to_csv("./5."+ontarget_seq+".nobulge.p6.v3.csv", index=False)
    
    
    ##bulge part
    print "bulge"
    #init
    start = 0
    end = len(glob.glob(bulge_file))-1
    print start, end
    job_server.set_ncpus(min(ncpus,(end-start+1)))
    parts = job_server.get_ncpus()
    print "Starting with", job_server.get_ncpus(), "workers"
    
    #scatter by bulge_files
    
    index_split = numpy.array_split(range(start,end+1),parts)
    file_name = "./2."+ontarget_seq+".bulge"
    
    
    #offtarget6
    print "offtarget6"
    #scatter by filenum
    offtarget6_jobs = []
    for index_lst in index_split:
        print index_lst
        offtarget6_jobs.append(job_server.submit(offtarget6, (file_name, index_lst,), (FormattingCasoffinderOutput,FillColumns,CleavageOffset,), ("pandas","numpy",)))
    
    print "concat_offtargets"
    offtargets = pandas.concat([offtarget_job() for offtarget_job in offtarget6_jobs])
    
    job_server.print_stats()
    
    #drop_duplicate
    print "drop_dup"
    offtargets_dedup = offtargets.sort_values(['penalty','bulge_len']).drop_duplicates(['ontarget_sgRNAseq','chroms','clv_start']).copy()
    
    #processbulge
    print "processbulge"
    processbulge_jobs = []
    
    #scatter by offtarget_dedup
    start = 0
    end = len(offtargets_dedup)-1
    print start, end
    df_split = numpy.array_split(range(start,end+1), parts)
    
    for df_lst in df_split:
        starti = df_lst[0]
        endi= df_lst[-1]+1
        print starti
        processbulge_jobs.append(job_server.submit(processbulge, (offtargets_dedup, starti, endi,),(ProcessBulge,),("pandas","numpy",)))
    
    print "concat_processbulge"
    offtargets_dedup2 = pandas.concat([processbulge_job() for processbulge_job in processbulge_jobs])
    
    job_server.print_stats()
    
    #check
    #print "check"
    #assert(numpy.all(offtargets_dedup2.Position == offtargets_dedup.Position))

    #write
    print "write"
    offtargets_dedup2.to_csv("./4."+ontarget_seq+".bulge.p6.v3.csv", index=False)
    
    
    





