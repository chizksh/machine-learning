
# coding: utf-8

# In[ ]:

#!/usr/bin/python


# In[ ]:

import os
import sys
import subprocess
from sys import argv
import pysam
import time
import argparse


# In[ ]:

#default parameter
args = None
chrom_list=[]
path = './'


# In[ ]:

#fuc_steps    
steps = []
fuc_result = ''
def checkTime(func):
    def wrapped(*args, **kwargs):
        try:
            sys.stdout.write('Running: %s...\n' % func.__name__)
            start_time = time.time()
            func(*args, **kwargs)
            sys.stdout.write('Finished.\n')
            end_time = time.time()
            sys.stdout.write("Elapsed Time: %0.2f\n\n"%(end_time - start_time))
        except KeyboardInterrupt:
            sys.stdout.write("KeyboardInterrupt\n")
            raise SystemExit()
        except IOError:
            sys.stdout.write("IOError\n")
            raise SystemExit()
        except Exception as e:
            print str(e)
            sys.stdout.write("%s\n"%(fuc_result))
            
            raise SystemExit()
    steps.append(wrapped)
    return wrapped            


# In[ ]:

#step1
#input: ontarget sequences from user
#output: cas-offinder input file 
@checkTime
def gen_casoffinder_input():
    cmd = 'python {0}casoffinder_gen_inputfile.py {1} {2} {3}'
    fuc_result = subprocess.check_output(cmd.format(path, args.ontarget_seq, '1.'+args.ontarget_seq+'.input', args.ref_genome_path),
                                         stderr=subprocess.STDOUT,shell=True)
    print fuc_result


# In[ ]:

#step2
#input: cas-offinder inputfile
#output: cas-offinder outputfile (bulge,nobulge) ??
@checkTime
def run_casoffinder():
    cmd = 'python cas-offinder-bulge {0}{1} C {0}{2}'
    fuc_result = subprocess.check_output(cmd.format(path, '1.'+args.ontarget_seq+'.input', '2.'+args.ontarget_seq+'.bulge'),
                                         stderr=subprocess.STDOUT,shell=True)
    print fuc_result
    print 'run_casoffinder(assumption)'


# In[ ]:

#step3
#input: cas-offinder outputfile (bulge,nobulge)
#output: cas-offinder outputfile (no bulge only)
@checkTime
def casoffinder_nobulge():
    cmd1 = 'head -n 1 {0}2.{1}.bulge > {0}3.{1}.nobulge'
    cmd2 = 'grep -i "^X" {0}2.{1}.bulge >> {0}3.{1}.nobulge'
    fuc_result1 = subprocess.check_output(cmd1.format(path, args.ontarget_seq),
                                         stderr=subprocess.STDOUT,shell=True)
    print "header writed"
    fuc_result2 = subprocess.check_output(cmd2.format(path, args.ontarget_seq),
                                         stderr=subprocess.STDOUT,shell=True)
    
    print fuc_result1
    print fuc_result2


# In[1]:

#step4
#input: cas-offinder outputfile(bulge)
#output: splitted cas-offinder outputfile (each 4000k line)
@checkTime
def split_casoffinder_bulge():
    input_name = "{0}2.{1}.bulge".format(path, args.ontarget_seq)
    #output_name : input_name.00~99
    #expected exception : file name larger than 99
    split_count = 4000000
    f_in = open(input_name,'r')
    print "load total lines",
    lines = f_in.readlines()
    print "load complete"
    f_in.close()
    header = lines[0]; detail = lines[1:]
    print "total:",
    print len(detail),
    n_chunks = len(detail)/split_count
    print n_chunks
    for i in range(0, n_chunks+1):
        print i,
        f_out = open(input_name+".%02d"%i,'w')
        f_out.write(header)
        for each_line in detail[i*split_count:min(len(detail),(i+1)*split_count)]:
            f_out.write(each_line)
        f_out.close()


# In[ ]:

#step5
#input: cas-offinder nobulge, splited bulge files
#output: bulge,nobulge .p6.v3.csv
@checkTime
def deduplicate_casoffinder():
    cmd = 'python {0}data_preparing_test_new.py {1}'
    fuc_result = subprocess.check_output(cmd.format(path, args.ontarget_seq),
                                         stderr=subprocess.STDOUT,shell=True)
    print fuc_result


# In[ ]:

#step6
#input: ??
#output: mt11.test set, mt11.training set & cleavage_ratio
@checkTime
def labeling_score_file():
    cmd = 'python {0}offtarget_labeling_test.py {1}'
    fuc_result = subprocess.check_output(cmd.format(path, args.ontarget_seq),
                                         stderr=subprocess.STDOUT,shell=True)
    print fuc_result


# In[2]:

#step7
#input: 101rgen.NRG.bulge.mt11.test.csv 
#output: 101rgen.NRG.bulge.mt11.test.msparse_rep
@checkTime
def generate_feature():
    cmd1 = 'python {0}translate_to_features_new.py {0}{1} {0}{2}'
    cmd2 = 'python {0}translate_to_features_new.py {0}{1} {0}{2}'
    
    #cmd1 = 'cat data/101rgen.NRG.bulge.mt11.test.csv | python translate_to_features.py --with-label > data/101rgen.NRG.bulge.mt11.test.msparse_rep'
    #cmd2 = 'cat data/101rgen.NRG.bulge.mt11.train.csv | python translate_to_features.py --with-label > data/101rgen.NRG.bulge.mt11.train.msparse_rep'
    
    fuc_result1 = subprocess.check_output(cmd1.format(path, "6."+args.ontarget_seq+".test.csv", "8."+args.ontarget_seq+".test.msparse_rep"),
                                         stderr=subprocess.STDOUT,shell=True)
    fuc_result2 = subprocess.check_output(cmd2.format(path, "7."+args.ontarget_seq+".train.csv", "9."+args.ontarget_seq+".train.msparse_rep"),
                                         stderr=subprocess.STDOUT,shell=True)
    print fuc_result1
    print fuc_result2


# In[ ]:

#step8
#input: 101RGEN.offtarget.training.ipynb
#output:
@checkTime
def prediction():
    cmd = 'python {0}offtarget.training.py {0}{1} {0}{2} {0}{3} {0}{4} {0}{5}'
    fuc_result = subprocess.check_output(cmd.format(path, "9."+args.ontarget_seq+".train.msparse_rep", "8."+args.ontarget_seq+".test.msparse_rep",                                                   "7."+args.ontarget_seq+".train.csv", "6."+args.ontarget_seq+".test.csv", args.ontarget_seq),
                                         stderr=subprocess.STDOUT,shell=True)
    print fuc_result


# In[ ]:

def machine_learning_run():
    start_time = time.time()
    
    print "start:", args.step_start, "end:", args.step_end
    print range(args.step_start-1, args.step_end)
    
    for st in range(args.step_start-1, args.step_end):
        steps[st]()
    end_time = time.time()

    print('Total Elasped time:%0.2f sec.'%(end_time - start_time))


# In[ ]:

#main
if __name__ == '__main__':
    print 'machine-learning'
    print 'Usage: machine-learning -h or --help'

    #parser
    parser = argparse.ArgumentParser(prog="machine-learning")
    
    #required parameter
    parser.add_argument("ontarget_seq", type=str,
                        help="Specify ontarget_seq file")
    parser.add_argument("ref_genome_path", type=str,
                        help="Specify reference genome(e.g. hg19) path")
    #optional parameter
    parser.add_argument("-p", "--prefix", type=str,
                        help="prefix")
    parser.add_argument("-s", "--step_start", type=int, default=1,
                        help="step start (default: 1)")
    parser.add_argument("-e", "--step_end", type=int, default=8,
                        help="step end (default: 8)")
    #free parameter? 
        
    args = parser.parse_args()
    
    #prefix default
    if args.prefix != None:
        args.prefix = "./"+ args.prefix + "_"
    else:
        args.prefix = "./"        
    
    #machine-learning
    try:
        machine_learning_run()
    except KeyboardInterrupt as detail:
        print "Quitting machine_learning_run.", detail
        raise SystemExit()

