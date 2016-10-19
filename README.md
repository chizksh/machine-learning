Our program should be run under OpenCL enabled Linux system.
You might need "sudo" right to install neccessary libraries and run our program.

1. Dependency

- python 2.x version (https://www.python.org/downloads/) 
- python libraries: You can install using "pip" command

  1. numpy (http://www.scipy.org/scipylib/download.html)
  2. f2py (https://sysbio.ioc.ee/projects/f2py2e/#download)
  3. cython (http://cython.org/#download)
  4. scipy (https://www.scipy.org/install.html)
  5. pandas (https://pypi.python.org/pypi/pandas/0.19.0/#downloads)
  5. lasagne (https://github.com/Lasagne/Lasagne)
  6. sklearn (http://scikit-learn.org/stable/install.html)
  7. nolearn (https://pypi.python.org/pypi/nolearn)
  8. python-Levenshtein (https://pypi.python.org/pypi/python-Levenshtein/0.12.0)
  9. pp (http://www.parallelpython.com/content/view/18/32/)
  10. pysam (https://pypi.python.org/pypi/pysam)

  
- Cas-offinder binary (https://sourceforge.net/projects/cas-offinder/files/Binaries/)
  * Cas-offiner requirement : download and install a proper OpenCL SDK to install runtime API
   (AMD: http://developer.amd.com/tools-and-sdks/heterogeneous-computing/amd-accelerated-parallel-processing-app-sdk/downloads/)
   (Intel: http://software.intel.com/en-us/vcsource/tools/opencl-sdk)
   (NVidia: https://developer.nvidia.com/cuda-downloads)
  ** For successful installation of OpenCL SDK, recent version of lsb-core package of specific linux is required.

- fasta(*.fa) file of each chromosomes of hg19 (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)
  

2. Installation & Instruction

1. install python 2.x
2. install python libararies
3. unzip master.zip file into any directory
  -master.zip should include these files
    (main source codes)
    machine-learning-run.py
    casoffinder_gen_inputfile.py
    casoffinder_gen_inputfile.sh
    data_preparing_test_new.py
    machine-learning-run.py
    offtarget.training.py
    offtarget_labeling_test.py
    rgen_offtarget_learning.py
    sequence.py
    translate_to_features_new.py
    cas-offinder bulge

    (data files)
    101rgen.score
    guide-seq.score
    ontargets.default

    (input file example)
    test.seq
4. download cas-offinder binary and move into the same directory
5. download and unzip hg19 chromosome files in the same directory


3. Usage
Usage: python machine-learning-run.py -h or --help
usage: python machine-learning-run.py [-h] [-p PREFIX] [-s STEP_START] [-e STEP_END]
                        ontarget_seq ref_genome_path
(e.g. python machine-learning-run.py ./test.seq /data/hg19_chromosome)
                        
                        
positional arguments:
  ontarget_seq          Specify ontarget_seq file
  ref_genome_path       Specify reference genome(e.g. hg19) path

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        prefix
  -s STEP_START, --step_start STEP_START
                        step start (default: 1)
  -e STEP_END, --step_end STEP_END
                        step end (default: 8)

4. Expectation
as an example, if user's linux server has
Intel(R) Xeon(R) CPU E5-4610 v2 @ 2.30GHz (32 cores)
64 Giga byte of memory, this source code require
20hrs ~ 30hrs to finish the prediction




