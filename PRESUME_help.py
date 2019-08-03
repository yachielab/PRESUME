helpdoc="""
Usage:
    PRESUME.py 
    [-v] [--version] [-h] [--help] [-n sequence_number] [-L sequence_length] [-s standard_deviation]
    [-e extinction_probability] [--gtrgamma model_parameters] [-m mean_substitution_rate]
    [--timeind substitution_probability] [--qsub] [--output directory_path] [-f input_file]
    [--load file_name] [-u sequences_number] [--debug] [--bar] [--save] 
    [-r max_retrial_number] [--seed random_seed] [--limit time_limit] 

Options:
    -v --version
      Print PRESUME version; ignore all the other parameters
    -h --help
      Print the usage of PRESUME; ignore all the other parameters
    -n <Integer>
      Number of sequences to be generated. Default: 100
    -L <Integer>
      Length of sequences to be generated. Default: 1000
    -s <Float>
      Standard deviation of propagation speed. Default: 0
    -e <Float>
      Probability of extinction. Default: 0
    --gtrgamma <String>
      Execute GTR-Gamma model to simulate sequence diversification with its parameters
        Parameter format： --gtrgamma GTR{A-C/A-G/A-T/C-G/C-T/G-T}+FU{piA/piC/piG/piT}+G{alpha}
        For more details, see https://github.com/yachielab/PRESUME/blob/master/SubstitutionModelDetails.PRESUME.pdf
        or you can use the default parameter set by 
          --gtrgamma default
        which is equivalent to
          --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{10000}
    -m <Float>
        Mean of gamma distribution for relative substitution rates of different sequence 
          Positions. Default: 1
    --constant <Float>
　　　　　　Execute time-independent model to simulate sequence diversification with a parameter of
          constant substitution probability per generation of every sequence position
    --qsub
        Execute the distributed computing mode
    --output <String>
        Output directory path. PRESUME creates a directory unless exists. Default: current 
          directory
    -f <String>
        Input FASTA file name　for the root sequence. Random sequence will be generated 
          unless specified
    --param <String>
        Basic parameter values of PRESUME can be input as a csv file (the format of the file is same 
          for which is output by PRESUME when -–save is specified)
    -u <Integer>
        Maximum number of sequences to be generated. Default: 1000000000
    --debug
        Output intermediate files
    --bar
        Activate the monitoring of simulation progress with Python tqdm module
    --save
        Output basic parameter values used
    -r <Integer>
        Maximum number of retrials of simulation when all sequences are extinct
    --seed <Integer>
        Seed number for generation of random values. Default: 0
    --monitor <float>
　　　　　　Monitoring paramater. Default: 1
"""

def help_description():
  return helpdoc

if __name__=="__main__":
  pass
