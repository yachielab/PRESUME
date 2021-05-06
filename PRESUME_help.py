helpdoc="""
Usage:
    PRESUME.py 
    [-v] [--version] [-h] [--help] [-n sequence_number] [-L sequence_length]
    [-s standard_deviation] [-e extinction_probability] 
    [--gtrgamma model_parameters] [-m mean_substitution_rate]
    [--constant substitution_probability] [--qsub] [--output directory_path]
    [-f input_file] [--load file_name] [-u sequences_number] [--debug] [--bar] [--save] [-r max_retrial_number] 
    [--seed random_seed] [--limit time_limit] 

Options:
    -v --version
      Print PRESUME version. Ignore all of the other parameters
    -h --help
      Print the usage of PRESUME. Ignore all of the other parameters
    --output <String>
　　　 Output directory path. PRESUME creates a directory unless exists. Default: current directory
    -n <Integer>
      Number of sequences to be generated. Default: 100
    --tree <String>
　　　 Input Newick format file path if a template tree is given. Ignore -s -e -f -u -r --constant –-qsub –-load 
        –-debug –-bar –-save –-seed and --limit
    -s <Float>
      Standard deviation of propagation speed. Default: 0
    -e <Float>
      Probability of extinction. Default: 0
    -f <String>
　　　 Input FASTA file path for the root sequence. Random sequence will be generated unless specified. Ignore -L 
    -L <Integer>
      Length of sequences to be generated. Default: 1000 
    --gtrgamma <String>
      GTR-Gamma model parameters
        Format： --gtrgamma GTR{A-C/A-G/A-T/C-G/C-T/G-T}+FU{piA/piC/piG/piT}+G{alpha}
        For more details, see https://github.com/yachielab/PRESUME/blob/master/SubstitutionModelDetails.PRESUME.pdf
        or you can use the default parameter set by 
          --gtrgamma default
        which is equivalent to
          --gtrgamma GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}+FU{0.25/0.25/0.25/0.25}+G4{10000}
    -m <Float>
      Mean of gamma distribution for relative substitution rates of different sequence positions. Default: 1
    --inprob <String>       
      File path for insertion probability at each sequence position per generation
        (Default: Indels are inactivated)
    --inlength <String>
      File path for insertion length distribution
    --delprob <String>
      File path for deletion probability at each sequence position per generation
    --dellength <String>
      File path for deletion length distribution
    --editprofile <String>
      File path for substitution probabilities at each sequence position per generation
    --qsub
　　　 Execute the distributed computing mode. PRESUME is executed with only a single node unless specified
    --dop <String>
      Option of qsub command when distributed computing mode is enabled with –qsub. Default: ""
    -u <Integer>
　　　 Maximum number of sequences to be generated. Default: 1000000000
    --ud <float>
      Upper limit of doubling time. Default: 10^10
    --ld <float>
      Lower limit of doubling time. Default: 10^(-5)
    --bar
　　　 Activate the monitoring of simulation progress with Python tqdm module
    --save
　　　 Output a CSV file for parameter values used for the simulation. Default: False
    --param <String>
　　　 CSV file for parameter values. This file can be obtained from a previous simulation run executed with
        –-save option.
    --seed <Integer>
　　　 Seed value for generation of random values. Default: 0
    -r <Integer>
　　　 Maximum number of retrials of simulation when all sequences are extinct
      (Default: 100000)    
"""
def help_description():
  return helpdoc

if __name__=="__main__":
  pass
