helpdoc="""
   shellv$ python3 PRESUME.py -h
   usage: PRESUME.py [-h] [--load FileName] [-V] [--limit TimeLimit] [-n NumberOfSequences] 
   		              [-L SequenceLength] [-f FASTAfile] [-d DoublingTime] 
   		              [-s StandardDeviation] [-e DeletionProbability] 
   		              [--model ModelParameter] [-u NumberOfSequences] [-m SubstitutionRate]
   		              [--delta SubstitutionFrequency] [--qsub] [--debug] [--bar] [--save]
                     [--no_retry]
   
   PRESUME.py
    
    optional arguments:
    -h, --help           shows this help message and exit.
    --load FILENAME			 loads argument from the specified file.
    										 (The format is same as the csv file PRESUME outputs when --option attached.)
    -V, --version        shows version of PRESUME.
    --limit LIMIT        specifies the time limit (default=1)
    -n N                 specifies required number of sequences: when you specified this
                         option, the time limit becomes the time point when N                   
                         sequences are present for the first time. (default=None)
    -L L                 specifies the length of a sequence (default=1000)
    -f F                 specifies the fasta file name　of the initial sequence.
                         (default: The initial sequence will be randomly generated.)
    -d D                 specifies the doubling time of the initial sequence (default=1)
    -s S                 specifies the standard deviation of doubling time of the initial      
                         sequence (default=0)
    -e E                 specifies the random deletion probability: any sequence is
                         randomly deleted at the probability of E when it 
    										 emerges (default=0)
    --output OUTPUT      specifies the output folder (default:current working directory)
    --model MODEL        specifies parameters for substitution rate matrix 
    										 Format： --model GTR{A-C/A-G/A-T/C-G/C-T/G-T}+FU{piA/piC/piG/piT}+G{shape of gamma distribution}
    										 For details, please also refer another document on substitution model (<link>)
                         Or, you can use default parameters by 
                         --model default
                         (default: GTR{GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333}
                         +FU{0.25/0.25/0.25/0.25}+G4{10000} )
    -u U                 specifies the upper limit of number of sequences for safety (default=2^20)
    -m M                 specifies the mean of relative substitution rate following gamma distribution (default=1)
    --delta DELTA        specifies the fixed substitution frequency of every site in time-independent model (default=None)
    --qsub               activates the distributed computing mode (defalt=inactivated)
    --debug              inactivates the deletion of intermediate files
    --bar                activates showing progress bar (python tqdm module is required.)
    --save               saves the argument information in args.csv file
    --no_retry           inactivates retrying simulation when all sequences are deleted.
    --seed							 spacifies the seed of random number
    --polyC              set the initial sequence as poly-C.
    -r R								 specifies maximum iteration number:
                         If all sequences are deleted, PRESUME.py repeats the simulation from the beginning at most R times.
"""

def help_description():
  return helpdoc

if __name__=="__main__":
  pass