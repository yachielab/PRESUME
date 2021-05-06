helpdoc="""
usage: PRESUME.py [-h] [--param PARAM] [-V] [--monitor MONITOR] [-n N]
                  [--tree TREE] [-L L] [-f F] [--polyC] [-d D] [-s S] [-T T]
                  [-e E] [--gtrgamma GTRGAMMA] [-u U] [--ud UD] [--ld LD]
                  [-m M] [--constant CONSTANT] [--output OUTPUT] [--qsub]
                  [--idANC IDANC] [--tMorigin TMORIGIN] [--debug] [--bar]
                  [--viewANC] [--save] [--CV] [-r R] [--seed SEED]
                  [--inprob INPROB] [--inlength INLENGTH] [--delprob DELPROB]
                  [--dellength DELLENGTH] [--indels INDELS] [--dop DOP]
                  [--dist DIST] [--editprofile EDITPROFILE]

PRESUME.py

optional arguments:
  -h, --help            show this help message and exit
  --param PARAM         load argument file(csv file)
  -V, --version
  --monitor MONITOR     time limit (default=None)
  -n N                  required number of sequences: if you specified this
                        parameter, timelimit will be postponed until the
                        number of the sequence reach the specified number
                        (default=1)
  --tree TREE           file name of a guide tree in Newick format.
  -L L                  length of sequence (default=1000)
  -f F                  fasta file name　of the common ancestor sequence.
                        (default: poly-C)
  --polyC               use polyC sequence as root
  -d D                  doubling time of origin sequence (default=1)
  -s S                  sigma of doubling time of origin sequence (default=0)
  -T T                  Threashold of doubling time to be deleted (default =
                        1000)
  -e E                  random deletion probability (default=0)
  --gtrgamma GTRGAMMA   parameters for substitution rate matrix
                        GTR{A-C/A-G/A-T/C-G/C-T/G-T} +FU{piA/piC/piG/piT}
                        +G4{shape of gamma distribution} Or, you can use
                        default parameters by "--gtrgamma default" default:
                        GTR{0.3333/0.3333/0.3333/0.3333/0.3333/0.3333}
                        +FU{0.25/0.25/0.25/0.25} +G4{10000}
  -u U                  upper limit of number of sequences (default=2^20)
  --ud UD               upper limit of doubling time (default=10^10)
  --ld LD               lower limit of doubling time (default=10^(-5))
  -m M                  mean of relative subtitution rate according to gamma
                        distribution (default=1)
  --constant CONSTANT   fixed mutation rate of each site (default=None)
  --output OUTPUT       output folder (default:current directory)
  --qsub                activate preparation for distributed processes
                        (defalt=inactivated)
  --idANC IDANC         corresponging ancestral sequence (in upstream tree),
                        in case of distributed computing (default=None)
  --tMorigin TMORIGIN   birth time of origin sequence
  --debug               inactivate deletion of intermediate files
  --bar                 deactivate unstable functions
  --viewANC             generate fasta of ancestoral sequences
  --save                generate args.csv
  --CV                  sigma use as CV(Coefficient Variance) of Normal
                        Distribution
  -r R                  limit of retrying simulation (default=100000)
  --seed SEED           random seed used to initialize the pseudo-random
                        number generator
  --inprob INPROB       file name of insertion probability for each position
  --inlength INLENGTH   file name of insertion length distribution
  --delprob DELPROB     file name of deletion probability for each position
  --dellength DELLENGTH
                        file name of insertion length distribution
  --indels INDELS       file name of indels accumulated before simualtion (for
                        distributed computing mode)
  --dop DOP             Option of qsub for downstream simulation (for
                        distributed computing mode)
  --dist DIST           Distribution of d or 1/d (permissive values: 'norm',
                        'lognorm', 'gamma', 'gamma2') (default: 'gamma2)
  --editprofile EDITPROFILE
                        file name of a base editing profile.
"""
def help_description():
  return helpdoc

if __name__=="__main__":
  pass
