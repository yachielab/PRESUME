<h2>PRESUME Installation and User Manual</h2>

- [Overview of PRESUME](#overview-of-presume)
- [Supported Environment](#supported-environment)
- [Software Dependency](#software-dependency)
- [Software installation](#software-installation)
- [Sample Codes](#sample-codes)
- [PRESUME Usage](#presume-usage)
### Overview of PRESUME

**PRESUME** is a software tool that simulates cell division or speciation and diversification of DNA sequences in the growing population. By employing a distributed computing platform, PRESUME progressively
generates a large number of sequences that accumulate substitutions with their lineage history information. In this process, daughter sequences are duplicated at a certain speed which is incompletely inherited from that of the maternal sequence under a stochastic model (Figure 1). The substitution probability at different positions in each sequence is defined in a time-dependent manner using GTR-Gamma model or set to a certain rate. The software allows the user to simulate various types of sequence diversification processes with different sets of input parameters.

<img src=images/presume_concept.jpg width=50%>

**Figure 1.** Schematic diagram of PRESUME. PRESUME simulates the propagation and diversification of sequences that accumulate substitutions and generates a large set of descendant sequences with lineage information. *m* refers to a maternal sequence, and *d<sub>1</sub>* and *d<sub>2</sub>* refers to two daughter sequences derived from *m*. In this simulation, the doubling times of the two daughter sequences (*t<sub>d1</sub>* and *t<sub>d2</sub>*) are incompletely inherited from the doubling time of the mother sequence (*t<sub>m</sub>*). This occurs under a stochastic model, in which 1/*t<sub>d1</sub>* and 1/*t<sub>d2</sub>* follow a normal distribution where the mean and variance are 1/*t<sub>m</sub>* and *&sigma;*<sup>2</sup> respectively. Additionally, sequence extinction is set at a random rate (&epsilon;) and also occurs when the sequence doubling speed reaches a negative value. The substitution probabilities at different positions in each sequence of length *L* are defined in a time-dependent manner using GTR-Gamma model with parameters *Q*, *&alpha;* and *&mu;*, or set to a certain rate *&phi;* (see [SubstitutionModelDetails.PRESUME.pdf](https://github.com/yachielab/PRESUME/blob/master/SubstitutionModelDetails.PRESUME.pdf)).

### Supported Environment

1. PRESUME can be executed on MacOS or Linux.
2. The distributed computing mode of PRESUME requires UGE (Univa Grid Engine) 

### Software Dependency

1. Python3 (newer than 3.7.0) with Biopython module *required* and tqdm module *optional; if you want to visualize a simulation progress*

### Software Installation

##### Installation of PRESUME
Each step of installation takes less than 1 min.

1. Download PRESUME by

   ```
   git clone https://github.com/yachielab/PRESUME
   ```

2. Add the absolute path of PRESUME directory to $PATH

3. Make PRESUME executable

   ```
   chmod u+x PRESUME.py
   ```

##### Installation of [Anaconda](https://www.anaconda.com/distribution/) (required)

1. Execute the following commands

   ```
   wget https://repo.anaconda.com/archive/Anaconda3-2018.12-Linux-x86_64.sh
   bash Anaconda3-2018.12-Linux-x86_64.sh
   ```

2. Set $PATH to anaconda3/bin

##### Installation of [Biopython](https://anaconda.org/anaconda/biopython) 1.76 (required)

1. Install Biopython by

   ```shell
   conda install -c anaconda biopython
   ```

##### Installation of [tqdm](https://anaconda.org/conda-forge/tqdm) 4.43.0 (optional)

1. Install tqdm by

   ```
   conda install -c conda-forge tqdm
   ```

### Sample Codes

The software functions can be tested by the following example commands:

**Example 1**

Generation of ~100 sequences using GTR-Gamma model with the default parameter set without distributed computing. The computation will take several minutes.

```
PRESUME.py -n 100 --gtrgamma default --save
```

Output: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/tree/master/example/example_1/PRESUMEout) containing the following files will be created in your working directory:

1. `PRESUMEout.fa` : FASTA file for generated descendant sequences

2. `root.fa` : FASTA file describing the root sequence used for the simulation

3. `PRESUMEout.nwk`: Newick format file for the lineage history of the generated sequences

4. `args.csv`: CSV file containing basic patameters used for the simulation (enabled by --save). 

**Example 2**

Generation of ~100 sequences using a time-independent model with the substitution frequency of 5% per site per generation along with a highly unbalanced lineage trajectory (*&sigma;* of 10). The computation will take several minutes.

```
PRESUME.py -n 100 --constant 0.05 -s 10
```

Output data: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/tree/master/example/example_2/PRESUMEout) will be created in your working directory.

 **Example 3**

Generation of ~10,000 sequences using GTR-Gamma model with a defined parameter set with distributed computing. The computation will take several minutes.

```
PRESUME.py -n 10000 --gtrgamma GTR{0.927000/2.219783/1.575175/0.861651/4.748809/1.000000}+FU{0.298/0.215/0.304/0.183}+G{0.553549} --qsub 
```

Output data: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/tree/master/example/example_3/PRESUMEout) will be created in your working directory.

**Example 4**

Generation of ~100 sequences using GTR-Gamma model and an original indel model with a defined parameter set with distributed computing. The computation will take ~1 minute.

```shell
PRESUME.py -n 100 --gtrgamma default --inprob prob.txt --inlength length.txt --delprob prob.txt --dellength length.txt
```

Input data: [`prob.txt`](https://github.com/yachielab/PRESUME/example/example_4/example_4/prob.txt) defines indel probability per generation for each initial sequence postion and [`length.txt`](https://github.com/yachielab/PRESUME/example/example_4/length.txt) defines the distribution of each indel 

Output data: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/example/example_4/PRESUMEout) will be created in your working directory.

**Example 5**

Generation of ~100 sequences using time-independent model with the substitution frequency of each site specified in `editprofile.txt`. The computation will take ~1 minute.

```shell
PRESUME.py -n 100  --editprofile editprofile.txt
```
Input data: [`editprofile.txt`]() defines substitution matrix of each site.

Output data: a directory [`PRESUMEout`]() will be created in your working directory.

See [SubstitutionModelDetails.PRESUME.pdf](https://github.com/yachielab/PRESUME/blob/master/SubstitutionModelDetails.PRESUME.pdf) for more details of how to specify the GTR-Gamma model parameters.

***

Note that as the number of sequences are only sporadically monitored during the simulation, the number of generated descendant sequences can be fluctuated and differed from the number of sequences *N* required to be generated by -n.

In the distributed computing mode, the number of jobs will be around √N; PRESUME first generates ~√N number of sequences in a single node, each of which is then subjected to the further downstream process in a distributed computing node.


#### PRESUME Usage

```
usage: PRESUME.py [-h] [--param PARAM] [-V] [--monitor MONITOR] [-n N]
                  [--tree TREE] [-L L] [-f F] [--polyC] [-d D] [-s S] [-T T]
                  [-e E] [--gtrgamma GTRGAMMA] [-u U] [--ud UD] [--ld LD]
                  [-m M] [--constant CONSTANT] [--output OUTPUT] [--qsub]
                  [--idANC IDANC] [--tMorigin TMORIGIN] [--debug] [--bar]
                  [--viewANC] [--save] [--CV] [-r R] [--seed SEED]
                  [--inprob INPROB] [--inlength INLENGTH] [--delprob DELPROB]
                  [--dellength DELLENGTH] [--indels INDELS] [--chunks CHUNKS]
                  [--dop DOP] [--dist DIST] [--editprofile EDITPROFILE]

PRESUME.py

optional arguments:
  -h, --help            show this help message and exit
  --param PARAM         load argument file(csv file)
  -V, --version
  --monitor MONITOR     time limit (default=None)
  -n N                  required number of sequences: if you specified this
                        parameter, timelimit will be postponed until the
                        number of the sequence reach the specified number
                        (default=None)
  --tree TREE           file name of a guide tree in Newick format(for debug,
                        unsupported.)
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
                        file name of a base editing profile(for debug,
                        unsupported.)

```

### Contact

1. Keito Watano (The University of Tokyo) [watano.k10.yachielab@gmail.com](watano.k10.yachielab@gmail.com)
2. Naoki Konno (The University of Tokyo) [naoki@bs.s.u-tokyo.ac.jp](mailto:naoki@bs.s.u-tokyo.ac.jp)
3. Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)

