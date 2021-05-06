<h2>PRESUME Installation and User Manual</h2>

- [Overview of PRESUME](#overview-of-presume)
- [Supported Environment](#supported-environment)
- [Software Dependency](#software-dependency)
- [Software installation](#software-installation)
- [Sample Codes](#sample-codes)
- [PRESUME Usage](#presume-usage)
### Overview of PRESUME

**PRESUME** enables the simulation of various sequence diversification processes based on user-defined models (Figure 1). It is executed with the input information of a root sequence, a target number of sequences to be generated *n*, as well as a parameter *σ* that determines the branch unbalancedness of the generating tree, and parameters that determine how to introduce mutations in sequences. A template tree with defined branch lengths (or generation times) can also be given, such that the mutational process is simulated along with the provided tree.

<img src=images/github_PRESUME.jpg width=60%>

**Figure 1. Schematic diagram of PRESUME.**
<span style="font-size: 60%; color: black;">**Unless a tree of defined topology is provided by the user, PRESUME first simulates a template tree topology for the proliferation of propagating units (PUs) (i.e., sequences or cells as units to harbor multiple sequences) using a single parameter *σ*. Following a linear time *T* progression, PRESUME progressively proliferates PUs, in which each PU duplicates when its generation time *d* has passed since its birth. The generation time *d* is given to each new PU, so its reciprocal (doubling speed) follows a gamma probability distribution, wherein the mean and variance are 1 and *σ*<sup>2</sup>, respectively.**</span>

Upon a template tree is provided by the simulation or by the user, the mutational processes of sequences are simulated along with it. In the present implementation of PRESUME, substitutions are simulated either with time-dependent probability functions assigned for different sequence positions using the GTR-Gamma model or time-independent probabilities per branch (or generation) defined for independent sequence positions by the user.

### Supported Environment

1.	PRESUME can be executed on MacOS or Linux.
2.	The distributed computing mode of PRESUME requires UGE (Univa Grid Engine)

### Software Dependency

1. 1.	Python3 (newer than 3.7.0) with Biopython module *required* and tqdm module *optional; for visualization of simulation progress*

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

Generation of approximately 100 sequences using the GTR-Gamma model with the default parameter set without distributed computing. This computing will take several minutes.

```
PRESUME.py -n 100 --gtrgamma default --save
```

Output: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/tree/master/example/example_1/PRESUMEout) containing the following files will be created in your working directory:

1. `PRESUMEout.fa` : FASTA file for generated descendant sequences

2. `root.fa` : FASTA file describing the root sequence used for the simulation

3. `PRESUMEout.nwk`: Newick format file for the lineage history of the generated sequences

4. `args.csv`: CSV file containing basic patameters used for the simulation (enabled by --save). 

**Example 2**

Generation of approximately 100 sequences using a time-independent model with the substitution frequency of 5% per nucleotide position per generation along with a highly unbalanced lineage trajectory (*σ* of 10). This will take several minutes.

```
PRESUME.py -n 100 --constant 0.05 -s 10
```

Output data: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/tree/master/example/example_2/PRESUMEout) will be created in your working directory.

 **Example 3**

Generation of approximately 10,000 sequences using the GTR-Gamma model with a user-defined parameter set with distributed computing. This will take several minutes.

```
PRESUME.py -n 10000 --gtrgamma GTR{0.927000/2.219783/1.575175/0.861651/4.748809/1.000000}+FU{0.298/0.215/0.304/0.183}+G{0.553549} --qsub 
```

Output data: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/tree/master/example/example_3/PRESUMEout) will be created in your working directory.

**Example 4**

Generation of approximately 128 sequences using GTR-Gamma model and an indel model with user-defined parameters with distributed computing. This will take approximately 1 minute.

```shell
PRESUME.py -n 100 --gtrgamma default --inprob prob.txt --inlength length.txt --delprob prob.txt --dellength length.txt
```

Input data: [`prob.txt`](https://github.com/yachielab/PRESUME/example/example_4/example_4/prob.txt) defines indel probability per generation for each initial sequence postion and [`length.txt`](https://github.com/yachielab/PRESUME/example/example_4/length.txt) defines the distribution of each indel 

Output data: a directory [`PRESUMEout`](https://github.com/yachielab/PRESUME/example/example_4/PRESUMEout) will be created in your working directory.

**Example 5**

Generation of approximately 128 sequences using a time-independent model with the substitution frequency of each nucleotide character in each sequence position specified in transition_probability.txt. This will take approximately 1 minute.

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


```

### Contact

1. Keito Watano (The University of Tokyo) [watano.k10.yachielab@gmail.com](watano.k10.yachielab@gmail.com)
2. Naoki Konno (The University of Tokyo) [naoki@bs.s.u-tokyo.ac.jp](mailto:naoki@bs.s.u-tokyo.ac.jp)
3. Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)

