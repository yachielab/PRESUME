<h2>PRESUME Installation and User Manual</h2>

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

See [SubstitutionModelDetails.PRESUME.pdf](https://github.com/yachielab/PRESUME/blob/master/SubstitutionModelDetails.PRESUME.pdf) for more details of how to specify the GTR-Gamma model parameters.

***

Note that as the number of sequences are only sporadically monitored during the simulation, the number of generated descendant sequences can be fluctuated and differed from the number of sequences *N* required to be generated by -n.

In the distributed computing mode, the number of jobs will be around √N; PRESUME first generates ~√N number of sequences in a single node, each of which is then subjected to the further downstream process in a distributed computing node.


#### PRESUME Usage

```
Usage:
    PRESUME.py 
    [-v] [--version] [-h] [--help] [-n sequence_number] [-L sequence_length] [-s standard_deviation]
    [-e extinction_probability] [--gtrgamma model_parameters] [-m mean_substitution_rate]
    [--constant substitution_probability] [--qsub] [--output directory_path] [-f input_file]
    [--load file_name] [-u sequences_number] [--debug] [--bar] [--save] 
    [-r max_retrial_number] [--seed random_seed] [--limit time_limit] 

Options:
    -v --version
      Print PRESUME version; ignore all of the other parameters
    -h --help
      Print the usage of PRESUME; ignore all of the other parameters
    -n <Integer>
      Number of sequences to be generated. Default: 100
    -L <Integer>
      Length of sequences to be generated. Default: 1000
    -s <Float>
      Standard deviation of propagation speed. Default: 0
    -e <Float>
      Probability of extinction. Default: 0
    --gtrgamma <String>
      GTR-Gamma model parameters
        Format： --gtrgamma GTR{A-C/A-G/A-T/C-G/C-T/G-T}+FU{piA/piC/piG/piT}+G{alpha}
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
　　　　　 Output directory path. PRESUME creates a directory unless exists. Default: current directory
    -f <String>
　　　　　 Input FASTA file name for the root sequence. Random sequence will be generated unless specified
    -u <Integer>
　　　　　 Maximum number of sequences to be generated. Default: 1000000000
    --debug
　　　　　 Output intermediate files
    --bar
　　　　　 Activate the monitoring of simulation progress with Python tqdm module
    --save
　　　　　 Output a CSV file for parameter values used for the simulation
    --param <String>
　　　　　 CSV file for parameter values.
　　　　　 This file can be obtained from a previous simulation run executed with –-save option.
    -r <Integer>
　　　　　 Maximum number of retrials of simulation when all sequences are extinct
    --seed <Integer>
　　　　　 Seed value for generation of random values. Default: 0
    --monitor <float>
　　　　　 Stepper size parameter for monitoring of lineage generation. Default: 1
    --tree <String>
　　　　　 Input Newick format file name if a template tree is given.
　　　　　   The following parameters will be ignored:
　　　　　     -L -s -e -f -u -r --constant –-qsub –-load –-debug –-bar –-save –-seed --limit

```

### Contact

1. Keito Watano (The University of Tokyo) watano.k10.yachielab@gmail.com
2. Naoki Konno (The University of Tokyo) [naoki@bs.s.u-tokyo.ac.jp](mailto:naoki@bs.s.u-tokyo.ac.jp)
3. Nozomu Yachie (The University of Tokyo) [nzmyachie@gmail.com](mailto:yachie@synbiol.rcast.u-tokyo.ac.jp)

