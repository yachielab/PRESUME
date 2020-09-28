#!/usr/bin/env python
# -*- coding: utf-8 -*-
__version__ = "2.0.0"

__author__ = "\
            Keito Watano <watano.k10.yachielab@gmail.com>,\
            Naoki Konno <naoki@bs.s.u-tokyo.ac.jp>, \
            Nozomu Yachie <nzmyachie@gmail.com>"

__date__ = "2020/9/4"

import sys
import os
from Bio import Phylo, SeqIO
import subprocess
import numpy as np
import random
import argparse
import gzip

from submodule import nwk2fa_mutation
from submodule import args_reader

LOGO='''
######     ######     #######     #####     #     #    #     #    #######
#     #    #     #    #          #     #    #     #    ##   ##    #
#     #    #     #    #          #          #     #    # # # #    #
######     ######     #####       #####     #     #    #  #  #    #####
#          #   #      #                #    #     #    #     #    #
#          #    #     #          #     #    #     #    #     #    #
#          #     #    #######     #####      #####     #     #    #######


###### #####   ####  #    #    #    # ###### #    # #  ####  #    # 
#      #    # #    # ##  ##    ##   # #      #    # # #    # #   #  
#####  #    # #    # # ## #    # #  # #####  #    # # #      ####   
#      #####  #    # #    #    #  # # #      # ## # # #      #  #   
#      #   #  #    # #    #    #   ## #      ##  ## # #    # #   #  
#      #    #  ####  #    #    #    # ###### #    # #  ####  #    # 
                                                                    
                                                  
#####  ####     ######   ##    ####  #####   ##   
  #   #    #    #       #  #  #        #    #  #  
  #   #    #    #####  #    #  ####    #   #    # 
  #   #    #    #      ######      #   #   ###### 
  #   #    #    #      #    # #    #   #   #    # 
  #    ####     #      #    #  ####    #   #    # 

DEBUG MODE of PRESUME
from newick to fasta
github: https://github.com/yachielab/PRESUME
'''

def tabulate_names(tree):
    ### from https://biopython.org/wiki/Phylo_cookbook
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name is None:
            clade.name = str(idx)
    return tree

def topology_shaper(tree):
    internal_nodes = list(tree.find_clades(terminal=False, order='preorder'))
    terminal_nodes = list(tree.find_clades(terminal=True, order='preorder'))
    topology = {} # <mother>:[<daughter_L>, <daughter_R>]
    branch_length_dict = {} # <mother.name> : <branch_length>
    root = tree.root
    branch_length_dict[root.name] = root.branch_length if root.branch_length is not None else 1.0
    for clade in internal_nodes:
        if len(clade.clades) == 2:
            mother_name     = clade.name
            daughter_L_name = clade.clades[0].name
            daughter_R_name = clade.clades[1].name
            mother_branch_length  = clade.branch_length if clade.branch_length is not None else 1.0
            # define topology
            topology[mother_name] = [daughter_L_name, daughter_R_name]
            #topology[daughter_L_name] = []
            #topology[daughter_R_name] = []
            # define branch_length
            branch_length_dict[mother_name] = clade.branch_length
        else:
            raise ValueError('wrong clades! \n # of clades {}'.format(len(clade.clades)))
            pass

    for clade in terminal_nodes:
        topology[clade.name] = []
        if clade.name != root.name:
            branch_length_dict[clade.name] = clade.branch_length if clade.branch_length is not None else 1.0
    return topology, branch_length_dict, root.name

# Conduct a simulation
def translate_tree(topology_dict, branch_length_dict,name_of_root, initseq, parsed_args):
    if not initseq:
        initseq=''.join([np.random.choice(['A', 'G', 'C', 'T']) for i in range(1000)])
    init_clade = nwk2fa_mutation.Lineage(branch_length=branch_length_dict[name_of_root], name=name_of_root, seq= initseq, ROOT=True, parsed_args=parsed_args)
    newtree = Phylo.BaseTree.Tree(init_clade)
    stack=[init_clade]
    cnt = 0
    while (stack !=[]):
        clade=stack.pop()
        node_name=clade.name
        mother_seq=clade.seq
        if len(topology_dict[node_name]) ==2:
            children = [
                nwk2fa_mutation.Lineage(branch_length = branch_length_dict[topology_dict[node_name][0]], name=str(topology_dict[node_name][0]), seq=mother_seq, parsed_args=parsed_args),
                nwk2fa_mutation.Lineage(branch_length = branch_length_dict[topology_dict[node_name][1]], name=str(topology_dict[node_name][1]), seq=mother_seq, parsed_args=parsed_args)
            ]



            print(children)
            
            
            
            clade.clades.extend(children)
            stack.extend(children)
            cnt += 1
    return newtree

def nwk2fa_light(tree, initseq=False,parsed_args=None):
    # translate Phylo.BaseTree.Tree into Lineage
    topology_dict, branch_length_dict, name_of_root= topology_shaper(tree)
    print(topology_dict)
    lineage_tree = translate_tree(topology_dict, branch_length_dict, name_of_root, initseq,parsed_args=parsed_args)
    if len(tree.get_terminals()) != len(lineage_tree.get_terminals()):
        raise ValueError('something went wrong!!!')
    fasta_data = {item.name:item.seq for item in lineage_tree.get_terminals()}
    return fasta_data, tree

def fasta_writer_single(name_seq_dict, outfp):
    for name in name_seq_dict.keys():
        filename = "{0}/{1}.fasta".format(outfp, name)
        with open(filename, "w") as writer:
            writer.write(">{0}\n{1}\n".format(name, name_seq_dict[name]))
    return

def fasta_writer_multiple(name_seq_dict, outfp, filename):
    filename = "{0}/{1}.gz".format(outfp, filename)
    with gzip.open(filename, "wt") as writer:
        for name in name_seq_dict.keys():
            writer.write(">{0}\n{1}\n".format(name, name_seq_dict[name]))
    return

# decomp tree
def get_node2Ntips(tree):
    # Caution: all internal nodes and tips must be named!!!!
    node2Ntips = {}
    for tip in tree.get_terminals():
        node2Ntips[tip.name] = 1
    
    idx = 0
    for internal in tree.get_nonterminals():
        if internal.name is None:
            internal.name = "internal_{}".format(int(idx))
            idx += 1
    for internal in reversed(list(tree.get_nonterminals())):
        node2Ntips[internal.name] = sum([node2Ntips[child.name] for child in internal.clades])
        
    return node2Ntips

# func
def decompose(
    treefile,
    number,
    output
    ):
    # All nodes must be named!!!

    tree = Phylo.read(treefile, 'newick')

    # count downstream tips

    node2Ntips = get_node2Ntips(tree)

    # print(node2Ntips)

    # decompose tree

    stack = [tree.clade]
    counter = 1
    while counter != 0:
        node = stack.pop()
        if ( node2Ntips[node.name] > number ):
            stack.extend(node.clades)
            counter += 1
        else:
            newtree    = Phylo.BaseTree.Tree(node)
            Phylo.write(newtree, output + "/Down_"+node.name+".nwk", 'newick')
            node.clades = []
            counter -= 1
    Phylo.write(tree, output + "/deUp.nwk", 'newick')
    return

def shell_generator(shell_outfp, treefile_list, fastafile_list, tree_outfp, fasta_outfp, stdeo):
    PYTHON3 = (((
        subprocess.Popen('which python3', stdout=subprocess.PIPE, shell=True)
        .communicate()[0])
        .decode('utf-8'))
        .split('\n')
        )[0]

    NWK2FA = os.path.dirname(os.path.abspath(__file__))

    PATH = (((
        subprocess.Popen('echo $PATH', stdout=subprocess.PIPE,
                            shell=True)
        .communicate()[0])
        .decode('utf-8'))
        .split('\n'))[0]

    LD_LIBRARY_PATH = (((
        subprocess.Popen('echo $LD_LIBRARY_PATH', stdout=subprocess.PIPE,
                        shell=True)
        .communicate()[0])
        .decode('utf-8'))
        .split('\n'))[0]

    terminal_idx = 0
    for idx, (treefile, fastafile) in enumerate(zip(treefile_list, fastafile_list)):
        with open("{}/downstream_{}.sh".format(shell_outfp, idx + 1), 'w') as qf:
            qf.write("#!/bin/bash\n")
            qf.write("#$ -S /bin/bash\n")
            qf.write("#$ -cwd\n")
            qf.write("PATH={}\n".format(PATH))
            qf.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
            qf.write("pwd\n")

            python_command = PYTHON3 + " " + NWK2FA + "/nwk2fa.py "\
                + " --tree " + treefile\
                + " --fasta " + fastafile\
                + " --outdir_nwk " + tree_outfp\
                + " --outdir_fasta " + fasta_outfp\
                + " --filename " + "Down_{}".format(idx + 1)

            qf.write(python_command)
        terminal_idx = idx
    
    submit_command = "qsub -e {3} -o {3} -sync y -t 1-{0} {1}/nwk2fa_launcher.sh {2} &> /dev/null".format(
        terminal_idx + 1,
        NWK2FA,
        shell_outfp,
        stdeo
        )
    return submit_command

def nwk2fa_qsub(args, parsed_args):
    INFILE, OUTDIR =args.tree, args.output
    initseq = False
    # initial sequence specification
    if (args.f is not None):
        with gzip.open(args.f, 'rt') as f:
            sequences = SeqIO.parse(f, 'fasta')
            initseq = str(list(sequences)[0].seq)

    # prepare intermediates
    intermediate_path = "{}/intermediate".format(OUTDIR)
    os.makedirs(intermediate_path, exist_ok = True)

    decomp_nwk_path = "{}/decomp_tree".format(intermediate_path)
    os.makedirs(decomp_nwk_path, exist_ok = True)

    # decompose tree
    tree = Phylo.read(INFILE, "newick")
    tips_threshold = int(len(tree.get_terminals())**(1/2))
    decompose(INFILE, tips_threshold, decomp_nwk_path)
    filelist = os.listdir(decomp_nwk_path)
    terminals = []
    for file in filelist:
        if file.split("_")[0] == "Down":
            filename_without_ext = file.split(".")[0]
            terminal_name = "_".join(filename_without_ext.split("_")[1:])
            terminals.append(terminal_name)
    
    upper_tree = Phylo.read("{}/deUp.nwk".format(decomp_nwk_path), "newick")
    upper_terminals = [terminal.name for terminal in upper_tree.get_terminals()]

    # error check
    for terminal in terminals:
        if not (terminal in upper_terminals):
            print(terminal, " is not found!!")
    print("NWK2FA:created upper tree...")
    intermediate_fasta_path = "{}/fastas".format(intermediate_path)
    os.makedirs(intermediate_fasta_path, exist_ok = True)

    upper_fasta, tree = nwk2fa_light(upper_tree, initseq)
    fasta_writer_single(upper_fasta, intermediate_fasta_path)

    fasta_filelist = os.listdir(intermediate_fasta_path)
    nwk_filelist = os.listdir(decomp_nwk_path)

    for file in nwk_filelist:
        if file.split("_")[0] == "Down":
            filename_without_ext = file.split(".")[0]
            nodename = "_".join(filename_without_ext.split("_")[1:])
            if not ("{}.fasta".format(nodename) in fasta_filelist):
                print(nodename,".fasta not found!")

    downstream_fasta_path = "{}/down_fastas".format(intermediate_path)
    downstream_newick_path = "{}/down_newicks".format(intermediate_path)
    os.makedirs(downstream_fasta_path, exist_ok = True)
    os.makedirs(downstream_newick_path, exist_ok = True)

    treefile_list, fastafile_list = [], []
    for file_nwk in nwk_filelist:
        if file_nwk.split("_")[0] == "Down":
            filename_without_ext = file_nwk.split(".")[0]
            nodename = "_".join(filename_without_ext.split("_")[1:])
            file_fasta = "{}.fasta".format(nodename)

            path_of_fasta_in = "{}/{}".format(intermediate_fasta_path, file_fasta)
            path_of_newick_in = "{}/{}".format(decomp_nwk_path, file_nwk)

            fastafile_list.append(path_of_fasta_in)
            treefile_list.append(path_of_newick_in)

    shell_path = "{}/shell".format(intermediate_path)
    os.makedirs(shell_path, exist_ok = True)
    stdeo_path = "{}/stdeo".format(intermediate_path)
    os.makedirs(stdeo_path, exist_ok = True)    
    submit_command = shell_generator(shell_path, treefile_list, fastafile_list, downstream_newick_path, downstream_fasta_path ,stdeo_path)
    subprocess.call(submit_command, shell=True)
    print("bottom tree created!")

    command = "cat {}/* > {}/PRESUMEout.fasta.gz".format(downstream_fasta_path, OUTDIR)
    subprocess.call(command, shell=True)
    print("Done!")
    return

def nwk2fa_single(args, parsed_args):
    initseq = False
    # initial sequence specification
    if (args.f is not None):
        with gzip.open(args.f, 'rt') as f:
            sequences = SeqIO.parse(f, 'fasta')
            initseq = str(list(sequences)[0].seq)
    
    # read tree
    tree = Phylo.read(args.tree, "newick")

    # translate tree
    newtree_fasta, newtree = nwk2fa_light(tree, initseq, parsed_args=parsed_args)
    fasta_writer_multiple(newtree_fasta, args.outdir_fasta, "{}.fasta".format(args.filename))
    Phylo.write(newtree, "{}/{}.nwk".format(args.outdir_nwk, args.filename), "newick")

if __name__ == "__main__":
    # INFILE = "/Users/keitowatano/Desktop/tmp_PRESUME/in/test_10k_tips.nwk"
    # OUTDIR = "/Users/keitowatano/Desktop/tmp_PRESUME/out/prototype_nwk2fa"
    parser = argparse.ArgumentParser(description='PRESUME.py', add_help=True)
    parser.add_argument(
        "--tree",
        type=str
        )
    parser.add_argument(
        "--fasta",
        type=str
        )

    parser.add_argument(
        "--outdir_nwk",
        type=str
        )

    parser.add_argument(
        "--outdir_fasta",
        type=str
        )

    parser.add_argument(
        "--filename",
        type=str
        )

    args = parser.parse_args()

    # read sequence
    with open(args.fasta, "r") as reader:
        for seq_record in SeqIO.parse(reader, "fasta"):
            initseq = str(seq_record.seq)
    
    # read tree
    tree = Phylo.read(args.tree, "newick")

    # translate tree
    newtree_fasta, newtree = nwk2fa_light(tree, initseq)
    fasta_writer_multiple(newtree_fasta, args.outdir_fasta, "{}.fasta".format(args.filename))
    Phylo.write(newtree, "{}/{}.nwk".format(args.outdir_nwk, args.filename), "newick")

