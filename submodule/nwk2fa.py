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
import shutil

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

def rename_internals(tree):
    ### from https://biopython.org/wiki/Phylo_cookbook
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name is None:
            clade.name = "clade"+str(idx)
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
            if clade.branch_length==None: 
                print("Notice!: Input tree contains branches whose lengths are None.")
                clade.branch_length=10^(-10)             # if branch_length is not defined, PRESUME ignore mutations ccuring at the branch 
        else:
            raise ValueError('wrong clades! \n # of clades {}'.format(len(clade.clades)))
            pass

        
    for clade in terminal_nodes:
        topology[clade.name] = []
        if clade.name != root.name:
            branch_length_dict[clade.name] = clade.branch_length if clade.branch_length is not None else 1.0 
    
    return topology, branch_length_dict, root.name

# Conduct a simulation
def translate_tree(topology_dict, branch_length_dict,name_of_root, initseq, parsed_args): # should change func name
    
    #init_branch_length = parsed_args.dorigin
    init_branch_length = branch_length_dict[name_of_root]

    init_clade = nwk2fa_mutation.Lineage(branch_length=init_branch_length, name=name_of_root, seq= initseq, ROOT=True, parsed_args=parsed_args,indelsM=[])
    init_clade.indels=parsed_args.initindels
    newtree = Phylo.BaseTree.Tree(init_clade)
    stack=[init_clade]
    cnt = 0
    while (stack !=[]):
        clade=stack.pop()
        node_name=clade.name
        mother_seq=clade.seq
        if len(topology_dict[node_name]) == 2:
            children = [
                nwk2fa_mutation.Lineage(branch_length = branch_length_dict[topology_dict[node_name][0]], name=str(topology_dict[node_name][0]), seq=mother_seq, parsed_args=parsed_args, mother_clade = clade, indelsM=clade.indels),
                nwk2fa_mutation.Lineage(branch_length = branch_length_dict[topology_dict[node_name][1]], name=str(topology_dict[node_name][1]), seq=mother_seq, parsed_args=parsed_args, mother_clade = clade, indelsM=clade.indels)
            ]
            clade.clades.extend(children)
            stack.extend(children)
            cnt += 1
    return newtree

def nwk2fa_light(tree, initseq,parsed_args):
    # translate Phylo.BaseTree.Tree into Lineage
    topology_dict, branch_length_dict, name_of_root= topology_shaper(tree)
    lineage_tree = translate_tree(topology_dict, branch_length_dict, name_of_root, initseq,parsed_args=parsed_args)

    if len(tree.get_terminals()) != len(lineage_tree.get_terminals()):
        raise ValueError('something went wrong!!!')
    

    name2seq_without_indel  = {terminal.name:terminal.seq for terminal in lineage_tree.get_terminals()}

    if (parsed_args.CRISPR):
        name2seq        = {}
        name2alignedseq = {}
        name2indellist  = {}
        for terminal in lineage_tree.get_terminals():
            seq, alignedseq, indel_list    = terminal.get_seq_with_indel()
            name2seq[terminal.name]        = seq
            name2alignedseq[terminal.name] = alignedseq
            name2indellist[terminal.name]  = indel_list
    else:
        name2seq, name2alignedseq = name2seq_without_indel, name2seq_without_indel
        name2indellist            = None

    return name2seq_without_indel, name2seq, name2alignedseq, name2indellist, tree, lineage_tree

def fasta_writer_single(name_seq_dict, outfp):
    for name in name_seq_dict.keys():
        filename = "{0}/{1}.fasta".format(outfp, name)
        with open(filename, "w") as writer:
            writer.write(">{0}\n{1}\n".format(name, name_seq_dict[name]))
    return

def indel_writer_single(name_indellist_dict, outfp):
    for name in name_indellist_dict.keys():
        filename = "{0}/{1}.indel".format(outfp, name)
        with open(filename, "w") as writer:
            indellist = name_indellist_dict[name]
            for indel in indellist:
                if (indel[0] == "del") : writer.write(str(indel[0])+"\t"+str(indel[1])+"\t"+str(indel[2])+"\n")
                elif (indel[0] == "in"): writer.write(str(indel[0])+"\t"+str(indel[1])+"\t"+str(indel[2])+"\t"+str(indel[3])+"\n")
    return

def fasta_writer_multiple(name_seq_dict, outfp, filename):
    filename = "{0}/{1}.fa.gz".format(outfp, filename)
    with gzip.open(filename, "wt") as writer:
        for name in name_seq_dict.keys():
            writer.write(">{0}\n{1}\n".format(name, name_seq_dict[name]))
    return

def indel_writer_multiple(name_indellist_dict, outfp, filename):
    filename = "{0}/{1}.indel.gz".format(outfp, filename)
    with gzip.open(filename, "wt") as writer:
        for name in name_indellist_dict.keys():
            writer.write("{0}\t".format(name))
            
            for indel_idx, indel in enumerate(name_indellist_dict[name]):
                if (indel[0] == "del") : 
                    mid    = indel[1] # pos is the midpoint of deletion
                    length = indel[2]
                    writer.write(("D_mid"+str(mid)+"_len"+str(length)))
                elif (indel[0] == "in"): 
                    pos    = indel[1] # pos is the midpoint of deletion
                    seq    = indel[3]
                    writer.write(("I_"+str(pos+0.5)+"_"+seq))
                if (indel_idx < len(name_indellist_dict[name])-1):
                    writer.write(";")
            writer.write("\n")
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

def shell_generator(shell_outfp, treefile_list, fastafile_list, indelfile_list, tree_outfp, fasta_outfp, stdeo, args, parsed_args):
    PYTHON3 = (((
        subprocess.Popen('which python3', stdout=subprocess.PIPE, shell=True)
        .communicate()[0])
        .decode('utf-8'))
        .split('\n')
        )[0]

    NWK2FAdir  = os.path.dirname(os.path.abspath(__file__))
    PRESUMEdir = "/".join(NWK2FAdir.split("/")[:-1])

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
    for idx, (treefile, fastafile, indelfile) in enumerate(zip(treefile_list, fastafile_list, indelfile_list)):
        
        parent_name= treefile.split("/")[-1].split(".")[0]
        OUTPUT_DIR = fasta_outfp+"/{}".format(parent_name)

        with open("{}/downstream_{}.sh".format(shell_outfp, idx+1), 'w') as qf:
            qf.write("#!/bin/bash\n")
            qf.write("#$ -S /bin/bash\n")
            qf.write("#$ -cwd\n")
            qf.write("PATH={}\n".format(PATH))
            qf.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
            qf.write("mkdir "+OUTPUT_DIR+"\n")
            qf.write("cd "+OUTPUT_DIR+"\n")
            qf.write("pwd\n")

            python_command = PYTHON3 + " " + PRESUMEdir + "/PRESUME.py "\
                + " -f " + fastafile\
                + " --seed " + str(np.random.randint(0, 10000))\
                + " --tree " + treefile
            if parsed_args.CRISPR:
                python_command += \
                    " --inprob "     + args.inprob   +\
                    " --inlength "   + args.inlength +\
                    " --delprob "    + args.delprob  +\
                    " --dellength "  + args.dellength+\
                    " --indels "     + indelfile
            if (args.gtrgamma is not None):
                python_command += \
                    " --gtrgamma "+str(args.gtrgamma)
            if (args.constant is not None):
                python_command += \
                    " --constant "+str(args.constant)
            if (args.editprofile is not None):
                python_command += \
                    " --editprofile "+str(args.editprofile)
            if (args.debug):
                python_command += \
                    " --debug "

            qf.write(python_command)
        terminal_idx = idx
    
    #submit_command = "qsub -e {3} -o {3} -sync y -t 1-{0} {1}/nwk2fa_launcher.sh {2} &> /dev/null".format(
    submit_command = "qsub -e {3} -o {3} -sync y -t 1-{0} {1}/submodule/nwk2fa_launcher.sh {2} > {3}/nwk2fa_launcher.sh.out 2> {3}/nwk2fa_launcher.sh.err".format(
        terminal_idx + 1,
        PRESUMEdir,
        shell_outfp,
        stdeo
        )
    return submit_command

def count_mutations_per_branch(lineage, outfp):

    def seq_dist(str1, str2):
        if (len(str1)!=len(str2)):
            print("Error: seq_dist() len(str1)!=len(str2)")
            return None
        else:
            Ndiff = 0
            for i in range(len(str1)):
                if str1[i] != str2[i]:
                    Ndiff += 1
            return Ndiff

    for node in lineage.get_nonterminals():
        for child in node.clades:
            print(node.name, child.name, seq_dist(node.seq, child.seq), sep = ',', file = outfp)

def count_mutations_per_position(lineage, outfp):

    def seq_diff(str1, str2):
        if (len(str1)!=len(str2)):
            print("Error: seq_dist() len(str1)!=len(str2)")
            return None
        else:
            diffvec = np.zeros(len(str1))
            for i in range(len(str1)):
                if str1[i] != str2[i]:
                    diffvec[i] = 1
            return diffvec

    L = len(lineage.clade.seq)
    count_array = np.zeros(L)

    for node in lineage.get_nonterminals():
        for child in node.clades:
            count_array += seq_diff(node.seq, child.seq)
    
    for i in range(L):
        print(str(i),count_array[i], sep = ',', file = outfp)

def nwk2fa_qsub(args, parsed_args):
    INFILE, OUTDIR =args.tree, args.output
    # initial sequence specification
    initseq = parsed_args.initseq
    fasta_writer_multiple({"root":initseq}, args.output+"/PRESUMEout", "root")

    '''
    initseq = False
    # initial sequence specification
    if (args.f is not None):
        with gzip.open(args.f, 'rt') as f:
            sequences = SeqIO.parse(f, 'fasta')
            initseq = str(list(sequences)[0].seq)
    '''

    # prepare intermediates
    intermediate_path = "{}/PRESUMEout/intermediate".format(OUTDIR)
    os.makedirs(intermediate_path, exist_ok = True)
    decomp_nwk_path = "{}/decomp_tree".format(intermediate_path)

    # decompose tree

    # if a tree file was given
    if os.path.isfile (INFILE):
        tree = Phylo.read(INFILE, "newick")
        tree = rename_internals(tree)

        Phylo.write(tree, "{}/{}.nwk".format(args.output+"/PRESUMEout", "PRESUMEout"), "newick")

        tips_threshold = int(len(tree.get_terminals())**(1/2))
        os.makedirs(decomp_nwk_path, exist_ok = True)
        decompose(INFILE, tips_threshold, decomp_nwk_path)
    # if a directory containing decompsed tree files was given
    elif os.path.exists(INFILE):
        shutil.copytree(INFILE, decomp_nwk_path)
    else:
        print("Error: no such file or directory!", INFILE, file=sys.stderr)
        sys.exit()
        
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
    
    
    intermediate_fasta_path = "{}/fastas".format(intermediate_path)
    intermediate_indel_path = "{}/indels".format(intermediate_path)
    os.makedirs(intermediate_fasta_path, exist_ok = True)
    os.makedirs(intermediate_indel_path, exist_ok = True)

    # simulation of upper lineage
    print("NWK2FA: Simulating on upper trees...", end = '')
    upper_name2seq_without_indel, upper_name2seq, upper_name2alignedseq, upper_name2indellist, tree, lineage_tree = \
        nwk2fa_light(upper_tree, initseq, parsed_args, )   
    print("Finished!")

    #  mutation count per branch
    if (parsed_args.save_N_mutations):
        with open(args.output+"/PRESUMEout/mother_daughter_Nsubstitutions.up.csv",'a') as outfp:
            count_mutations_per_branch(lineage_tree, outfp)
        with open(args.output+"/PRESUMEout/position_Nsubstitutions.up.csv",'a') as outfp:
            count_mutations_per_position(lineage_tree, outfp)

    fasta_writer_single(upper_name2seq_without_indel, intermediate_fasta_path) # Seems tricky but "upper_name2seq_without_indel" shoule be appropriate here
    if (parsed_args.CRISPR): indel_writer_single(upper_name2indellist, intermediate_indel_path)

    fasta_filelist = os.listdir(intermediate_fasta_path)
    if(parsed_args.CRISPR): indel_filelist = os.listdir(intermediate_indel_path)
    nwk_filelist   = os.listdir(decomp_nwk_path)

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

    treefile_list, fastafile_list, indelfile_list = [], [], []
    for file_nwk in nwk_filelist:
        if file_nwk.split("_")[0] == "Down":
            filename_without_ext = file_nwk.split(".")[0]
            nodename = "_".join(filename_without_ext.split("_")[1:])
            file_fasta = "{}.fasta".format(nodename)

            path_of_fasta_in = "{}/{}".format(intermediate_fasta_path, file_fasta)
            path_of_newick_in = "{}/{}".format(decomp_nwk_path, file_nwk)

            fastafile_list.append(path_of_fasta_in)
            treefile_list.append(path_of_newick_in)
            
            if(parsed_args.CRISPR): 
                file_indel = "{}.indel".format(nodename)
                path_of_indel_in = "{}/{}".format(intermediate_indel_path, file_indel)
                indelfile_list.append(path_of_indel_in)
            else:
                indelfile_list.append(None)

    shell_path = "{}/shell".format(intermediate_path)
    os.makedirs(shell_path, exist_ok = True)
    stdeo_path = "{}/stdeo".format(intermediate_path)
    os.makedirs(stdeo_path, exist_ok = True)   
    print("NWK2FA: Simulating on bottom trees... ", end = '') 
    submit_command = shell_generator(
        shell_path, 
        treefile_list, 
        fastafile_list, 
        indelfile_list, 
        downstream_newick_path, 
        downstream_fasta_path,
        stdeo_path,
        args,
        parsed_args
        )
    subprocess.call(submit_command, shell=True)
    print("Finished!")

    command  = "cat {}/*/PRESUMEout/PRESUMEout.fa.gz > {}/PRESUMEout.fa.gz; ".format(downstream_fasta_path, OUTDIR+"/PRESUMEout")
    if(parsed_args.CRISPR):
        command += "cat {}/*/PRESUMEout/PRESUMEout.aligned.fa.gz > {}/PRESUMEout.aligned.fa.gz; ".format(downstream_fasta_path, OUTDIR+"/PRESUMEout")
        command += "cat {}/*/PRESUMEout/PRESUMEout.indel.gz > {}/PRESUMEout.indel.gz; ".format(downstream_fasta_path, OUTDIR+"/PRESUMEout")
    if(parsed_args.save_N_mutations):
        command += "cat {}/mother_daughter_Nsubstitutions.up.csv {}/*/PRESUMEout/mother_daughter_Nsubstitutions.csv > {}/mother_daughter_Nsubstitutions.csv; \
                    rm {}/mother_daughter_Nsubstitutions.up.csv; "\
                    .format(OUTDIR+"/PRESUMEout", downstream_fasta_path, OUTDIR+"/PRESUMEout",OUTDIR+"/PRESUMEout")
        command += "cat {}/position_Nsubstitutions.up.csv {}/*/PRESUMEout/position_Nsubstitutions.csv > {}/position_Nsubstitutions.csv; \
                    rm {}/position_Nsubstitutions.up.csv; "\
                    .format(OUTDIR+"/PRESUMEout", downstream_fasta_path, OUTDIR+"/PRESUMEout",OUTDIR+"/PRESUMEout")
    subprocess.call(command, shell=True)
    if (not args.debug) : shutil.rmtree(intermediate_path)
    print("Done!")
    return

def nwk2fa_single(args, parsed_args):
    # initial sequence specification
    initseq = parsed_args.initseq

    fasta_writer_multiple({"root":initseq}, args.output+"/PRESUMEout", "root")
    
    # read tree
    tree = Phylo.read(args.tree, "newick")
    tree = rename_internals(tree)

    # output
    upper_name2seq_without_indel, name2seq, name2alignedseq, name2indellist, newtree, lineage_tree = nwk2fa_light(tree, initseq, parsed_args=parsed_args)
    #  sequences
    fasta_writer_multiple(name2seq, args.output+"/PRESUMEout", "PRESUMEout")
    #  aligned sequences and indel list
    if parsed_args.CRISPR:
        fasta_writer_multiple(name2alignedseq, args.output+"/PRESUMEout", "PRESUMEout.aligned")
        indel_writer_multiple(name2indellist , args.output+"/PRESUMEout", "PRESUMEout")
    #  tree
    Phylo.write(newtree, "{}/{}.nwk".format(args.output+"/PRESUMEout", "PRESUMEout"), "newick")
    #  mutation count per branch
    if (parsed_args.save_N_mutations):
        with open(args.output+"/PRESUMEout/mother_daughter_Nsubstitutions.csv",'a') as outfp:
            count_mutations_per_branch(lineage_tree, outfp)
        with open(args.output+"/PRESUMEout/position_Nsubstitutions.csv",'a') as outfp:
            count_mutations_per_position(lineage_tree, outfp)