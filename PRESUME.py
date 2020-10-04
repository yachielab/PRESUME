#!/usr/bin/env python
# -*- coding: utf-8 -*-
__version__ = "1.0.0"

__author__ = "\
            Keito Watano <watano.k10.yachielab@gmail.com>,\
            Naoki Konno <naoki@bs.s.u-tokyo.ac.jp>, \
            Nozomu Yachie <nzmyachie@gmail.com>"

__date__ = "2020/6/14"

import os
os.environ["OMP_NUM_THREADS"] = "1"

import random
import numpy as np
import sys
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess
import shutil
import re
import csv
import copy
import gzip
from submodule import args_reader

LOGO = '''
######     ######     #######     #####     #     #    #     #    #######
#     #    #     #    #          #     #    #     #    ##   ##    #
#     #    #     #    #          #          #     #    # # # #    #
######     ######     #####       #####     #     #    #  #  #    #####
#          #   #      #                #    #     #    #     #    #
#          #    #     #          #     #    #     #    #     #    #
#          #     #    #######     #####      #####     #     #    #######

Version:     1.0.0
Last update: April 24, 2020
GitHub:      https://github.com/yachielab/PRESUME
'''

# set environment

# absolute path of python3 directory
PYTHON3 = (((
    subprocess.Popen('which python3', stdout=subprocess.PIPE, shell=True)
    .communicate()[0])
    .decode('utf-8'))
    .split('\n')
    )[0]

# absolute path of PRESUME directory
PRESUME = os.path.dirname(os.path.abspath(__file__))


#   for Exception
class PresumeException(Exception):
    """Base class for exceptions in this module."""
    pass


class ExtinctionError(PresumeException):
    def __init__(self):
        self.message = "PRESUME: An Error Occured !!!\n \
                        doubling time of initial SEQ is 0 !"

    def __str__(self):
        return self.message


class NoChildError(PresumeException):
    def __init__(self, seqid):
        self.seqid = seqid

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        No alive children of id:{}
        """.format(self.seqid)
        return message


class CreateNewickError(PresumeException):
    def __init__(self, e):
        self.occured_error = e

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        Something went wrong when we create newick...
        ERROR MESSAGE:{}
        """.format(self.occured_error)
        return message


class UpperLimitExceededError(PresumeException):
    def __init__(self, u):
        self.upper_limit = u

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        Upper limit exceeded !
        upper limit:{}
        """.format(self.upper_limit)
        return message


class TopologyError(PresumeException):
    def __init__(self, tree):
        self.tree = tree

    def __str__(self):
        num_of_internal_nodes = len(self.tree.get_terminals())
        num_of_terminal_nodes = len(self.tree.get_terminals())
        message = """
        PRESUME: An Error Occured !!!
        PRESUME generated invalid topology tree !
        # of internal nodes:{0}
        # of terminal nodes:{1}
        """.format(num_of_internal_nodes, num_of_terminal_nodes)
        return message


class OutputError(PresumeException):
    def __init__(self, fa_cnt, tip_cnt):
        self.fa_cnt = fa_cnt
        self.tip_cnt = tip_cnt

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        number of SEQs in fasta and in newick is different!
        # of sequences:{0}
        # of tip labels:{1}
        """.format(self.fa_cnt, self.tip_cnt)
        return message


class SimulationFailureError(PresumeException):
    def __init__(self, take_count):
        self.take_count = take_count

    def __str__(self):
        message = """
        PRESUME: An Error Occured !!!
        Recursive limit exceeded !
        Take:{}
        """.format(self.take_count)
        return message
    
class ExceptionOfCodeOcenan(PresumeException):
    def __init__(self):
        self.message="""
        PRESUME: Because Code Ocean does not support Univa Grid Engine, 
        distributed computing mode(the option of --qsub) is disabled.
        If you want to use it, please install PRESUME from http:/github.com/yachielab/PRESUME.
        """

    def __str__(self):
        return self.message


#   for Simulation
class SEQ():
    # idM & mseq means id & sequence of mother SEQ, respectively.
    def __init__(self, SEQid, idM, mseq, CVM, rM, dM, tM, indelsM, CV=False): # indels: list of [('in' or 'del', start_pos, length)]
        self.id = SEQid  # for example: 0,1,2,...
        self.idM = idM
        self.CV = max(np.random.normal(CVM, CVM*alpha), 0)  # should be > 0
        if CV:
            self.r = self.growing_rate_dist(rM, rM*self.CV)
        else:
            self.r = self.growing_rate_dist(rM, self.CV)
        
        self.is_alive = (self.r >= 0 and np.random.rand() > e)
        
        if(self.is_alive):
            if self.r == 0:
                self.d = float("inf")
            else:
                self.d = 1/self.r

            self.t      = tM + self.d  # time t of doubling of this SEQ
            self.seq    = self.daughterseq(str(mseq), dM)

            if (CRISPR): self.indels = indelsM + self.gen_indels() # CRISPR == True if an inprob file path is specified 
            else       : self.indels = None

            self.mutation_rate = compare_sequences(str(mseq), self.seq)

    def growing_rate_dist(self, mu, sigma):
        growing_rate = np.random.normal(mu, sigma)
        if growing_rate <= 0:
            growing_rate = -1
        return growing_rate

    # receive mother SEQ sequence, introduce mutations,
    # return daughter sequence.
    def daughterseq(self, seq, dM):
        dseq = ""
        for i in range(L):
            if(args.homoplasy is not None):
                dseq = dseq + self.homoplastic_mutation(seq, i, True)
            elif(args.constant is not None):
                dseq = dseq + self.time_independent_mutation(seq[i], mu[i])
            elif(args.gtrgamma is not None):
                dseq = dseq+self.time_dependent_mutation(seq[i], gamma[i], dM)
        return dseq

    # mutation of a site (NOT Jukes Cantor model.
    # directly define mutation matrix, not the mutation rate matrix
    # it's enough for calculate mutation of each duplication
    def homoplastic_mutation(self, c, mu):
        base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
        matrix = [
            [1-mu, mu/3, mu/3, mu/3],
            [mu/3, 1-mu, mu/3, mu/3],
            [mu/3, mu/3, 1-mu, mu/3],
            [mu/3, mu/3, mu/3, 1-mu]
            ]
        return random.choices(
            ['A', 'C', 'G', 'T'], k=1, weights=matrix[base[c]]
            )[0]
    
    def time_independent_mutation(self, c, mu):
        base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
        matrix = [
            [1-mu, mu/3, mu/3, mu/3],
            [mu/3, 1-mu, mu/3, mu/3],
            [mu/3, mu/3, 1-mu, mu/3],
            [mu/3, mu/3, mu/3, 1-mu]
            ]
        return random.choices(
            ['A', 'C', 'G', 'T'], k=1, weights=matrix[base[c]]
            )[0]
    
    # mutation of a site (GTR-Gamma model)
    def time_dependent_mutation(self, c, gamma, dM):
        base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        matrix = P(dM, gamma)

        # np.matrix[x] returns matrix, then the matrix is converted to array()
        return random.choices(
            ['A', 'C', 'G', 'T'], k=1, weights=np.array(matrix[base[c]])[0]
            )[0]

    def gen_indels(self): # pos2inprob, in_lengths, pos2delprob, del_lengths were defined from argument files

        def randomstr(alphabet, length):
            seq=""
            for _ in range(length):
                seq = seq + random.choice(alphabet)
            return seq

        seq_length = len(pos2inprob)
        
        shuffled_in_pos  = list(range(seq_length)); random.shuffle(shuffled_in_pos)
        shuffled_del_pos = list(range(seq_length)); random.shuffle(shuffled_del_pos)
        indel_order      = ['in', 'del'];                random.shuffle(indel_order)

        generated_indels = []

        for indel in indel_order:

            if (indel == 'in'):

                for pos in shuffled_in_pos:

                    if ( random.random() < pos2inprob[pos] ):

                        length = random.choices(in_lengths[0], k = 1, weights = in_lengths[1])[0]
                        
                        generated_indels.append( [ 'in',  pos, length, randomstr(['A','T','G','C'], length) ] )
            
            elif (indel == 'del' ):

                for pos in shuffled_del_pos:
                    
                    if ( random.random() < pos2delprob[pos] ):
                        
                        length = random.choices(del_lengths[0], k = 1, weights = del_lengths[1])[0]

                        generated_indels.append( [ 'del', pos, length ] )
            
        return generated_indels

def compare_sequences(motherseq, myseq):
    if len(motherseq) != len(myseq):
        print('Sequences should be same length')
        return None
    score = []  # mutation: ["1"], replication: ["0"]
    for position in range(len(motherseq)):
        its_a_match = motherseq[position] == myseq[position]
        score += [int(0)] if its_a_match else [int(1)]
    return score


def argument_saver(args):
    argname = dir(args)
    arg_pair = []
    csv_L1 = []
    csv_L2 = []
    while argname != []:
        arg = argname.pop()
        find_dunder = re.match('^_', arg)
        if not find_dunder:
            if args.__dict__[arg] is not None and not(arg in {"load", "save"}):
                csv_L1.append(arg)
                csv_L2.append(args.__dict__[arg])
        arg_pair = [csv_L1, csv_L2]

        with open("args.csv", "w") as fout:
            csvout = csv.writer(fout)
            csvout.writerows(arg_pair)


# when a child of SEQ of <SEQid> is dead,
# kill ancestors whose any descendant is dead.
def death(Lineage, SEQid):
    index = SEQid
    while(index != 0
            and Lineage[Lineage[index][1]] == ["dead"]
            and Lineage[Lineage[index][2]] == ["dead"]):
        idM = Lineage[index][0]
        Lineage[index] = ["dead"]
        index = idM


# if all SEQ died, create a file "all_SEQ_dead.out"
def all_dead(idANC):
    with open("all_SEQ_dead.out", 'w') as handle:
        handle.write(str(idANC)+"\n")


# count number & length of sequences in a specified fasta file
def count_sequence(in_fname):
    origin = gzip.open(in_fname, "rt")
    if(True):
    #with open(in_fname) as origin:
        input_itr = SeqIO.parse(origin, "fasta")
        # Build a list sequences:
        number_of_sequences = 0
        for sequ in input_itr:
            number_of_sequences += 1
    return number_of_sequences


def create_newick(Lineage, Timepoint, timelimit, upper_tree=False):
    def find_downstream_bifurcation(query, Lineage):
        if Lineage[query] != ['dead'] and len(Lineage[query]) == 1:
            return query
        both_children_are_alive = Lineage[Lineage[query][1]] != ['dead'] and Lineage[Lineage[query][2]] != ['dead']
        if both_children_are_alive:
            return query
        child = query
        while(not both_children_are_alive):
            if Lineage[Lineage[child][1]] != ['dead']:
                child = Lineage[child][1]
            else:
                child = Lineage[child][2]
            if len(Lineage[child])==1:
                break
            both_children_are_alive = Lineage[Lineage[child][1]] != ['dead'] and Lineage[Lineage[child][2]] != ['dead']
        return child

    if not upper_tree:
        for key in Timepoint.keys():
            if Timepoint[key] is None:
                pass
            elif Timepoint[key] > timelimit:
                Timepoint[key] = timelimit
    else:
        for key in Timepoint.keys():
            if Timepoint[key] is None:
                pass
            elif Timepoint[key] >= 2 * timelimit:
                Timepoint[key] = 3 * timelimit
    init_clade = Phylo.BaseTree.Clade(name="0")
    tree = Phylo.BaseTree.Tree(init_clade)
    list_of_dead = []
    stack = [init_clade]
    while(stack != []):
        clade = stack.pop()
        SEQid = int(clade.name)
        if(len(Lineage[SEQid]) == 3):  # if the SEQ has children
            # by this, every internal node will have 2 children
            while(True):
                both_children_are_alive = (
                    Lineage[Lineage[SEQid][1]] != ["dead"]
                    and Lineage[Lineage[SEQid][2]] != ["dead"])
                both_children_are_dead = (
                    Lineage[Lineage[SEQid][1]] == ["dead"]
                    and Lineage[Lineage[SEQid][2]] == ["dead"])
                if(both_children_are_alive):
                    child_L = find_downstream_bifurcation(Lineage[SEQid][1], Lineage)
                    child_R = find_downstream_bifurcation(Lineage[SEQid][2], Lineage)
                    children = [
                        Phylo.BaseTree.Clade(
                            name=str(child_L),
                            branch_length= Timepoint[child_L] - Timepoint[int(clade.name)]
                            ),
                        Phylo.BaseTree.Clade(
                            name=str(child_R),
                            branch_length= Timepoint[child_R] - Timepoint[int(clade.name)]
                            )
                        ]
                    clade.clades.extend(children)
                    stack.extend(children)
                    break
                elif(both_children_are_dead):
                    raise NoChildError(SEQid)
                    return
                # if either of the children is dead, renew SEQid
                # to be the id of alive child
                elif(Lineage[Lineage[SEQid][1]] != ["dead"]):
                    SEQid = Lineage[SEQid][1]
                elif(Lineage[Lineage[SEQid][2]] != ["dead"]):
                    SEQid = Lineage[SEQid][2]
                # if the clade of renewed SEQid is a terminal,
                # rename the clade
                if(len(Lineage[SEQid]) == 1):
                    clade.name = str(SEQid)
                    break

    # only in downstream tree of distributed computing,
    # the name of terminals is <upstream id>_<downstream id>
    if args.idANC != 0:
        for clade in tree.get_terminals():
            clade_name_prefix = str(hex(args.idANC)).split("x")[1]
            clade_name_suffix = str(hex(int(clade.name))).split("x")[1]
            new_clade_name = "{0}_{1}". \
                format(clade_name_prefix, clade_name_suffix)
            clade.name = new_clade_name

    # for error check
    if (len(tree.get_nonterminals()) + 1 != len(tree.get_terminals())):
        raise TopologyError(tree)
        return

    # file write in default
    if args.idANC == 0:
        Phylo.write(tree, "PRESUMEout.nwk", 'newick')

    else:
        # file write in case of downstream lineage
        Phylo.write(tree, "SEQ_"+str(args.idANC)+".nwk", 'newick')
    
    return len(tree.get_terminals()), tree, list_of_dead

    # except Exception as e:
    #     raise CreateNewickError(e)


def fasta_writer(name, seq, indels, file_name, overwrite_mode, Nchunks, filepath2writer=None, indelseq=True):

    if overwrite_mode:
        writer_mode = "a"
    else:
        writer_mode = "w"

    Lchunk = len(seq) // Nchunks

    is_dead = True
    true_indels = [] # List of indels 
    seq_list = []
    aligned_seq_list = []
    for chunkidx in range(Nchunks):

        chunk = seq[ Lchunk*chunkidx:Lchunk*(chunkidx+1) ]
        chunk_default = copy.deepcopy(chunk)

        # extract edits which occurred in the chunk
        ext_indels = []
        if (indels != None):
            for indel in indels:
                pos = indel[1]
                if ( pos//Lchunk == chunkidx ):
                    ext_indel    = copy.deepcopy(indel)
                    ext_indel[1] = pos % Lchunk
                    ext_indels.append((ext_indel,indel))

        #with open(chunk_file_name, writer_mode) as writer:
        if (True):

            if (indels != None):
                
                refpos2pos      = {}
                pos2refposindel = list(range(len(chunk)))
                for pos, refpos in enumerate(pos2refposindel):
                    refpos2pos[refpos] = pos

                for indel, global_indel in ext_indels:

                    if ( len (chunk) > 0 ):

                        if (indel[0] == 'del'):

                            if( indel[1] in refpos2pos.keys() ):
                                true_indels.append(global_indel)
                                mid    = refpos2pos[indel[1]] # pos is the midpoint of deletion
                                length = indel[2]
                                start  = max ( 0, mid - length//2 )
                                end    = min ( len(chunk) - 1, mid - length//2 + length - 1) 
                                chunk  = chunk[:start] + chunk[(end+1):]
                                pos2refposindel = pos2refposindel[:start] + pos2refposindel[(end+1):]
                                    
                                refpos2pos = {}
                                pos=0
                                for pos, refposindel in enumerate(pos2refposindel):
                                    if(type(refposindel)==int):
                                        refpos2pos[refposindel] = pos
                                    elif(refposindel=="in"):
                                        None

                        elif ( indel[0] == "in" ):

                            if( indel[1] in refpos2pos.keys() ):
                                true_indels.append(global_indel)

                                start  = refpos2pos[indel[1]]
                                length = indel[2]
                                chunk    = chunk[:start+1] + indel[3] + chunk[start+1:]
                                    
                                pos2refposindel = pos2refposindel[:start+1] + ["in"]*length + pos2refposindel[start+1:]

                                refpos2pos = {}
                                pos=0
                                for pos, refposindel in enumerate(pos2refposindel):
                                    if(type(refposindel)==int):
                                        refpos2pos[refposindel] = pos
                                    elif(refposindel=="in"):
                                        None
                        
                    else:
                        
                        dropout = True
                        break # True means the sequence died
                
            
            dropout = (len (chunk) <= 0)

            if(not indelseq):
                chunk = chunk_default

            if ( dropout ):
                aligned_chunk = "-"*Lchunk
            else:
                is_dead = False

                if (indelseq and indels != None):
                    aligned_chunk = "-"*Lchunk
                    for refpos in refpos2pos.keys():
                        aligned_chunk = aligned_chunk[:refpos] + chunk[refpos2pos[refpos]] + aligned_chunk[(refpos+1):]
            
            seq_list.append(chunk)
            if (indels != None and indelseq):
                aligned_seq_list.append(aligned_chunk)

    if (not is_dead):
    
        for chunkidx in range(Nchunks):

            # set output FASTA name
            if ( Nchunks > 1 ):
                chunk_file_name = file_name.split(".fa")[0]+"."+str(chunkidx)+".fa"
            else:
                chunk_file_name = file_name

            if filepath2writer == None:
                writer         = gzip.open(chunk_file_name + ".gz", writer_mode+"b")
            else:
                writer         = filepath2writer[chunk_file_name + ".gz"]
            
            # write unaligned seq.
            #SEQ_seq = SeqRecord(Seq(seq_list[chunkidx]))
            #SEQ_seq.id = str(name)
            #SEQ_seq.description = ""
            #SeqIO.write(SEQ_seq, writer, "fasta")
            writer.write((">"+str(name)+"\n").encode())
            writer.write((seq_list[chunkidx]+"\n").encode())
            if filepath2writer == None:
                writer.close()

            # write aligned seq.
            if (indels != None and indelseq):
                if ( Nchunks > 1 ):
                    aligned_chunk_file_name = file_name.split(".fa")[0]+"."+str(chunkidx)+".aligned.fa"
                else:
                    aligned_chunk_file_name = file_name.split(".fa")[0]+".aligned.fa"
                if filepath2writer == None:
                    aligned_writer = gzip.open(aligned_chunk_file_name + ".gz", writer_mode+"b")
                else:
                    aligned_writer = filepath2writer[aligned_chunk_file_name + ".gz"]
                #SEQ_seq = SeqRecord(Seq(aligned_seq_list[chunkidx]))
                #SEQ_seq.id = str(name)
                #SEQ_seq.description = ""
                #SeqIO.write(SEQ_seq, aligned_writer, "fasta")
                aligned_writer.write((">"+str(name)+"\n").encode())
                aligned_writer.write((aligned_seq_list[chunkidx]+"\n").encode())
                if filepath2writer == None:
                    aligned_writer.close()
            
            is_dead = False
    
    return true_indels, is_dead

def survey_all_dead_lineages(Lineage):
    try:
        command = "cat intermediate/DOWN/*/PRESUMEout/all_SEQ_dead.out \
            > intermediate/all_dead.out 2> /dev/null; \
            rm intermediate/DOWN/*/PRESUMEout/all_SEQ_dead.out 2> /dev/null"
        subprocess.call(command, shell=True)

    except Exception as e:
        print("no lineages extinct")
        return

    with open("intermediate/all_dead.out", 'r') as handle:
        lines = handle.readlines()
        for line in lines:
            dead_SEQ_id = int(line.split("\n")[0])
            idM = Lineage[dead_SEQ_id][0]
            Lineage[dead_SEQ_id] = ["dead"]
            death(Lineage, idM)


def CombineTrees():
    top_tree = Phylo.read('PRESUMEout.nwk', 'newick')
    terminals = top_tree.get_terminals()
    TOP_CELL_CNT = len(terminals)
    index = 0
    for tip in terminals:
        index += 1
        if(index % 100 == 0):
            message = "tree assembly:{}percent finished". \
                format(index * 100 / (TOP_CELL_CNT))
            sys.stdout.write(message)
            sys.stdout.flush()
        bottom_tree_filepath = \
            "intermediate/DOWN/esu_{0}/PRESUMEout/SEQ_{0}.nwk".format(tip.name)
        bottom_tree = Phylo.read(bottom_tree_filepath, 'newick')
        tip.clades.extend(bottom_tree.clade.clades)
        newick_from = \
            "intermediate/DOWN/esu_{0}/PRESUMEout/SEQ_{0}.nwk".\
            format(tip.name)
        newick_to = \
            "intermediate/DOWN/esu_{0}/PRESUMEout/SEQ_{0}_attached.nwk".\
            format(tip.name)
        shutil.move(newick_from, newick_to)
    Phylo.write(top_tree, 'PRESUMEout_combined.nwk', 'newick')
    return len(top_tree.get_terminals())


def progress_bar(pbar, current_value):
    if pbar.n <= current_value:
        pbar.n = current_value
    else:
        pbar.n = 0
    pbar.refresh()


def mut_rate_log_writer(L_dict, dead_lst):
    event = len(L_dict)
    seq_length = len(L_dict[next(iter(L_dict))])
    mutation_counter = np.zeros(seq_length)
    list_of_death=dead_lst
    for name in L_dict.keys():
        if name not in list_of_death:
            score = np.array(L_dict[name])
            mutation_counter += score

    print("[DEBUG]dead sequences:", len(dead_lst))
    print("[DEBUG]event:", event)
    print("[DEBUG]seq_length:", seq_length)

    with open("mut_rate_log.csv", "w") as f:
        csvout = csv.writer(f)
        csvout.writerow(["position", "mutation_rate"])
        for idx, item in enumerate(mutation_counter):
            csvout.writerow([idx+1, item])


def unbalance(clade):
    if len(clade.clades) == 2:
        child1 = len(
            list(Phylo.BaseTree.Tree(clade.clades[0]).get_terminals())
            )
        child2 = len(
            list(Phylo.BaseTree.Tree(clade.clades[1]).get_terminals())
            )
        return min(child1, child2)/max(child1, child2)
    else:
        return 0


def Um_analyzer(tree):
    quene = [tree.clade]
    list_of_u = []  # Um
    while(quene != []):
        m_clade = quene.pop()
        if len(Phylo.BaseTree.Tree(m_clade).get_terminals()) > 10:
            Um = unbalance(m_clade)
            list_of_u.append(Um)
            quene.extend(m_clade.clades)
    return list_of_u

def make_list(txtfile, datatype, column):
    List=[]
    with open(txtfile, 'r') as handle:
        for line in handle:
            if(datatype=='float'): List.append(float(line.split("\n")[0].split("\t")[column]))
            elif(datatype=='int'): List.append(int  (line.split("\n")[0].split("\t")[column]))
    return List

# main
def main(timelimit):
    '''
    main function
            commandline argument: "python3 PRESUME.py <timelimit> <options...>"
    '''
    # note: tqdm should be ver. 4.19 (conda install tqdm=4.19)
    if args.bar:
        from tqdm import tqdm
    C = args.n  # C : number of SEQs you want to get (default == None)
    if(C is not None):  # if C is specified
        if(args.qsub):  # for distributed computing
            C = int(C**(1/2))  # update C

    if args.bar:
        pbar = tqdm(range(C))

    # for safety: program definitely stops when the number of SEQs
    # exceeds UPPER_LIMIT
    UPPER_LIMIT = args.u
    delta_timelimit = timelimit
    inittimelimit = timelimit

    # next id
    i = 0

    # current existing SEQ
    SEQqueue = []

    # Lineage[j]=[k,l,m] means an SEQ whose id is j is a daughter of Lineage[k]
    # and is a mother of Lineage[l], Lineage[m]
    Lineage = [[]] * UPPER_LIMIT
    Timepoint = {}

    # DEBUG: for gathering mutation rate of each site
    if args.debug:
        mut_rate_log = {}

    # First of all, there exits only 1 SEQ.
    if args.CV:
        SEQqueue.append(
            SEQ(i, -1, initseq, sigma_origin, growing_rate_origin,
                dorigin, args.tMorigin, initindels, True))
    else:
        SEQqueue.append(
            SEQ(i, -1, initseq, sigma_origin, growing_rate_origin,
                dorigin, args.tMorigin, initindels, True))
    SEQqueue[0].seq = initseq
    SEQqueue[0].indels = initindels
    if args.idANC == 0:
        Timepoint[0] = SEQqueue[0].t if SEQqueue[0].is_alive else None
    else:
        Timepoint[0] = args.tMorigin + dorigin - (timelimit / 2)
    i += 1
    c  = 1  # current number of SEQs

    # SEQs propagation
    while(True):
        if (c == 0):
            all_dead(args.idANC)
            try:
                os.remove("ancestoral_sequences.fasta")
            except FileNotFoundError:
                pass
            if args.save:
                argument_saver(args)
            return 1

        # SEQs divide until t (time) of all of them exceed timelimit
        k = 0  # iterator of SEQqueue

        # while there remains some SEQs whose time is less than timelimit
        while(k < c):

            if (not SEQqueue[k].is_alive):
                # Note: Variable name "esu" means Evolutionally Sequence Unit
                # as container of SEQ object
                esu = SEQqueue.pop(k)
                Lineage[esu.id] = ["dead"]
                if(esu.id != 0):
                    death(Lineage, esu.idM)
                c -= 1
                if args.bar:
                    progress_bar(pbar, c)

            elif(SEQqueue[k].t < timelimit):
                esu = SEQqueue.pop(k)
                # save ancestoral sequences
                if args.viewANC:
                    true_indels, is_dead = fasta_writer(esu.id, esu.seq, esu.indels, "ancestral_sequences.fasta", True, Nchunks = args.chunks)
                # duplication
                if args.CV:
                    daughter = [SEQ(i, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t, esu.indels, True),
                                SEQ(i+1, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t, esu.indels, True)]
                else:
                    daughter = [SEQ(i, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t, esu.indels),
                                SEQ(i+1, esu.id, esu.seq, esu.CV,
                                    esu.d, esu.r, esu.t, esu.indels)]

                SEQqueue.extend(daughter)

                if args.debug:
                    for sister in daughter:
                        if sister.is_alive:
                            mut_rate_log[sister.id]=sister.mutation_rate

                # [<mother>, <daughter1>, <daughter2> ]
                Lineage[esu.id] = [esu.idM, i, i+1]
                if daughter[0].is_alive:
                    Timepoint[daughter[0].id] = daughter[0].t
                if daughter[1].is_alive:
                    Timepoint[daughter[1].id] = daughter[1].t

                try:
                    Lineage[i] = [esu.id]
                    Lineage[i+1] = [esu.id]  # Terminal ESUs

                except IndexError:
                    print("UPPER LIMIT EXCEEDED!!!")
                    return 1

                i += 2
                c += 1
                if args.bar:
                    progress_bar(pbar, c)

            else:
                k += 1

            if(c > UPPER_LIMIT):  # for safety
                raise UpperLimitExceededError(UPPER_LIMIT)
                return 1

        if (c == 0):
            all_dead(args.idANC)
            try:
                os.remove("ancestoral_sequences.fasta")
            except FileNotFoundError:
                pass
            if args.save:
                argument_saver(args)
            return 1

        # if the required number of SEQs is not specified (default)
        if(C is None):
            break
        # if the required number of SEQs specified
        else:
            if(c < C):
                delta_timelimit = inittimelimit/(timelimit/inittimelimit)
                timelimit += delta_timelimit
            else:
                print("Number of generated sequences reached "+str(C))
                break


    # in case of "sequential computing"
    # or "downstream SEQ simulation of distributed computing"
    if(not args.qsub):
        # output initial sequence
        fasta_writer("root", initseq, None, "root.fa", False, Nchunks = args.chunks)
        # create fasta
        print("Generating a FASTA file...")

        # open .gz file to write
        filepath2writer={}
        if (args.chunks == 1):
            filepath2writer["PRESUMEout.fa.gz"]                = gzip.open("PRESUMEout.fa.gz"           , 'wb')
            filepath2writer["PRESUMEout.aligned.fa.gz"]        = gzip.open("PRESUMEout.aligned.fa.gz"   , 'wb')
        else:
            for i in range(args.chunks):
                filepath2writer["PRESUMEout."+str(i)+".fa.gz"]         = gzip.open("PRESUMEout."+str(i)+".fa.gz", 'wb')
                filepath2writer["PRESUMEout."+str(i)+".aligned.fa.gz"] = gzip.open("PRESUMEout."+str(i)+".aligned.fa.gz", 'wb')

        esuname2zerolength={}
        for esu in SEQqueue:
            if(args.idANC is None):
                esu_name = str(esu.id)
            else:
                esu_name_prefix = str(hex(args.idANC)).split("x")[1]
                esu_name_suffix = str(hex(esu.id)).split("x")[1]
                new_esu_name = "{}_{}".\
                    format(esu_name_prefix, esu_name_suffix)
                esu_name = new_esu_name
            
            true_indels, esu_zero_length = fasta_writer(esu_name, esu.seq, esu.indels, "PRESUMEout.fa", True, filepath2writer=filepath2writer, Nchunks=args.chunks) # fasta_writer() returns True if the seq. length <= 0
            esu.indels = true_indels
            esuname2zerolength[esu_name] = esu_zero_length

            if (esu_zero_length):
                Lineage[esu.id] = ["dead"]
                if(esu.id != 0):
                    death(Lineage, esu.idM)
                c -= 1

        if (c == 0):
            all_dead(args.idANC)
            fa_count = 0
            tip_count = 0

        else: 

            # record indels
            if(CRISPR):
                handle = gzip.open("indel.txt.gz", "wb")
                if (True):
                #with open("indel.txt", 'w') as handle:
                    for esu_idx, esu in enumerate(SEQqueue):

                        if(args.idANC is None):
                            esu_name = str(esu.id)
                        else:
                            esu_name_prefix = str(hex(args.idANC)).split("x")[1]
                            esu_name_suffix = str(hex(esu.id)).split("x")[1]
                            new_esu_name = "{}_{}".\
                                format(esu_name_prefix, esu_name_suffix)
                            esu_name = new_esu_name

                        if (not esuname2zerolength[esu_name]):
                            handle.write((esu_name+"\t").encode())
                            
                            for indel_idx, indel in enumerate(esu.indels):
                                if (indel[0] == "del") : 
                                    mid    = indel[1] # pos is the midpoint of deletion
                                    length = indel[2]
                                    handle.write(("D_mid"+str(mid)+"_len"+str(length)).encode())
                                elif (indel[0] == "in"): 
                                    pos    = indel[1] # pos is the midpoint of deletion
                                    seq    = indel[3]
                                    handle.write(("I_"+str(pos+0.5)+"_"+seq).encode())
                                if (indel_idx < len(esu.indels)-1):
                                    handle.write(";".encode())
                            handle.write("\n".encode())
                handle.close()

            # create newick
            del(SEQqueue)
            print("Generating a Newick file......")
            tip_count, returned_tree, list_of_dead = create_newick(Lineage, Timepoint, timelimit)

            if args.debug:
                mut_rate_log_writer(mut_rate_log, list_of_dead)

            # close FASTA
            for writer in list(filepath2writer.values()):
                writer.close()

    # in case of distributed computing
    else :
        # Root sequence
        fasta_writer("root", initseq, None, "root.fa", False, Nchunks = args.chunks)
        # preparation for qsub
        os.mkdir("intermediate")
        os.mkdir("intermediate/DOWN")
        os.mkdir("intermediate/fasta")
        os.mkdir("intermediate/shell")

        if(CRISPR): os.mkdir("intermediate/indel")

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

        itr = 0
        for esu in SEQqueue:

            fasta_file_path = \
                "intermediate/fasta/{}.fa".\
                format(str(esu.id))
            true_indels, esu_zero_length = fasta_writer(esu.id, esu.seq, esu.indels, fasta_file_path, True,  Nchunks=1, indelseq=False)
            esu.indels = true_indels

            if (esu_zero_length): # if the sequence length <= 0
                Lineage[esu.id] = ["dead"]
                if(esu.id != 0):
                    death(Lineage, esu.idM)
                c -= 1
            
            else:
                itr += 1
                if(CRISPR):
                    indel_file_path =\
                        "intermediate/indel/{}.txt".\
                        format(str(esu.id))
                    with open(indel_file_path, 'w') as handle:
                        for indel in esu.indels:
                            if (indel[0] == "del") : handle.write(str(indel[0])+"\t"+str(indel[1])+"\t"+str(indel[2])+"\n")
                            elif (indel[0] == "in"): handle.write(str(indel[0])+"\t"+str(indel[1])+"\t"+str(indel[2])+"\t"+str(indel[3])+"\n")

                with open("intermediate/shell/esu_"+str(itr)+".sh", 'w') as qf:
                    qf.write("#!/bin/bash\n")
                    qf.write("#$ -S /bin/bash\n")
                    qf.write("#$ -cwd\n")
                    qf.write("PATH={}\n".format(PATH))
                    qf.write("LD_LIBRARY_PATH={}\n".format(LD_LIBRARY_PATH))
                    qf.write("mkdir intermediate/DOWN/esu_"+str(esu.id)+"\n")
                    qf.write("cd intermediate/DOWN/esu_"+str(esu.id)+"\n")
                    qf.write("pwd\n")

                    # divide until time point of (2 * timelimit)
                    python_command = PYTHON3 + " " + PRESUME + "/PRESUME.py "\
                        "--monitor " + str(2*timelimit)\
                        + " -L " + str(L)\
                        + " -f " + "../../fasta/"+str(esu.id)+".fa.gz"\
                        + " -d " + str(esu.d)\
                        + " -s " + str(sigma_origin)\
                        + " -T " + str(T)\
                        + " -e " + str(e)\
                        + " -u " + str(UPPER_LIMIT)\
                        + " --idANC "    + str(esu.id)\
                        + " --tMorigin " + str(esu.t - esu.d)\
                        + " --seed "     + str(np.random.randint(0, 10000))\
                        + " --chunks "   + str(args.chunks)
                    if args.CV:
                        python_command += " --CV"
                    
                    if CRISPR:
                        python_command += \
                            " --inprob "     + args.inlsprob   +\
                            " --inlength "   + args.inlength +\
                            " --delprob "    + args.delprob  +\
                            " --dellength "  + args.dellength+\
                            " --indels "     + "../../indel/"+str(esu.id)+".txt"

                    qf.write(python_command)
                    if (args.gtrgamma is not None):
                        qf.write(" --gtrgamma "+str(args.gtrgamma)+"\n")
                    if (args.constant is not None):
                        qf.write(" --constant "+str(args.constant)+"\n")

        del(SEQqueue)
        # submit job script to grid engine
        print("\ncreating bottom trees by qsub ...")
        submit_command = "qsub " + args.dop + " -sync y -t 1-{0} \
            {1}/exe_PRESUME.sh &> intermediate/qsub.out".\
            format(str(itr), PRESUME)

        subprocess.call(submit_command, shell=True)

        # finalize

        # remove extinct downstream lineages
        survey_all_dead_lineages(Lineage)

        tip_count, returned_tree, list_of_dead = create_newick(Lineage, Timepoint, timelimit, True)
        if args.debug:
            mut_rate_log_writer(mut_rate_log, list_of_dead)

        command = "cat PRESUME.e*.* > intermediate/err 2> /dev/null; \
                   cat PRESUME.o*.* > intermediate/out 2> /dev/null; \
                   rm PRESUME.*;"
        #subprocess.call(command, shell=True)
        if args.f is None:
            if(CRISPR):
                command = "cat intermediate/DOWN/*/PRESUMEout/indel.txt.gz > indel_combined.txt.gz 2> /dev/null;"
                
            if (args.chunks == 1):
                command += "cat intermediate/DOWN/*/PRESUMEout/PRESUMEout.fa.gz > PRESUMEout.fa.gz;"
                if(CRISPR):
                    command += "cat intermediate/DOWN/*/PRESUMEout/PRESUMEout.aligned.fa.gz > PRESUMEout.aligned.fa.gz 2> /dev/null;"
            else:
                for i in range(args.chunks):
                    command += "cat intermediate/DOWN/*/PRESUMEout/PRESUMEout."+str(i)+".fa.gz > PRESUMEout."+str(i)+".fa.gz;"
                    if(CRISPR):
                        command += "cat intermediate/DOWN/*/PRESUMEout/PRESUMEout."+str(i)+".aligned.fa.gz > PRESUMEout."+str(i)+".aligned.fa.gz 2> /dev/null;"
            subprocess.call(command, shell=True)  # combine fasta

            
        tip_count = CombineTrees()  # Combine trees
        # shutil.move("PRESUMEout.nwk", "intermediate")
        # os.rename("PRESUMEout_combined.nwk", "PRESUMEout.nwk")
        if (not args.debug):
            shutil.rmtree("intermediate")

    # error check
    if(args.chunks == 1):
        fa_count = count_sequence("PRESUMEout.fa.gz")
        if (fa_count!=tip_count):
            raise OutputError(fa_count, tip_count)
            return 0

    # remove PRESUME.aligned.fa.gz when indel mode is not active
    if (not CRISPR):
        if (os.path.isfile("PRESUMEout.aligned.fa.gz")):
            os.remove("PRESUMEout.aligned.fa.gz")

    # finish
    if args.save:
        argument_saver(args)

    print("\n=====================================================")
    print("Simulation end time point:         "+str(2*timelimit))
    print("Number of generated sequences:     "+str(tip_count))
    print("Seed for random number generation: "+str(seed))
    print("=====================================================\n")
    return 0


def recursive_main(timelimit, limit, main_func, repeated=1):
    result = main_func(timelimit)
    if limit < repeated:
        raise SimulationFailureError(repeated)
    if result == 0:
        return repeated
    else:
        return recursive_main(timelimit, limit, main_func, repeated=repeated+1)


#   interface
if __name__ == "__main__":
    # interface
    parser = argparse.ArgumentParser(description='PRESUME.py', add_help=True)
    parser.add_argument(
        "--param",
        help="load argument file(csv file)",
        type=str
        )
    parser.add_argument(
        "-V",
        "--version",
        action="store_true",
        default=False
        )

    parser.add_argument(
        "--monitor",
        help="time limit (default=1)",
        type=float,
        default=1
        )

    parser.add_argument(
        "-n",
        help="required number of sequences: if you specified this parameter,\
            timelimit will be postponed until the number of the sequence reach\
            the specified number (default=None)",
        type=int
        )

    parser.add_argument(
        "-L",
        help="length of sequence (default=1000)",
        type=int,
        default=1000
        )

    parser.add_argument(
        "-f",
        help="fasta file name　of the common ancestor sequence.\
            (default: poly-C)",
        type=str
        )

    parser.add_argument(
        "-d",
        help="doubling time of origin sequence (default=1)",
        type=float,
        default=1
        )

    parser.add_argument(
        "-s",
        help="sigma of doubling time of origin sequence (default=0)",
        type=float,
        default=0
        )

    parser.add_argument(
        "-a", help="CV of CV of doubling time of origin sequence (default=0)",
        type=float,
        default=0
        )

    parser.add_argument(
        "-T",
        help="Threashold of doubling time to be deleted (default=1000)",
        type=float,
        default=1000
        )

    parser.add_argument(
        "-e",
        help="random deletion probability (default=0)",
        type=float,
        default=0
        )

    parser.add_argument(
        "--gtrgamma",
        help="parameters for substitution rate matrix\n \
            GTR{A-C/A-G/A-T/C-G/C-T/G-T} \
            +FU{piA/piC/piG/piT} \
            +G4{shape of gamma distribution}\n \
            Or, you can use default parameters by \"--gtrgamma default\"\n \
            default: \
            GTR{0.3333/0.3333/0.3333/0.3333/0.3333/0.3333} \
            +FU{0.25/0.25/0.25/0.25} \
            +G4{10000}",
        type=str
        )

    parser.add_argument(
        "-u",
        help="upper limit of number of sequences (default=2^20)",
        type=int,
        default=np.power(2, 20)
        )

    parser.add_argument(
        "-m",
        help="mean of relative  \
        subtitution rate according to gamma distribution \
            (default=1)",
        type=float,
        default=1
        )

    parser.add_argument(
        "--constant",
        help="fixed mutation rate of each site (default=None)",
        type=float)

    parser.add_argument(
        "--output",
        help="output folder (default:current directory)",
        type=str,
        default=os.getcwd()
        )

    # for distributed computing
    parser.add_argument(
        "--qsub",
        help="activate preparation for distributed processes\
            (defalt=inactivated)",
        action="store_true",
        default=False
        )

    # for distributed computing
    parser.add_argument(
        "--idANC",
        help="corresponging ancestral sequence (in upstream tree), \
            in case of distributed computing (default=None)",
        type=int,
        default=0
        )

    # for distributed computing
    parser.add_argument(
        "--tMorigin",
        help="birth time of origin sequence",
        type=float,
        default=0
        )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="inactivate deletion of intermediate files"
        )

    parser.add_argument(
        "--bar",
        help="deactivate unstable functions",
        action="store_true",
        default=False
        )

    # for debug
    parser.add_argument(
        "--viewANC",
        help="generate fasta of ancestoral sequences",
        action="store_true",
        default=False
        )

    parser.add_argument(
        "--save",
        help="generate args.csv",
        action="store_true",
        default=False
        )

    parser.add_argument(
        "--CV",
        help="sigma use as CV(Coefficient Variance) of Normal Distribution",
        action='store_true',
        default=False
        )

    parser.add_argument(
        "-r",
        help="limit of retrying simulation (default=100000)",
        type=int,
        default=100000
        )

    parser.add_argument(
        "--seed",
        help="random seed used to initialize \
            the pseudo-random number generator",
        type=int,
        default=None
        )
    '''
    parser.add_argument(
        "-h", "--help",
        help="print help document",
        action='store_true',
        default=False
        )
    '''

    parser.add_argument(
        "--polyC",
        help="use polyC sequence as root ",
        action='store_true',
        default=False
        )

    parser.add_argument(
        "--tree",
        help="file name of a guide tree in Newick format(for debug, unsupported.)",
        type=str,
        default=None,
    )


    parser.add_argument(
        "--inprob",
        help="file name of insertion probability for each position",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--inlength",
        help="file name of insertion length distribution",
        type=str,
        default=None,
    )
    
    parser.add_argument(
        "--delprob",
        help="file name of deletion probability for each position",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--dellength",
        help="file name of insertion length distribution",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--indels",
        help="file name of indels accumulated before simualtion (for distributed computing mode)",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--chunks",
        help="number of chunks for each sequence unit",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--dop",
        help="Option of qsub for downstream simulation (for distributed computing mode)",
        type=str,
        default="",
    )

    parser.add_argument(
        "--homoplasy",
        help="Option of enable homoplasy model",
        action="store_true",
        default=None
    )

    args = parser.parse_args()
    
    # "Code Ocean" specific process
    # if args.qsub:
    #     raise ExceptionOfCodeOcenan
    #     quit()
    
    #   to show help
    '''
    if args.help:
        import PRESUME_help as ph
        print(ph.help_description())
        exit()
    '''

    if args.version:
        print(LOGO)
        exit()

    OUTDIR = args.output
    
    if args.tree:
        if (os.path.exists(os.getcwd() + "/" + args.tree)):
            args.tree      = os.getcwd() + "/" + args.tree
        processed_args = args_reader.PARSED_ARGS(args)

        import nwk2fa as n2f
        print("tree mode!")
        os.chdir(OUTDIR)
        os.makedirs("PRESUMEout", exist_ok=True)
        os.chdir("PRESUMEout")

        if (not args.qsub):
            n2f.nwk2fa_single(args, processed_args)
        else:
            n2f.nwk2fa_qsub(args, processed_args)
        exit()
    
    # read argument from input CSV
    if args.param:
        with open(args.param, "rt") as fin:
            cin = csv.reader(fin)
            arglist = [row for row in cin]
            valid_format = len(arglist[0]) == len(arglist[1])
            if not valid_format:
                print("ERROR:invalid format!")
                quit()
            for position in range(len(arglist[0])):
                if len(arglist[1][position].split('.')) != 1:
                    loaded_arg = float(arglist[1][position])
                    args.__dict__[arglist[0][position]] = loaded_arg
                elif arglist[1][position] in ["True", "False"]:
                    if arglist[1][position] == "True":
                        args.__dict__[arglist[0][position]] = True
                    else:
                        args.__dict__[arglist[0][position]] = False
                elif len(arglist[1][position].split('/')) != 1:
                    loaded_arg = str(arglist[1][position])
                    args.__dict__[arglist[0][position]] = loaded_arg
                else:
                    loaded_arg = int(arglist[1][position])
                    args.__dict__[arglist[0][position]] = loaded_arg
            print("CSV loaded!")

    # --gamma xor --constant should be specified
    if(args.gtrgamma is None and args.constant is None):
        print("please specify --gtrgamma or --constant")
        sys.exit(1)
    if(args.gtrgamma is not None and args.constant is not None):
        print("please don't specify both --gtrgamma and --constant")
        sys.exit(1)

    # initialize the pseudo-random number generator
    if args.seed is not None:
        seed = args.seed
    elif args.seed == "rand":
        seed = np.random.randint(0, args.r)
    elif args.seed is None:
        seed = 0
    else:
        seed = int(args.seed)
    np.random.seed(int(seed))
    random.seed(int(seed))


    # parameters###### (corresponding to the Figure.2a)
    L = args.L
    dorigin = args.d
    if dorigin == 0:
        print("fatal error: doubling time of initial SEQ is 0!")
        sys.exit(1)

    growing_rate_origin = 1 / args.d
    sigma_origin = args.s
    alpha = args.a
    T = args.T
    e = args.e
    m = args.m

    # In case expected number of mutation is independent
    # on the doubling time of the SEQ
    if args.homoplasy is not None:
        None
    elif(args.constant is not None):
        mu = [args.constant] * L

    # substitution rate matrix (GTR-Gamma model)
    # -> define substitution matrix function
    if(args.gtrgamma is not None):
        if(args.gtrgamma == "default"):
            model = "GTR{0.03333/0.03333/0.03333/0.03333/0.03333/0.03333} \
            +FU{0.25/0.25/0.25/0.25} \
            +G4{10000}"
        else:
            model = args.gtrgamma
        models_str = str(model).split("+")
        abcdef = (models_str[0].split("{"))[1].split("}")[0].split("/")
        piACGT = (models_str[1].split("{"))[1].split("}")[0].split("/")
        gamma_str = (models_str[2].split("{"))[1].split("}")[0].split("/")
        a1 = float(abcdef[0])  # A <-> C
        a2 = float(abcdef[1])  # A <-> G
        a3 = float(abcdef[2])  # A <-> T
        a4 = float(abcdef[3])  # C <-> G
        a5 = float(abcdef[4])  # C <-> T
        a6 = float(abcdef[5])  # G <-> T
        piA = float(piACGT[0])
        piC = float(piACGT[1])
        piG = float(piACGT[2])
        piT = float(piACGT[3])
        if (abs(piA+piC+piG+piT-1) > 0.001):
            print("error piA+piC+piG+piT not equal to 1!")
            sys.exit(1)
        # substitution rate matrix
        R = np.matrix([
            [-(a1*piC+a2*piG+a3*piT), a1*piC, a2*piG, a3*piT],
            [a1*piA, -(a1*piA+a4*piG+a5*piT), a4*piG, a5*piT],
            [a2*piA, a4*piC, -(a2*piA+a4*piC+a6*piT), a6*piT],
            [a3*piA, a5*piC, a6*piG, -(a3*piA+a5*piC+a6*piG)]])
        print("Substitution rate matrix:")
        print(R)
        # Al: eigen values (np.array),
        # U: eigen vectors matrix :
        # R = U * diag(Al) * U^(-1)
        Al, U = np.linalg.eig(R)

        # return transition matrix
        def P(t, gamma):
            exp_rambda = np.diag(
                    np.array([
                        np.exp(Al[0] * t * gamma),
                        np.exp(Al[1] * t * gamma),
                        np.exp(Al[2] * t * gamma),
                        np.exp(Al[3] * t * gamma)]
                    )
                )
            return np.dot(np.dot(U, exp_rambda), np.linalg.inv(U))
        # model of site heterogeneity:
        # calculate relative substitution rate gamma for each site
        shape = float(gamma_str[0])  # shape of gamma distribution
        gamma = np.random.gamma(shape, m / shape, L)  # mean is args.m

    # set indel parameters from raw lab experiments
    CRISPR      = False
    if (args.inprob != None): 
        args.inprob    = os.path.abspath(args.inprob)
        args.inlength  = os.path.abspath(args.inlength)
        args.delprob   = os.path.abspath(args.delprob)
        args.dellength = os.path.abspath(args.dellength)
        pos2inprob  =   make_list(args.inprob  , 'float', column = 0)
        in_lengths  = [ make_list(args.inlength, 'int'  , column = 0)   ,
                        make_list(args.inlength, 'float', column = 1)   ]
        pos2delprob =   make_list(args.delprob , 'float', column = 0)
        del_lengths = [ make_list(args.dellength,'int'  , column = 0)   ,
                        make_list(args.dellength,'float', column = 1)   ]
        CRISPR      = True

    # initial sequence specification
    if (args.f is not None):

        handle = gzip.open(args.f, 'rt')
        sequences = SeqIO.parse(handle, 'fasta')
        initseq = str(list(sequences)[0].seq)
        L = len(initseq)
        handle.close()

    elif (args.polyC):
        initseq = 'C' * L  # initial sequence
    else:
        initseq = ''.join([
            np.random.choice(['A', 'G', 'C', 'T']) for i in range(L)
            ])
    
    # initial indel specification (for distributed computing mode)
    if(CRISPR):
        if (args.indels is not None):
            initindels = []
            with open(args.indels, 'r') as handle:
                for line in handle:
                    chunks = line.split()
                    if (chunks[0] == "del") : initindels.append([chunks[0], int(chunks[1]), int(chunks[2])])
                    elif (chunks[0] == "in"): initindels.append([chunks[0], int(chunks[1]), int(chunks[2]), chunks[3]])
        else:
            initindels=[]
    else:
        initindels=None

    # setup directory
    if not os.path.exists(OUTDIR):
        os.makedirs(OUTDIR)
    os.chdir(OUTDIR)
    os.makedirs("PRESUMEout", exist_ok=True)
    os.chdir("PRESUMEout")

    #   excecute main
    print(LOGO)
    if args.r == 1:
        return_main = main(args.monitor)
    else:
        counter = recursive_main(args.monitor, args.r, main)
#       print("Number of retrials:", counter)