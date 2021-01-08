import csv
import sys
import numpy as np
import random
import os
import gzip
from Bio import SeqIO

class PARSED_ARGS():
    ###########################################################################
    # idM & mseq means id & sequence of mother SEQ, respectively.
    ###########################################################################
    def __init__(self, args): # indels: list of [('in' or 'del', start_pos, length)]

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
        if(args.gtrgamma is None and args.constant is None and args.editprofile is None):
            print("please specify --gtrgamma, --constant, or --editprofile")
            sys.exit(1)
        if(args.gtrgamma is not None and args.constant is not None):
            print("please don't specify both --gtrgamma and --constant")
            sys.exit(1)

        if args.seed is not None:
            self.seed = args.seed
        elif args.seed == "rand":
            self.seed = np.random.randint(0, args.r)
        elif args.seed is None:
            self.seed = 0
        else:
            self.seed = int(args.seed)

        # parameters###### (corresponding to the Figure.2a)
        self.L = args.L
        self.dorigin = args.d
        if self.dorigin == 0:
            print("fatal error: doubling time of initial SEQ is 0!")
            sys.exit(1)

        self.growing_rate_origin = 1 / args.d
        self.sigma_origin = args.s
        self.T = args.T
        self.e = args.e
        self.m = args.m
        self.n = args.n
        self.constant = args.constant
        self.gtrgamma = args.gtrgamma
        self.save_N_mutations = args.debug
        self.editprofile      = args.editprofile
        self.chunks = 1

        # In case expected number of mutation is independent
        # on the doubling time of the SEQ
        if args.editprofile is not None:
            self.sub_prob_mtx_list = self.sub_mat_parser(args.editprofile)
            self.L = len(self.sub_prob_mtx_list)
            if (os.path.exists(os.getcwd() + "/" + args.editprofile)):
                args.editprofile       = os.getcwd() + "/" + args.editprofile
        elif(args.constant is not None):
            self.mu = [args.constant] * self.L

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
            self.R = np.matrix([
                [-(a1*piC+a2*piG+a3*piT), a1*piC, a2*piG, a3*piT],
                [a1*piA, -(a1*piA+a4*piG+a5*piT), a4*piG, a5*piT],
                [a2*piA, a4*piC, -(a2*piA+a4*piC+a6*piT), a6*piT],
                [a3*piA, a5*piC, a6*piG, -(a3*piA+a5*piC+a6*piG)]])
            print("Substitution rate matrix:")
            print(self.R)
            # Al: eigen values (np.array),
            # U: eigen vectors matrix :
            # R = U * diag(Al) * U^(-1)
            self.Al, self.U = np.linalg.eig(self.R)

            # model of site heterogeneity:
            # calculate relative substitution rate gamma for each site
            shape = float(gamma_str[0])  # shape of gamma distribution
            self.gamma = np.random.gamma(shape, args.m / shape, self.L)  # mean is args.m

        # set indel parameters from raw lab experiments
        self.CRISPR      = False
        if (args.inprob != None): 
            def make_list(txtfile, datatype, column):
                List=[]
                with open(txtfile, 'r') as handle:
                    for line in handle:
                        if(datatype=='float'): List.append(float(line.split("\n")[0].split("\t")[column]))
                        elif(datatype=='int'): List.append(int  (line.split("\n")[0].split("\t")[column]))
                return List
            args.inprob    = os.path.abspath(args.inprob)
            args.inlength  = os.path.abspath(args.inlength)
            args.delprob   = os.path.abspath(args.delprob)
            args.dellength = os.path.abspath(args.dellength)
            self.pos2inprob  =   make_list(args.inprob  , 'float', column = 0)
            self.in_lengths  = [ make_list(args.inlength, 'int'  , column = 0)   ,
                                 make_list(args.inlength, 'float', column = 1)   ]
            self.pos2delprob =   make_list(args.delprob , 'float', column = 0)
            self.del_lengths = [ make_list(args.dellength,'int'  , column = 0)   ,
                                 make_list(args.dellength,'float', column = 1)   ]
            self.CRISPR      = True

        # initial indel specification (for distributed computing mode)
        if(self.CRISPR):
            if (args.indels is not None):
                self.initindels = []
                with open(args.indels, 'r') as handle:
                    for line in handle:
                        chunks = line.split()
                        if (chunks[0] == "del") : self.initindels.append([chunks[0], int(chunks[1]), int(chunks[2])])
                        elif (chunks[0] == "in"): self.initindels.append([chunks[0], int(chunks[1]), int(chunks[2]), chunks[3]])
            else:
                self.initindels=[]
        else:
            self.initindels=None

        # for nwk2fa specification
        if args.tree:
            self.tree = args.tree
                
        # other specification
        self.bar = args.bar
        self.qsub = args.qsub
        self.u = args.u
        self.debug = args.debug
        self.CV = args.CV
        self.tMorigin = args.tMorigin


        # initial sequence specification
        if (args.f is not None):
            if args.f.split(".")[-1]=='gz':
                handle = gzip.open(args.f, 'rt')
            else:
                handle = open(args.f, 'r')
            sequences = SeqIO.parse(handle, 'fasta')
            self.initseq = str(list(sequences)[0].seq)
            self.L = len(self.initseq)
            handle.close()

        else:
            self.initseq = ''.join([
                np.random.choice(['A', 'G', 'C', 'T']) for i in range(args.L)
                ])

    ###########################################################################
    # some functions
    ###########################################################################
    def P(self, gamma, t):
        exp_rambda = np.diag(
                np.array([
                    np.exp(self.Al[0] * t * gamma),
                    np.exp(self.Al[1] * t * gamma),
                    np.exp(self.Al[2] * t * gamma),
                    np.exp(self.Al[3] * t * gamma)]
                )
            )
        return np.dot(np.dot(self.U, exp_rambda), np.linalg.inv(self.U))

    def sub_mat_parser(self,infile):

        # substitution matrix
        # = [
        # A: [AtoA, AtoT, AtoG, AtoC]
        # C: [CtoA, CtoT, CtoG, CtoC]
        # G: [GtoA, GtoT, GtoG, GtoC]
        # T: [TtoA, TtoT, TtoG, TtoC]
        # ] * L

        base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        with open(infile, 'r') as handle:
            
            row = 0
            for line in handle:

                if (row == 0):
                    
                    substitution_prob_matrix = np.zeros((self.L, 4, 4))
                
                else:

                    transition = line.split(",")[0] # eg. 'CT'
                    prob_list  = [float(prob_str) for prob_str in line.split(",")[1:]]

                    for pos, prob in enumerate(prob_list):

                        if (pos < self.L):
                            substitution_prob_matrix[pos][base[transition[0]]][base[transition[1]]] = prob
                        else:
                            print("args_reader.py: Edit profile assumes longer sequence than the actual sequence length")
                            sys.exit(1)

                row       += 1
        
        for pos in range(self.L):
            
            for i in range(4):

                change_prob = 0
                for j in range(4):
                    if (i != j):
                        change_prob += substitution_prob_matrix[pos, i, j] 

                if (change_prob > 1):
                    print("args_reader.py: Editing probability > 1!")
                    sys.exit(1)

                substitution_prob_matrix[pos, i, i] = 1 - change_prob
        
        #for matrix in substitution_prob_matrix:
        #    print (matrix)
        return substitution_prob_matrix