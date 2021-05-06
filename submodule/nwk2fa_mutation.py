import random
import numpy as np
from Bio import Phylo, SeqIO

# nwk2fa light
class Lineage(Phylo.BaseTree.Clade):
    def __init__(self, branch_length=1.0, name=None, clades=None, confidence=None, color=None, width=None, 
    seq=None, mother_name=None, ROOT=False, parsed_args=None, indelsM=None, mother_clade=None):
        super(Lineage, self).__init__(branch_length, name, clades, confidence, color, width)

        self.mother_clade  = mother_clade
        self.mother_name   = mother_name
        self.branch_length = branch_length
        self.seq = seq if ROOT else self.mutation(seq, parsed_args)
        self.substitution_list = [] 

        if (parsed_args.CRISPR): self.indels = indelsM + self.generate_indels(parsed_args) # CRISPR == True if an inprob file path is specified 
        else                   : self.indels = None

    # receive mother SEQ sequence, introduce mutations,
    def mutation(self, seq, parsed_args):
        dseq=""
        for i in range(parsed_args.L):
            #if(parsed_args.homoplasy is not None):
            #    dseq = dseq + self.homoplastic_mutation(seq, i, True)
            if(parsed_args.constant is not None):
                dseq = dseq + self.time_independent_mutation(seq[i], parsed_args.mu[i])
            elif(parsed_args.gtrgamma is not None):
                #dseq = dseq + self.time_dependent_mutation(seq[i], parsed_args.gamma[i], self.mother_clade.branch_length, parsed_args) ### Important: equal to original PRESUME
                dseq = dseq + self.time_dependent_mutation(seq[i], parsed_args.gamma[i], self.branch_length, parsed_args) ### Important: different from original PRESUME
            elif(parsed_args.editprofile is not None):
                dseq = dseq + self.base_editing(seq[i], i, parsed_args.sub_prob_mtx_list)
        return dseq
    
    # mutation of a site (NOT Jukes Cantor model.
    # directly define mutation matrix, not the mutation rate matrix
    # it's enough for calculate mutation of each duplication
    def base_editing(self, c, position, matrix_list):

        base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        matrix      = matrix_list[position]

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
    def time_dependent_mutation(self, c, gamma, dM, parsed_args):
        base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        matrix = parsed_args.P(gamma,dM)

        # np.matrix[x] returns matrix, then the matrix is converted to array()
        return random.choices(
            ['A', 'C', 'G', 'T'], k=1, weights=np.array(matrix[base[c]])[0]
            )[0]

    def generate_indels(self,parsed_args): # pos2inprob, in_lengths, pos2delprob, del_lengths were defined from argument files

        def randomstr(alphabet, length):
            seq=""
            for _ in range(length):
                seq = seq + random.choice(alphabet)
            return seq

        seq_length = len(parsed_args.pos2inprob)
        
        shuffled_in_pos  = list(range(seq_length)); random.shuffle(shuffled_in_pos)
        shuffled_del_pos = list(range(seq_length)); random.shuffle(shuffled_del_pos)
        indel_order      = ['in', 'del'];           random.shuffle(indel_order)

        generated_indels = []

        for indel in indel_order:

            if (indel == 'in'):

                for pos in shuffled_in_pos:

                    if ( random.random() < parsed_args.pos2inprob[pos] ):

                        length = random.choices(parsed_args.in_lengths[0], k = 1, weights = parsed_args.in_lengths[1])[0]
                        
                        generated_indels.append( [ 'in',  pos, length, randomstr(['A','T','G','C'], length) ] )
            
            elif (indel == 'del' ):

                for pos in shuffled_del_pos:
                    
                    if ( random.random() < parsed_args.pos2delprob[pos] ):
                        
                        length = random.choices(parsed_args.del_lengths[0], k = 1, weights = parsed_args.del_lengths[1])[0]

                        generated_indels.append( [ 'del', pos, length ] )
            
        return generated_indels
    
    # return a sequence with indel
    def get_seq_with_indel(self): # seqstr: str
        seqstr          = self.seq

        indel_list      = self.indels
        refpos2pos      = {i:i for i in range(len(seqstr))}
        initL           = len(seqstr)
        true_indel_list = []
        pos2refpos      = list(range(len(seqstr)))
        dropout         = False

        for indel in indel_list:
            if (len(seqstr)==0):
                dropout=True
                break

            if (indel[0] == 'del'):

                if( indel[1] in refpos2pos.keys() ):
                    true_indel_list.append(indel)
                    mid    = refpos2pos[indel[1]] # pos is the midpoint of deletion
                    length = indel[2]
                    start  = max ( 0, mid - length//2 )
                    end    = min ( len(seqstr) - 1, mid - length//2 + length - 1) 
                    seqstr     = seqstr[:start]     + seqstr[(end+1):]
                    pos2refpos = pos2refpos[:start] + pos2refpos[(end+1):]

                    if (len(seqstr) != len(pos2refpos)):
                        print("len(seqstr) != len(pos2refpos)!")
                        return
                        
                    refpos2pos = {}
                    for pos, refposindel in enumerate(pos2refpos):
                        if(type(refposindel)==int):
                            refpos2pos[refposindel] = pos
                        elif(refposindel=="in"):
                            None
            elif ( indel[0] == "in" ):
                if( indel[1] in refpos2pos.keys() ):
                    true_indel_list.append(indel)
                    start      = refpos2pos[indel[1]]
                    length     = indel[2]
                    seqstr     = seqstr[:start+1]     + indel[3]      + seqstr[start+1:]         
                    pos2refpos = pos2refpos[:start+1] + ["in"]*length + pos2refpos[start+1:]
                    refpos2pos = {}
                    for pos, refposindel in enumerate(pos2refpos):
                        if(type(refposindel)==int):
                            refpos2pos[refposindel] = pos
                        elif(refposindel=="in"):
                            None
            
        
        if ( dropout ):
            aligned_seqstr = "-"*initL
        else:
            aligned_seqstr = "-"*initL
            for refpos in refpos2pos.keys():
                aligned_seqstr = aligned_seqstr[:refpos] + seqstr[refpos2pos[refpos]] + aligned_seqstr[(refpos+1):]
        
        #print(seqstr, true_indel_list)
        return seqstr, aligned_seqstr, true_indel_list 