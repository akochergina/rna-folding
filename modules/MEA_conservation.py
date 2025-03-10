import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.simplecount import simple_count_cost
from modules.MEA import *
from modules.MTA import *

def gamma_cov(indices, aligned_sequences):
    i,j=indices
    count=0
    for seq1 in aligned_sequences:
        for seq2 in aligned_sequences:
            if simple_count_cost(seq1[i], seq1[j]) and simple_count_cost(seq2[i], seq2[j]):
                count+= int(seq1[i]!=seq2[i])
                count+= int(seq1[j]!=seq2[j])
    return count/2

def gamma_inc(indices, aligned_sequences, gamma):
    count=0
    i,j=indices
    for seq in aligned_sequences:
        if seq[i]=="-" and seq[j]=="-":
            count-=0.25
        elif not simple_count_cost(seq[i], seq[j]):
            count-=1

    return gamma*count

def constfunction_conservation(constante, list_of_bpp, aligned_sequences, gamma):
    return lambda indices, list_of_struc : paired_expected_acc(indices, list_of_bpp, constante) + gamma_cov(indices, aligned_sequences) + gamma_inc(indices, aligned_sequences, gamma)

#def constfunction_conservation(list_of_bp, list_of_bpp, aligned_sequences, gamma, constante):
#    return lambda indices,list_of_struc : constfunctionMTA(list_of_bp, constante) + gamma_cov(indices, aligned_sequences) + gamma_inc(indices, aligned_sequences, gamma)

def MEA_conservation(list_of_seq,constante, gamma,m=2):
    print(list_of_seq)
    list_of_bpp, list_of_struc=list_of_bp_proba_and_struc(list_of_seq)
    print(list_of_struc)
    cost_fun=constfunction_conservation(constante, list_of_bpp, list_of_seq, gamma)
    init_cost=initcostMEA(list_of_bpp)
    nussMat=MTA_nussinov_matrix(list_of_struc, cost_fun, init_cost, m)
    print(nussMat)
    base_pairs= nussinov_TB_MTA(list_of_struc,nussMat, cost_fun, init_cost, m)
    return (dot_bracket_string(len(list_of_struc[0]),base_pairs,0))
