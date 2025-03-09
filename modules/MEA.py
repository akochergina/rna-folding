import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.MTA import *
import RNA

def list_of_bp_proba_and_struc(list_of_sequences):
    '''Computes, for each structure, the matrix of base pair probabilities. 
    Will be useful then to compute unpaired_expected_accuracy and paired_expected accuracy.
    Arg : list_of_sequences : a list of sequences
    Returns : list of the bpp matrixes
                    list_bpp[i] is the base pair proba matrix of sequence list_of_seq[i]
            list_of_struc a list of the associated structues
    '''
    #list of the matrixes
    bpp_list=[]
    #list of the structures
    list_of_struc=[]
    # create model details
    md = RNA.md()
    # activate unique multibranch loop decomposition
    md.uniq_ML = 1
    for struc in list_of_sequences:
        # create fold compound object
        fc = RNA.fold_compound(struc, md)
        # compute MFE
        (ss, mfe) = fc.mfe()
        # rescale Boltzmann factors according to MFE; rescaling avoids numerical problems for long sequences
        fc.exp_params_rescale(mfe)
        # compute partition function to fill DP matrices
        fc.pf()
        bpp = fc.bpp()
        bpp_list.append(bpp)
        list_of_struc.append(ss)
    return bpp_list, list_of_struc

def unpaired_expected_acc(i, list_of_bpp):
    ''' Computes the cost function unpaired_expected_acc for the MEA
    Arg :
        i the int we are computing the cost of
        list_of_bpp : list of the bpp matrixes associated to the structure
    Returns : float : unpaired_expected_acc(i) 
    '''
    #print(i)
    somme=0
    for k in range (len(list_of_bpp)):
        paired_proba=0
        for j in range (len(list_of_bpp[0])):
            #print(i,j,k)
            #print(list_of_bpp[k])
            paired_proba+=list_of_bpp[k][j][i]
            paired_proba+=list_of_bpp[k][i][j]
        somme+= (1-paired_proba)
    return somme
    

def paired_expected_acc(indices,list_of_bpp, constante):
    ''' Computes the cost function paired_acc
    Arg : list_of_struc : a list of structures, defined as a list of tuples of paired bases
        indices : the tuple of indices we are looking for the cost of
    Returns : int : paired_act(i)
    '''
    i,j=indices
    sum=0
    for k in range(len(list_of_bpp)):
        sum+=list_of_bpp[k][i][j]
        sum+=list_of_bpp[k][j][i] #only one of them is non null
    return 2*constante*sum

def constfunctionMEA(constante, list_of_bpp):
    return lambda indices, list_of_struc : paired_expected_acc(indices, list_of_bpp, constante)

def initcostMEA(list_of_bpp):
    return lambda i, list_of_struc : unpaired_expected_acc(i, list_of_bpp)


def MEA(list_of_seq, constante,m=2):
    print(list_of_seq)
    list_of_bpp, list_of_struc=list_of_bp_proba_and_struc(list_of_seq)
    print(list_of_struc)
    cost_fun=constfunctionMEA(constante, list_of_bpp)
    init_cost=initcostMEA(list_of_bpp)
    nussMat=MTA_nussinov_matrix(list_of_struc, cost_fun, init_cost, m)
    print(nussMat)
    base_pairs= nussinov_TB_MTA(list_of_struc,nussMat, cost_fun, init_cost, m)
    return (dot_bracket_string(len(list_of_struc[0]),base_pairs,0))

#print (MEA (["-GCC-AAA-GGC","GGGC-AUU-GCC", "-ACGGAAUCCGU"], 1))