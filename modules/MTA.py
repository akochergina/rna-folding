"""
This module implements the Maximum Total Accuracy
It uses a variant of the nussinov algorithm alrady used to construct the single rna structure
"""

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.rna_structure2 import *
from modules.rna_sequence2 import *
from modules.single_RNA_structure import *
from modules.map_structure_sequence import *
from modules.simplecount import *

def unpaired_acc(i,list_of_struc):
    ''' Computes the cost function unpaired_acc
    Arg : list_of_struc : a list of structures, defined as a dot brackets
        i the int we are computing the cost of
    Returns : int : unpaired_acc(i)
    '''
    count=0
    for struc in list_of_struc:
        if struc[i]=='.':
            count+=1
    return count

def paired_acc(indices, list_of_struc, constante):
    ''' Computes the cost function paired_acc
    Arg : list_of_struc : a list of structures, defined as dot brackets
        indices : the tuple of indices we are looking for the cost of
    Returns : int : paired_act(i)
    '''
    count=0
    for struc in list_of_struc:
        #print(struc)
        base_pairs=parse_RNA_structure(struc)
        if indices in base_pairs:
            count+=1
    return 2*constante*count

def constfunctionMTA(constante):
    return lambda indices,list_of_struc : paired_acc(indices, list_of_struc, constante)



def MTA_nussinov_matrix(list_of_structures, cost_function, init_cost, m=2):
    """
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of indices, list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : a cost fonction to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    n=len(list_of_structures[0])
    nussinovMatrix = [[0 for i in range (n)]for j in range(n)]

    #initialize for all i-j<=m
    for k in range(m+1):
        i=0
        j=i+k
        while (i<n and j<n):
            if j==i:
                nussinovMatrix[i][j]=init_cost(i, list_of_structures)
            else :
                nussinovMatrix[i][j]=nussinovMatrix[i][j-1]+init_cost(j, list_of_structures)
            i+=1
            j+=1

    #complete matrix
    for k in range(m, n):
        i=0
        j=i+k
        while (i<n and j<n):
            if (i+1<n and j-1>=0):
                alpha0= nussinovMatrix[i+1][j-1] + cost_function((i,j), list_of_structures)
            else :
                alpha0=0
            alpha1 = max(nussinovMatrix[i][k-1] + nussinovMatrix[k][j]
                         for k in range(i, j+1))
            nussinovMatrix[i][j]=max(alpha0, alpha1)
            i+=1
            j+=1
    return nussinovMatrix


def nussinov_MTA_traceback_rec(list_of_structures,nussMat, costFun, initCost, indices,m):
    """
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of indices, list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : a cost fonction to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    previous_bases=[]
    i,j=indices
    found=False
    if j-i<=m:
        return []
    else :
        n=len(nussMat[i])
        cost=nussMat[i][j]
        if (i+1<n and j-1>=0):
            alpha0= nussMat[i+1][j-1] + costFun(indices, list_of_structures)
            if cost==alpha0:
                previous_bases=nussinov_MTA_traceback_rec(list_of_structures,nussMat, costFun, initCost,(i+1,j-1),m)
                previous_bases.append((i,j))
                found=True
                #print("happens")
        if not found:
            for k in range(i+1,j+1):
                if cost== (nussMat[i][k-1] + nussMat[k][j]) and not found:
                    found=True
                    previous_bases=nussinov_MTA_traceback_rec(list_of_structures,nussMat, costFun, initCost,(i,k-1),m) + nussinov_MTA_traceback_rec(list_of_structures,nussMat, costFun, initCost,(k,j),m)
        return previous_bases
    

def nussinov_TB_MTA(list_of_structures,nussMat, costFun,initCost,m):
    """"
    The traceback algorithm for the nussinov matrix
    Args : list_of_structures : list of the aligned structures
        const_function function -> tuple of indices, list of structures -> int returns the cost of adding a base pair ij with regard to the list of aligned structures
        init_cost : a cost fonction to initialize the matrix
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    n=len(nussMat)
    return nussinov_MTA_traceback_rec(list_of_structures,nussMat, costFun, initCost, (0,n-1),m)



def maximum_total_accuracy_db_str(list_of_structures,m):
    """ Returns the consensus dbstring of the structure 
    of the list of the already aligned structures
    Args : list_of_structures : list of the aligned structures
            m : int, minimal loop length
    Returns : string : the consensus structure in form of dot bracket string
    """
    nussMat=MTA_nussinov_matrix(list_of_structures, constfunctionMTA(1), unpaired_acc, m)
    #print(nussMat)
    base_pairs= nussinov_TB_MTA(list_of_structures,nussMat, constfunctionMTA(1),unpaired_acc,m)
    return (dot_bracket_string(len(list_of_structures[0]),base_pairs,0))



def MTA(aligned_strings_list,m=2):
    #reconstruc the list of unaligned RNA sequences
    list_of_seq=[]
    for algstr in aligned_strings_list:
        list_of_seq.append(core_sequence(algstr))
    #Reconstruct each structure of each unaligned sequence, then realign the structure
    list_of_struc=[]
    i=0
    for sequence in list_of_seq:
        bp=(RNA_consensus_structure(sequence, simple_count_cost, m)) #the base pairs of the structure from the unaligned sequence
        db=dot_bracket_string(len(sequence),bp,0)#the dot-bracket string of the structure from the unaligned sequence
        print(aligned_strings_list[i], db)
        aligned_db=reverse_projection(aligned_strings_list[i], db) #the aligned dot bracket
        list_of_struc.append(aligned_db)
        i+=1
    return maximum_total_accuracy_db_str(list_of_struc,m)

#print(MTA (["-GCC-AAA-GGC","GGGC-AUU-GCC", "-ACGGAAUCCGU"]))
"""

test=["GCAGUCGUGGCCGAGU---GGUUAAGGCGUCUGACUCGAAAUCAGAUUCCCUCUGGGAGCGUAGGUUCGAAUCCUACCGGCUGCG",
"GCGGGGGUGCCCGAGCCUGGCCAAAGGGGUCGGGCUCAGGACCCGAUGGCGUAGGCCUGCGUGGGUUCAAAUCCCACCCCCCGCA",
"UGGAGUAUAGCCAAG--UGG--UAAGGCAUCGGUUUUUGGUACCG---------GCAUGCAAAGGUUCGAAUCCUUUUACUCCAG",
"CGGAAAGUAGCUUAGCUUGG--UAGAGCACUCGGUUUGGGACCGA---------GGGGUCGCAGGUUCGAAUCCUGUCUUUCCGA",
"GCCGGGGUGGGGUAGUGGCCAUCCUGG---GGGACUGUGGAUCCC----------CUGACCCGGGUUCAAUUCCCGGUCCCGGCC",
"GUAAACAUAGUUUA------AUCAAAACAUUAGAUUGUGAAUCUAA----------CAAUAGAGGCUCGAAACCUCUUGCUUACC",
"AGUAAAGUCAGCUA------AAAAAGCUUUUGGGCCCAUACCCCAA----------ACAUGUUGGUUAAACCCCUUCCUUUACUA"]

print(MTA (test))"""