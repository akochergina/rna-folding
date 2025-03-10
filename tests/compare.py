import RNA
# The RNA sequence
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import sys
import os

parent_dir = os.path.abspath(os.path.join(os.getcwd()))
sys.path.append(parent_dir)

from modules.MEA import *
from modules.MTA import *
from modules.MEA_conservation import *

def alifold_alignment(alignment):
    # create a new model details structure
    md = RNA.md()
    # optionally one could change some parameters
    # md.temperature = 25.0 # 25 Deg Celcius
    # md.dangles = 1 # keep default 2 for compatibility with partition folding
    # create a fold compound
    fc = RNA.fold_compound(alignment, md)
    # predict the  "Alifold" Minmum Free Energy and the corresponding secondary structure
    (ss, mfe) = fc.mfe()
    conservation_score = fc.eval_covar_structure(ss)
    print("%s\n%s [ %6.2f, %6.2f ]\n" % ('\n'.join(alignment), ss, mfe, conservation_score))
    return ss


def consensus_structures(aligned_seq, constante, gamma,m):
    result=[]
    print("Running MTA...")
    result.append(MTA(aligned_seq, constante, m))
    print("Running MEA")
    result.append(MEA(aligned_seq, constante, m))
    print("Running MEA conservation")
    result.append(MEA_conservation(aligned_seq, constante, gamma, m))
    print("Running RNA Alifold")
    result.append(alifold_alignment(aligned_seq))
    return result

def similarity (struc1, struc2):
    count=0
    for i in range(len(struc1)):
        if struc1[i]=='.' and struc2[i]=='.':
            count+=1
    bp_list1=parse_RNA_structure(struc1)
    bp_list2=parse_RNA_structure(struc2)
    for tuple in bp_list1:
        if tuple in bp_list2:
            count+=2
    return count/len(struc1)


example_alignement = ["GGAGGAUUAGCUCAGCUGGGAGAGCAUCUGCCUUACAAGCAGAGGG-----------UCGGCGGUUCGAGCCCGUCAUCCUCC",
"GCCUUCCUAGCUCAG-UGGUAGAGCGCACGGCUUUUAACCGUGUGG-----------UCGUGGGUUCGAUCCCCACGGAAGGC",
"GCCUUUAUAGCUUAG-UGGUAAAGCGAUAAACUGAAGAUUUAUUUA-----------CAUGUAGUUCGAUUCUCAUUAAGGGC",
"GCGGAUAUAACUUAGGGGUUAAAGUUGCAGAUUGUGGCUCUGAAAA------------CACGGGUUCGAAUCCCGUUAUUCGC",
"GGAAAAUU-GAUCAUCGGCAAGAUAAGUUAUUUACUAAAUAAUAGGAUUUAAUAACCUGGUGAGUUCGAAUCUCACAUUUUCC"
]


#print(similarity("(((((((..((((........))))((((((.......))))))...............(((((.......))))))))))))",
#"(((((((..((((........))))(((((.........)))))...............(((((.......))))))))))))"))

def similarity_matrix(aligned_seq, constante=1, gamma=1,m=2):
    structures= consensus_structures (aligned_seq, constante, gamma,m)
    print(structures)
    matrix=np.zeros((4,4))
    for i in range (4):
        for j in range(i,4):
            sim_score=similarity(structures[i], structures[j])
            matrix[i][j]=sim_score
            matrix[j][i]=sim_score
    return(matrix)

def plot_similarity_matrix(matrix, method_names):
    plt.figure(figsize=(8, 6))
    sns.heatmap(matrix, annot=True, cmap="coolwarm", xticklabels=method_names, yticklabels=method_names)
    plt.title("Matrice de Similarité des Structures Consensus")
    plt.xlabel("Méthodes")
    plt.ylabel("Méthodes")
    plt.show()


method_names = ["MTA", "MEA", "MEA conservation", "Alifold"]
#sim_matrix = similarity_matrix(example_alignement)
#plot_similarity_matrix(sim_matrix, method_names)