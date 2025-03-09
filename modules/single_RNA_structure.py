"""
This module implements simple nussinov and nussinov traceback 
to construct the structure of a single rna sequence
"""

complementary_bases=["AT","CG","AU","UG"]

def nussinov_matrix(alignstr, cost_function, m):
    """
    Args : alignstr : string of the sequence e.g. "AG-UGGUA"
        const_function function -> car, car -> int returns the cost of the distance between two bases
        m : int, minimal loop length
    Returns : matrix (list of list) of nussinov:
            where M[i][j] = max number of base pair in the sequence between base i and base j
    """
    n=len(alignstr)
    list_of_bases=list(alignstr.upper())
    nussinovMatrix = [[0 for i in range (n)]for j in range(n)]

    #complete matrix
    for k in range(m, n):
        i=0
        j=i+k
        while (i<n and j<n):
            if (i+1<n and j-1>=0):
                alpha0= nussinovMatrix[i+1][j-1] + cost_function(list_of_bases[i],list_of_bases[j])
            else :
                alpha0=0
            alpha1 = max(
                nussinovMatrix[i][k-1] + nussinovMatrix[k][j]
                for k in range(i, j+1)
            )
            nussinovMatrix[i][j]=max(alpha0, alpha1)
            i+=1
            j+=1
    return nussinovMatrix

def nussinov_traceback_rec(alignstr,nussMat, costFun, indices,m):
    '''returns the list of aligned base pairs from i to j 
    args : alignstr : The rna string
            nussMat : the nussinov matrix
            indices : tuple (i,j)
            costFun : the base cost function
            m = the minimum loop
    Returns : list of tuples corresponding to aligned strings between i and j
    '''
    previous_bases=[]
    i,j=indices
    #print("Looking at",i,j)
    list_of_bases=list(alignstr.upper())
    if (j-i)<=m:
        return []
    else :
        n=len(nussMat[i])
        cost=nussMat[i][j]
        if (i+1<n and j-1>=0):
            alpha0= nussMat[i+1][j-1] + costFun(list_of_bases[i],list_of_bases[j])
            if cost==alpha0:
                previous_bases=nussinov_traceback_rec(alignstr,nussMat,costFun, (i+1,j-1),m)
                if list_of_bases[i]+list_of_bases[j] in complementary_bases or list_of_bases[j]+list_of_bases[i] in complementary_bases :
                    previous_bases.append((i,j))
                return previous_bases
        for k in range(i+1,j):
            if cost== (nussMat[i][k-1] + nussMat[k][j]) :
                previous_bases=nussinov_traceback_rec(alignstr,nussMat,costFun, (i,k-1),m) + nussinov_traceback_rec(alignstr,nussMat,costFun, (k,j),m)
                return previous_bases
        print("Not Found precedent",i,j)
        return []
    
def nussinov_TB(alignstr,nussMat, costFun,m):
    """
    The Nussinov Traceback algorithm
    """
    n=len(nussMat)
    return nussinov_traceback_rec(alignstr,nussMat, costFun, (0,n-1),m)

def RNA_consensus_structure(alignstr, cost_function, m=2):
    """
    takes the rna string and returns the consensus structure
    Args : alignstr : string of the sequence e.g. "AG-UGGUA"
        const_function function -> car, car -> int returns the cost of the distance between two bases
        m : int, minimal loop length
    Returns : list of base pairs
    """
    nussMat=nussinov_matrix(alignstr, cost_function, m)
    return nussinov_TB(alignstr,nussMat, cost_function,m)