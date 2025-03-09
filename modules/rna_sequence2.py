"""
This module concerns all functions that processes RNA sequences
"""

def core_sequence(algstr:str) -> str:
    """
    Core sequence of an aligned RNA (ex: C-AT -> CAT)

    Args:
      algstr : alignment string of sequence 
    Returns:
      core sequence
    """

    algstr=list(algstr.upper())
    result=[]
    for i in range (len(algstr)):
        if algstr[i]!='-':
            result.append(algstr[i])
    return ''.join(result)

def sequence_identity(xp:str, yp:str) -> str:
    """
    Sequence identity of two aligned sequences

    Args:
      xp: alignment string of sequence x
      yp: alignment string of sequence y
    
    Returns:
      Percentage of identical alignment columns
    """
    identical_count=0
    #transform str to lists so we can modify them
    xp=list(xp.upper())
    yp=list(yp.upper())
    for i in range (len(xp)):
        if xp[i]=='T':
            xp[i]='U'
        if yp[i]=='T':
            yp[i]='U'
        if xp[i]==yp[i]:
            identical_count+=1
    #transform lists back to str
    xp= ''.join(xp)
    yp = ''.join(yp)
    #print(xp)
    return identical_count/len(xp)

def edges_alignement(xp:str, yp:str) -> list:
    """
    returns the list of the edges in an alignment

    Args:
      xp: alignment string of sequence x
      yp: alignment string of sequence y
    
    Returns:
      List of edges, an edge being a tuple of the indexes of the aligned bases
    """

    indexx=1
    indexy=1
    edges=[]
    for i in range(len(xp)):
        if xp[i]!="-":
            if yp[i]!="-":
                edges.append((indexx, indexy))
                indexx+=1
                indexy+=1
            else :
                indexx+=1
        else :
            if yp[i]!="-":
                indexy+=1
    return edges
