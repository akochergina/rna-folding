"""
This module concerns all functions that processes RNA structures
dbstring, parse_RNA
"""

if __name__ == "__main__":
    print("Module rna_structure2 importé avec succès !")


def dot_bracket_string(seq_length, base_pairs, decalage=-1):
    """
    returns the dot bracket string of the structure defined by the list of base pairs

    Args:
      seq_length : length of the sequence
      base_pairs : set or list of base pairs
    
    Returns:
      String that represents the dot bracket of the structure
    """
    db=["." for i in range (seq_length+1+decalage)]
    for pair in base_pairs:
        ind0, ind1=pair
        db[ind0+decalage]="("
        db[ind1+decalage]=")"
    db=''.join(db)
    return db

def parse_RNA_structure(dbstring):
    """
    returns the list of base pairs of the structure defined by the dot bracket string

    Args:
      dbstring : String that represents the dot bracket of the structure
    
    Returns:
      A list of tuples that represent the base pairs
    """
    base_pairs=[]
    list_db=list(dbstring)
    n=len(list_db)
    for i in range(len(list_db)):
        if dbstring[i] == "(":
            #find associated closing parenthesis
            index=0
            j=i+1
            found=False
            while j<n and not found:
                if list_db[j]=="(":
                    index+=1
                elif list_db[j] == ")":
                    if index==0:
                        found=True
                        base_pairs.append((i,j))  
                    else :
                        index-=1
                j+=1
            if not found :
                raise NameError('Invalid dbstring')
    return base_pairs
