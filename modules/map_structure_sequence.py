"""
This module is used to align structures on aligned sequences
"""

def get_rid_brackets_op(list_of_car, index_opening):
    """deletes in the list_of_car the opening bracket at index index_opening and replace associated closing bracket by "."
    returns the new list and the indix of the closing bracket
    """
    list_of_car.pop(index_opening)#erases opening bracket
    index=0
    for i in range(index_opening,len(list_of_car)):
        if list_of_car[i]=="(":
            index+=1
        elif list_of_car[i] == ")":
            if index==0:
                list_of_car[i]="."
                return list_of_car
            else :
                index-=1
    raise NameError('No associated closing bracket. Invalid dbstring')



def get_rid_brackets_cl(list_of_car, index_closing):
    """deletes in the list_of_car the closing bracket at index index_closing and replace associated opening bracket by "."
    returns the new list and the indix of the opening bracket
    """
    list_of_car.pop(index_closing)#erases opening bracket
    index=0
    for i in range(index_closing,-1,-1): #go back in the list to find the opening bracket
        if list_of_car[i]==")":
            index+=1
        elif list_of_car[i] == "(":
            if index==0:
                list_of_car[i]="."
                return list_of_car
            else :
                index-=1
    raise NameError('No associated opening bracket. Invalid dbstring')



def projection(alignmentstr,dbstring):
    """
    returns the projected sequence and the projected dot bracket string of the alignement alignmentstr

    Args:
      alignmentstr : string representing the aligned sequence
      dbstring : dbstring of the consensus sequence
    
    Returns:
      proj_seq : string of the projected sequence
      proj_db: Sstring of the projected dot bracket string
    """
    proj_seq=[]
    proj_db=list(dbstring)
    indexdbstring=0 #we need to index how we move in the dbstring because removing prenthesis will change the index
    for i in range (len(alignmentstr)):
        if alignmentstr[i]=="-":
            if proj_db[indexdbstring]=="(":
                proj_db=get_rid_brackets_op(proj_db,indexdbstring)
            elif proj_db[indexdbstring]==")":
                proj_db=get_rid_brackets_cl(proj_db,indexdbstring)
            else :
                proj_db.pop(indexdbstring)
        else :
            proj_seq.append(alignmentstr[i])
            indexdbstring+=1
    return ''.join(proj_seq), ''.join(proj_db)



def reverse_projection(alignementstr, dbstr):
    """
    returns the alignend structure (dbstr) of the aligned sequence
    ex : "C-AG, (.) -> (-.)"

    Args:
      alignmentstr : string representing the aligned sequence
      dbstring : dbstring of the core sequence (projected sequence)
    
    Returns:
      align_struc : the aligned structre
    """
    align_strc=[]
    count=0
    for car in alignementstr:
        if car=="-":
            align_strc.append("-")
        else:
            align_strc.append(dbstr[count])
            count+=1
    return ''.join(align_strc)