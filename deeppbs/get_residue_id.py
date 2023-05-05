def getResidueID(residue):
    cid = residue.get_parent().get_id()
    _, num, ins = residue.get_id()
    
    if(len(ins) == 0):
        ins = ' '
    
    return "{}.{}.{}".format(cid, num, ins)
