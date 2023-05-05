# built-in modules
import os
import logging
import subprocess

# third party modules
import numpy as np

# geobind modules
from .io_utils import moveFile

def stripChainID(pqr_file, prefix=None):
    if prefix is None:
        prefix = pqr_file.split('.')[0] + "_nochain"
    
    new_file_name = "%s.pqr" % prefix
    new_file = open(new_file_name, "w")
    
    # remove chain ID from PQR file 
    for line in open(pqr_file):
        if line[0:4] == "ATOM":
            line = list(line)
            line[21] = ""
            line = "".join(line)
        new_file.write(line)
    new_file.close()
    
    return new_file_name

def parseVTK(file_name, keep_vtk=True):
    if isinstance(file_name, str):
        file_handle = open(file_name)
    
    lines = file_handle.readlines()[5:] # skip header
    file_handle.close()
    
    _, npoints, __ = lines[0].strip().split()
    npoints = int(npoints)
    
    # get coordinate info
    sptr = 1
    eptr = sptr + npoints
    points = np.array([_.strip().split() for _ in lines[sptr:eptr]], dtype=np.float32)
    
    # get face info
    _, nfaces, __  = lines[eptr].strip().split()
    nfaces = int(nfaces)
    sptr = eptr + 1
    eptr = sptr + nfaces
    faces = np.array([_.strip().split()[1:] for _ in lines[sptr:eptr]], dtype=int)
    
    # get point data
    sptr = eptr + 4
    eptr = sptr + npoints
    pdata = np.array(lines[sptr:eptr], dtype=np.float32)
        
    # get normal derivative
    sptr = eptr + 2
    eptr = sptr + npoints
    ndata = np.array(lines[sptr:eptr], dtype=np.float32)
    
    # get normals
    sptr = eptr + 2
    eptr = sptr + npoints
    normals = np.array([_.strip().split() for _ in lines[sptr:eptr]], dtype=np.float32)
    
    assert npoints == len(points) == len(pdata) == len(ndata)
    assert nfaces == len(faces)
    
    if not keep_vtk:
        os.remove(file_name)
    
    return points, faces, normals, pdata, ndata

def runTABIPB(pqr, command="tabipb", prefix=None, basedir='.', mesh_type="skin", keep_vtk=True, keep_pqr=True, clean=True, strip_chain=True):
    if strip_chain:
        pqr = stripChainID(pqr)
        keep_pqr = False
    
    # run TABI-PB
    input_file = """mol               {}
mesh              {}
sdens             1
srad              1.4
pdie              2
sdie              78.0
bulk              0.15
temp              310.0
tree_degree       2
tree_max_per_leaf 50
tree_theta        0.8
precondition      OFF
outdata           VTK""".format(
        pqr,
        mesh_type,
    )
    inFile = os.path.join(basedir, "{}.in".format(prefix))
    FH = open(inFile, "w")
    FH.write(input_file)
    FH.close()
    
    logging.info("Running TABI-PB on input file: %s", inFile)
    outpt = subprocess.check_output([command, inFile], stderr=subprocess.STDOUT)
    
    fname = "output.vtk"
    if os.path.exists(fname):
        if prefix is not None:
            os.rename("output.vtk", "%s.vtk" % prefix)
            fname = "%s.vtk" % prefix
        
        if basedir != ".":
            fname = moveFile(fname, basedir)
    
    if clean:
        os.remove(inFile)
    
    if not keep_pqr:
        os.remove(pqr)
    
    return parseVTK(fname, keep_vtk)

