# standard packages
import os
import subprocess

# third party packages
import numpy as np

# geobind packages
from .io_utils import moveFile
from .mesh import Mesh

def runMSMS(atoms, file_prefix='mesh', basedir='.', 
        clean=True, quiet=True, hydrogens=True, area_only=False, mesh_kwargs={}, **kwargs
    ):
    # generate coordinate file
    if not isinstance(atoms, list):
        atoms = atoms.get_atoms()
    coordFile = "{}_coords.xyzr".format(file_prefix)
    FH = open(coordFile, "w")
    for atom in atoms:
        if (not hydrogens) and atom.element == "H":
            continue
        atmn = atom.name
        acoords = atom.get_coord()
        radius = atom.xtra["radius"]
        FH.write("{:7.4f} {:7.4f} {:7.4f} {:3.2f}\n".format(acoords[0], acoords[1], acoords[2], radius))
    FH.close()
    
    # set MSMS options
    msms_opts = {
        'probe_radius': 1.5,
        'density': 1.0,
        'hdensity': 3.0,
        'surface': 'tses'
    }
    msms_opts.update(kwargs)
    if area_only:
        msms_opts['surface'] = 'ases'
    
    # run MSMS and generate vertex and face file
    args = [
        "msms",
        "-probe_radius", str(msms_opts['probe_radius']),
        "-density", str(msms_opts['density']),
        "-hdensity", str(msms_opts['hdensity']),
        "-if", coordFile,
        "-of", file_prefix,
        "-af", "{}.area".format(file_prefix),
        "-surface",  msms_opts["surface"]
    ]
    if quiet:
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
    
    if area_only:
        # delete/move files and return path to the area file
        if clean:
            os.remove("{}_coords.xyzr".format(file_prefix))
        else:
            moveFile("{}_coords.xyzr".format(file_prefix), basedir)
        af = "{}.area".format(file_prefix)
        moveFile(af, basedir)
        
        return os.path.join(os.path.abspath(basedir), af)
    else:
        # Get vertices
        vertData = open("{}.vert".format(file_prefix)).readlines()[2:]
        vertexs = []
        normals = []
        for line in vertData[1:]:
            line = line.strip().split()
            vertex = np.array(line[0:3], dtype=np.float32)
            normal = np.array(line[3:6], dtype=np.float32)
            
            vertexs.append(vertex)
            normals.append(normal)
        vertexs = np.array(vertexs, dtype=np.float32)
        normals = np.array(normals, dtype=np.float32)
        
        # Get faces
        faceData = open("{}.face".format(file_prefix)).readlines()[3:]
        faces = []
        for line in faceData:
            line = line.strip().split()
            i = int(line[0])-1
            j = int(line[1])-1
            k = int(line[2])-1
            faces.append([i,j,k])
        faces = np.array(faces, dtype=np.int32)
        
        # clean up MSMS files
        if clean:
            os.remove("{}.face".format(file_prefix))
            os.remove("{}.vert".format(file_prefix))
            os.remove("{}_coords.xyzr".format(file_prefix))
            os.remove("{}.area".format(file_prefix))
        else:
            moveFile("{}.face".format(file_prefix), basedir)
            moveFile("{}.vert".format(file_prefix), basedir)
            moveFile("{}_coords.xyzr".format(file_prefix), basedir)
            moveFile("{}.area".format(file_prefix), basedir)
        
        # return mesh
        return Mesh(vertices=vertexs, faces=faces, vertex_normals=normals, name=file_prefix, **mesh_kwargs)
