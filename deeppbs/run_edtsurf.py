# builtin modules
import os
import subprocess

# geobind modules
from geobind_dna import moveFile
from .mesh import Mesh

def runEDTSurf(pdbfile, file_prefix='mesh', basedir='.', clean=True, quiet=True, mesh_kwargs={}, **kwargs):
    # run EDTSurf and generate .PLY file
    edtsurf_args = {
        "-t": "2",
        "-s": "3",
        "-c": "1",
        "-p": "1.4",
        "-f": "2.0",
        "-h": "2",
        "-o": file_prefix
    }
    edtsurf_args.update(kwargs)
    
    args = [
        "EDTSurf",
        "-i",
        pdbfile
    ]
    for key in edtsurf_args:
        args.append(key)
        args.append(edtsurf_args[key])
    if(quiet):
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
    
    meshfile = "{}.ply".format(file_prefix)
    mesh = Mesh(handle=meshfile, name=file_prefix, **mesh_kwargs)
    # clean up EDTSurf files
    files = [
        "{}-cav.pdb".format(file_prefix),
        meshfile,
        pdbfile
    ]
    if clean:
        for f in files:
            if os.path.exists(f):
                os.remove(f)
    else:
        for f in files:
            moveFile(f, basedir)
    
    return mesh
    
