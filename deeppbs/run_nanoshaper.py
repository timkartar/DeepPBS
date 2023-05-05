# builtin modules
import os
import subprocess

# geobind modules
from deeppbs import moveFile
from .mesh import Mesh

def runNanoShaper(atoms, file_prefix, basedir, 
        clean=True, quiet=True, hydrogens=True, pockets_only=False, mesh_kwargs={}, **kwargs
    ):
    # generate coordinate file
    if(not isinstance(atoms, list)):
        atoms = atoms.get_atoms()
    coordFile = "{}.xyzr".format(file_prefix)
    FH = open(coordFile, "w")
    for atom in atoms:
        if((not hydrogens) and atom.element == "H"):
            continue
        atmn = atom.name
        acoords = atom.get_coord()
        radius = atom.xtra["radius"]
        FH.write("{:7.4f} {:7.4f} {:7.4f} {:3.2f} {}\n".format(acoords[0], acoords[1], acoords[2], radius, atom.serial_number))
    FH.close()
    
    # run NanoShaper and generate .OFF file
    nanoshaper_args = {
        "grid_scale": 2.0,
        "grid_perfil": 90.0,
        "op_mode": "normal",
        "surface_type": "skin",
        "build_status_map": "true",
        "cavity_filling_volume": 100,
        "skin_surface_parameter": 0.45,
        "accurate_triangulation": "true",
        "smooth_mesh": "true",
        "blobbyness": -2.5,
        "radius_big": 3.0
    }
    nanoshaper_args.update(kwargs)
    
    if pockets_only:
        # change operation mode
        nanoshaper_args['op_mode'] = 'pockets'
    prm_template = """
# Global Parameters
Operative_Mode = {op_mode}
Grid_scale = {grid_scale}
Grid_perfil = {grid_perfil}
Number_thread = 32 
# Map Settings
Build_epsilon_maps = false
Build_status_map = {build_status_map}
# Surface Parameters
Surface = {surface_type}
Skin_Surface_Parameter = {skin_surface_parameter}
Blobbyness = {blobbyness}
Skin_Fast_Projection = false
Accurate_Triangulation = {accurate_triangulation}
Triangulation = true
Check_duplicated_vertices = true
Smooth_Mesh = {smooth_mesh}
# Pocket/Cavities settings
Pockets_And_Cavities = false
Cavity_Detection_Filling = true
Keep_Water_Shaped_Cavities = false
Conditional_Volume_Filling_Value = {cavity_filling_volume}
Num_Wat_Pocket = 2
Pocket_Radius_Big = {radius_big}
Pocket_Radius_Small = 1.4
# I/O Settings
XYZR_FileName = {file_prefix}.xyzr
Save_Cavities = false
Save_Status_map = false
Vertex_Atom_Info = false"""

    PRM = open("{}.prm".format(file_prefix), 'w')
    PRM.write(prm_template.format(file_prefix=file_prefix, **nanoshaper_args))
    PRM.close()
    args = [
        "NanoShaper",
        "{}.prm".format(file_prefix)
    ]
    if quiet:
        FNULL = open(os.devnull, 'w')
        subprocess.call(args, stdout=FNULL, stderr=FNULL)
        FNULL.close()
    else:
        subprocess.call(args)
    
    if pockets_only:
        # Do not return a mesh, exit from here
        if clean:
            os.remove("{}.prm".format(file_prefix))
            os.remove("{}.xyzr".format(file_prefix))
        return
    
    # rename mesh file
    if not os.path.exists("triangulatedSurf.off"):
        raise FileNotFoundError("NanoShaper did not produce any output! Check that your settings are valid...")
    os.rename("triangulatedSurf.off", "{}.off".format(file_prefix))
    
    # clean up NanoShaper files
    if clean:
        os.remove("{}.prm".format(file_prefix))
        os.remove("{}.xyzr".format(file_prefix))
    else:
        moveFile("{}.prm".format(file_prefix), basedir)
        moveFile("{}.xyzr".format(file_prefix), basedir)
    
    # move meshfile to base dir
    meshfile = "{}.off".format(file_prefix)
    if basedir != os.getcwd():
        moveFile(meshfile, basedir)
        meshfile = os.path.join(basedir, meshfile)
    
    return Mesh(handle=meshfile, name=file_prefix, **mesh_kwargs)
