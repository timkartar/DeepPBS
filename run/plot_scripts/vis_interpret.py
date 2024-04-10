from pymol import cmd
import numpy as np
#import pip
#pip.main(['install','matplotlib'])
import matplotlib.cm as cm
import sys
#from arrow import *
#pdbpath = "./"
#pdbpath = "../../test_md_pwms/2r5z_protonated_crystal_BSC1/"

#pdbpath = "/home/raktim/cameron/designed/requestionaboutmodelforpredictingpwmsfromcrystalst/splitted/"
fpath = sys.argv[2]
pdbpath = fpath + "/../process/pdb/"
pdb = sys.argv[3]
if pdb.endswith("pdb"):
    npz = pdb.replace("pdb","npz")
elif pdb.endswith("cif"):
    npz = pdb.replace("cif","npz")
else:
    print("Neither .pdb or .cif file extension, please check")

#cmap = cm.get_cmap('RdYlBu')
cmap = cm.get_cmap('hot')

npy_path = fpath + "/interpret_output/" + npz
v_prot_all = np.load(npy_path + "_v_prot.npy")
x_prot_all = np.load(npy_path +"_xprot.npy")
interface_atoms = np.load(npy_path +"_edge_index.npy")[0,:]
v_prot = v_prot_all[interface_atoms]
x_prot = x_prot_all[interface_atoms]

prefixes = ["interface"]
diffs = np.load(npy_path +"_diffs.npy")
diffs = diffs/diffs.max()
#colors = ["black","yellow","deepsalmon","violetpurple"]
colors = ["black","deepsalmon","thorium","yellow"]
'''
for i in v_prot_all:
    if i not in v_prot:
        name = cmd.get_unused_name(prefixes[0])
        cmd.pseudoatom(name, vdw=1, pos=i.tolist(), color="grey") 
        cmd.select("sele", name)
        cmd.set("transparency",0.5,"sele")
'''
for i in range(len(v_prot)):
    item = v_prot[i]
    unused_name = cmd.get_unused_name(prefixes[0])
    color = cmap(diffs[i])
    cmd.set_color(unused_name + "_color", color[:3])
    color = colors[x_prot[i,:].tolist().index(1)]
    cmd.pseudoatom(prefixes[0], vdw=0.75 + 1.5*diffs[i], pos=item.tolist(), color=color)#unused_name + "_color")


            
cmd.show("spheres")
#cmd.set("sphere_transparency", 0.3) 
#cmd.set("specular", 0)
#cmd.bg_color("white")
cmd.set("depth_cue", "off")
cmd.set("ray_trace_mode", 3)
cmd.set("ray_trace_color" , "grey")
cmd.orient()
cmd.ray()

#cmd.png("{}_{}.png".format(pdb, key))

cmd.load(pdbpath + pdb)#.split(":")[0] + ".pdb")
cmd.zoom()

#print("here")
'''
try:
    if sys.argv[4]=="sticks":
        cmd.do("@vis_sticks.pml")   
except:
    cmd.do("@vis.pml")
'''
cmd.do("@"+fpath+"/vis_sticks.pml")
cmd.save(fpath + "/interpret_output/"+pdb+"_pymol_session.pse")
'''
cmd.do("hide surface")
cmd.do("hide sticks")
cmd.do("orient (vis)")
cmd.do("show surface")
cmd.do("select DNA, resn DA+DC+DG+DT+PSD")
cmd.do("hide surface, DNA")
#print(cmd.get_object_list())
cmd.do("ray")
cmd.do("png dummy.png")
#cmd.do("@"+pdb.split("_")[0]+".pml")
#cmd.do("dssr_block block_color=N red | minor 0.9 | major yellow")
'''
