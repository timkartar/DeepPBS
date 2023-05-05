# builtin modules
import os
import re
import subprocess
import logging 

# geobind modules
from .interpolator import Interpolator

def runAPBS(pqr, command="apbs", prefix="tmp", basedir='.', quiet=True, clean=True, space=0.3, cfac=1.7, fadd=20, gmemceil=400, processors=8, keep_dx=False):
    """ run APBS and return potential """
    
    # run psize to get grid length parameters
    stdout = subprocess.getoutput("psize --space={} --cfac={} --fadd={} --gmemceil={} '{}'".format(space, cfac, fadd, gmemceil, pqr))
    
    cglenMatch = re.search('Coarse grid dims = (\d*\.?\d+) x (\d*\.?\d+) x (\d*\.?\d+) A', stdout, re.MULTILINE)
    cgx = cglenMatch.group(1)
    cgy = cglenMatch.group(2)
    cgz = cglenMatch.group(3)
    fglenMatch = re.search('Fine grid dims = (\d*\.?\d+) x (\d*\.?\d+) x (\d*\.?\d+) A', stdout, re.MULTILINE)
    fgx = fglenMatch.group(1)
    fgy = fglenMatch.group(2)
    fgz = fglenMatch.group(3)
    dimeMatch = re.search('Num. fine grid pts. = (\d+) x (\d+) x (\d+)', stdout, re.MULTILINE)
    dx = dimeMatch.group(1)
    dy = dimeMatch.group(2)
    dz = dimeMatch.group(3)
    paraMatch = re.search("Parallel solve required", stdout, re.MULTILINE)
    if paraMatch:
        pdimeMatch = re.search("Proc. grid = (\d+) x (\d+) x (\d+)", stdout, re.MULTILINE)
        px = pdimeMatch.group(1)
        py = pdimeMatch.group(2)
        pz = pdimeMatch.group(3)
        elec_type = "mg-para"
        ofrac = "ofrac 0.1"
        pdime = "pdime {} {} {}".format(px, py, pz)
    else:
        elec_type = "mg-auto"
        ofrac = ""
        pdime = ""
    
    # run APBS
    pot = os.path.join(basedir, prefix+"_potential")
    acc = os.path.join(basedir, prefix+"_access")
    input_file = """READ
    mol pqr {}
END

ELEC
    {}
    mol 1
    
    dime {}
    cglen {}
    fglen {}
    cgcent mol 1
    fgcent mol 1
    
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.0
    srfm smol
    chgm spl2
    sdens 10.00
    srad 1.40
    swin 0.30
    temp 310.0
    {}
    {}
    
    ion charge +1 conc 0.15 radius 2.0
    ion charge -1 conc 0.15 radius 1.8
    
    write pot dx {}
    write smol dx {}
END""".format(
        pqr,
        elec_type,
        " ".join([dx, dy, dx]),
        " ".join([cgx, cgy, cgz]),
        " ".join([fgx, fgy, fgz]),
        ofrac,
        pdime,
        pot,
        acc
    )
    inFile = os.path.join(basedir, "{}.in".format(prefix))
    FH = open(inFile, "w")
    FH.write(input_file)
    FH.close()
    
    logging.info("Running APBS on input file: %s", inFile)
    outpt = subprocess.check_output([command, inFile], stderr=subprocess.STDOUT)
    
    I1 = Interpolator("{}.dx".format(pot))
    I2 = Interpolator("{}.dx".format(acc))
    # cleanup
    if clean:
        if os.path.exists(os.path.join(basedir, "{}.pdb".format(prefix))):
            os.remove(os.path.join(basedir, "{}.pdb".format(prefix)))
        os.remove(inFile)
        if os.access('io.mc', os.R_OK):
            os.remove('io.mc')
    
    if not keep_dx:
        # remove .dx files
        os.remove("{}.dx".format(pot))
        os.remove("{}.dx".format(acc))
    
    return I1, I2 

