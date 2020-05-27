from builtins import next
from builtins import object
import numpy as np
import mdtraj as mdt
import mdtraj.formats as mdt_load
import mdtraj.formats
from mdtraj import FormatRegistry as mdt_FormatRegistry
from mdtraj.core.topology import Topology as mdt_Topology
from mdtraj.core.trajectory import Trajectory as mdt_Trajectory
from mdtraj.core.trajectory import _TOPOLOGY_EXTS as mdt_TOPOLOGY_EXTS
import ase
from ase import io as ase_io
import ase.io.formats
import logging
import os
import sys
import re
import math
import inspect
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    import MDAnalysis as mda
    assert str(w[-1].message)
import MDAnalysis.core.universe as mda_u
import MDAnalysis.coordinates as mda_c
import panedr
import parmed as pmd
import pymolfile as pym
from nomadcore.md_data_access import GromosTopoObjects as gto
from nomadcore.md_data_access.GromosTopoObjects import get_zipType, get_fileExtensions
import io
import struct
import copy

mda_coordinates_modules = tuple(x[1] for x in inspect.getmembers(mda_c,inspect.ismodule))
logger = logging.getLogger(__name__)

# The info in the table below is taken from MDAnalysis and definitions 
# are checked with the info in Amber, CHARMM, Gromacs, and GROMOS
# The units in the table are in atomic masses (u)
# Please see:
# https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/topology/tables.py
ELEMENTS_MASS_TABLE = {
        "Ac" : 227.028,   "Al" : 26.981539, "Am" : 243,       "Sb" : 121.757,
        "Ar" : 39.948,    "As" : 74.92159,  "At" : 210,       "Ba" : 137.327,
        "Bk" : 247,       "Be" : 9.012182,  "Bi" : 208.98037, "Bh" : 262,
        "B"  : 10.811,    "Br" : 79.90400,  "Cd" : 112.411,   "Ca" : 40.08000,
        "Cf" : 251,       "C"  : 12.01100,  "Ce" : 140.11600, "Cs" : 132.90000,
        "Cl" : 35.45000,  "Cr" : 51.9961,   "Co" : 58.9332,   "Cu" : 63.54600,
        "Cm" : 247,       "Db" : 262,       "Dy" : 162.5,     "Es" : 252,
        "Er" : 167.26,    "Eu" : 151.965,   "Fm" : 257,       "F"  : 18.99800,
        "Fr" : 223,       "Gd" : 157.25,    "Ga" : 69.723,    "Ge" : 72.61,
        "Au" : 196.96654, "Hf" : 178.49,    "Hs" : 265,       "He" : 4.00260,
        "Ho" : 164.93032, "H"  : 1.00800,   "In" : 114.82,    "I"  : 126.90450,
        "Ir" : 192.22,    "Fe" : 55.84700,  "Kr" : 83.8,      "La" : 138.9055,
        "Lr" : 262,       "Pb" : 207.2,     "Li" : 6.941,     "Lu" : 174.967,
        "Mg" : 24.30500,  "Mn" : 54.93805,  "Mt" : 266,       "Md" : 258,
        "Hg" : 200.59,    "Mo" : 95.94,     "N"  : 14.00700,  "Na" : 22.98977,
        "Nd" : 144.24,    "Ne" : 20.17970,  "Np" : 237.048,   "Ni" : 58.6934,
        "Nb" : 92.90638,  "No" : 259,       "Os" : 190.2,     "O"  : 15.99900,
        "Pd" : 106.42,    "P"  : 30.97400,  "Pt" : 195.08,    "Pu" : 244,
        "Po" : 209,       "K"  : 39.10200,  "Pr" : 140.90765, "Pm" : 145,
        "Pa" : 231.0359,  "Ra" : 226.025,   "Rn" : 222,       "Re" : 186.207,
        "Rh" : 102.9055,  "Rb" : 85.46780,  "Ru" : 101.07,    "Rf" : 261,
        "Sm" : 150.36,    "Sc" : 44.95591,  "Sg" : 263,       "Se" : 78.96,
        "Si" : 28.0855,   "Ag" : 107.8682,  "Na" : 22.989768, "Sr" : 87.62,
        "S"  : 32.06000,  "Ta" : 180.9479,  "Tc" : 98,        "Te" : 127.6,
        "Tb" : 158.92534, "Tl" : 204.3833,  "Th" : 232.0381,  "Tm" : 168.93421,
        "Sn" : 118.71,    "Ti" : 47.88,     "W"  : 183.85,    "U"  : 238.0289,
        "V"  : 50.9415,   "Xe" : 131.29,    "Yb" : 173.04,    "Y"  : 88.90585,
        "Zn" : 65.37000,  "Zr" : 91.224
        }  

TEXTCHARS = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})

def is_file_binary(fName, checkBytes=None):
    testin=None
    if checkBytes is None:
        checkBytes = 1024 
    try:
        with open(fName, 'rb') as fin:
            testin = fin.read(checkBytes)
    except(FileNotFoundError, IOError):
        pass
    if testin:
        if is_binary_string(testin):
            return True
        else:
            return False
    else:
        return None

def is_binary_string(inbytes):
    return bool(inbytes.translate(None, TEXTCHARS))

def bytes_per_record(record_type):

    if record_type == 'f': return 4
    if record_type == 'd': return 8
    if record_type == 'i': return 4
    if record_type == 's': return 1
    if record_type == 's4': return 4

    # Default we read 1 byte
    return 1

def read_start_marker(f):
    # Read Start Marker
    try:
        struct.unpack('i', f.read(4))
    except(struct.error):
        pass

def read_end_marker(f):
    # Read End Marker
    try:
        struct.unpack('i', f.read(4))
    except(struct.error):
        pass

def read_record(f, records, marker=False):

    if marker:
        # Read Start Marker
        read_start_marker(f)

    res=[]
    for record_type, length, smarker, emarker in records:
        if smarker:
            # Read Start Marker
            read_start_marker(f)
        x=[]
        if length<0:
            length= res[-length-1][0]
        if 's' in record_type:
            record_type = 's'
        for i in range(length):
            try:
                x.append(struct.unpack(record_type, f.read(bytes_per_record(record_type)))[0])
            except(struct.error):
                return None
        if 's' in record_type:
            x = ''.join([ s.decode('utf-8') for s in x])
        res.append(x)
        if emarker:
            # Read End Marker
            read_end_marker(f)

    if marker:
        # Read End Marker
        read_end_marker(f)

    return res

def get_any_element_name(atomname):
    elementlist = ELEMENTS_MASS_TABLE.keys()
    # check if the element is in list
    # but name is upper- or lower-case
    for name in elementlist:
        if atomname.lower() in name.lower():
            return name
    # check if the first two letters define
    # element name
    if len(atomname)>1:
        for name in elementlist:
            if atomname[0:2].lower() in name.lower():
                return name
    if len(atomname)>0:
        # check if the first letter defines
        # element name
        for name in elementlist:
            try:
                lowername = atomname[0:1].lower
                if lowername in name.lower:
                    return name
            except TypeError:
                pass
    return None

def get_element_name(atomname):
    elementlist = ELEMENTS_MASS_TABLE.keys()
    # check if the element is in list
    # but name is upper- or lower-case
    for name in elementlist:
        if atomname == name:
            return name
    # check if the first two letters define
    # element name
    if len(atomname)>1:
        for name in elementlist:
            if atomname[0:2] == name:
                return name
    if len(atomname)>0:
        # check if the first letter defines
        # element name
        for name in elementlist:
            if atomname[0:1] == name:
                return name
    if len(atomname)<3:
        return atomname[0:1]
    else:
        return atomname.upper()
    return None

def charmm_stream_reader(slines, topo=None):
    if topo is None:
        topo = pmd.charmm.CharmmParameterSet()

    lines = []
    title = []
    command = []

    return topo

def charmm_coor_reader(fname, ftextmem=None):
    coorDict = None
    if ftextmem is not None:
        title = []
        natoms = 0
        countatoms = 0
        atomlist = []
        for line in ftextmem:
            # If the line is a title line
            # add it to title list
            if(line.startswith('*') or 'TITLE>' in line.upper()):
                title.append(line)
            # else start readinf COOR format
            else:
                if natoms<1:
                    natoms = int(line.split()[0])
                else:
                    if countatoms<natoms:
                        atomlist.append(line.split())
                    countatoms+=1
        if coorDict is None:
            coorDict = {}
        coorDict.update({'binary':False})
        coorDict.update({'title':title})
        coorDict.update({'numatoms':natoms})
        coorDict.update({'atomlist':atomlist})
        coords=np.asarray(atomlist)
        if len(coords[0,:])>6:
            coorDict.update({'coords':np.asarray([[[float(x) for x in v] for v in coords[:,4:7]]])})
        else:
            coorDict.update({'coords':np.asarray([[[float(x) for x in v] for v in coords[:,3:6]]])})
    else:
        if is_file_binary(fname):
            try:
                with open(fname, 'rb') as fin:
                    hdr, icntrl = read_record(f=fin, records=[
                        ('s', 4, True, False), ('i', 20, False, True)])
                    if hdr == 'COOR':
                        ntitl = read_record(f=fin, records=[('i', 1, True, False)])
                        ntitl=ntitl[0][0]
                        for t in range(int(ntitl)):
                            title = read_record(f=fin, records=[('s', 80, False, False)])
                        endtitle = read_record(f=fin, records=[('i', 1, True, False)])
                        natoms = read_record(f=fin, records=[('i', 1, False, False)])
                        natoms = natoms[0][0]
                        xt = read_record(f=fin, records=[('i', 1, False, False)])
                        wx = read_record(f=fin, records=[('f', natoms, True, True)])
                        wy = read_record(f=fin, records=[('f', natoms, True, True)])
                        wz = read_record(f=fin, records=[('f', natoms, True, True)])
                        ww = read_record(f=fin, records=[('f', natoms, True, True)])
                        if ww is None:
                            ww = [0.0 for atm in range(natoms)]
                        if coorDict is None:
                            coorDict = {}
                        coords=[]
                        for r in range(len(wx)):
                            coords.append(np.column_stack([wx[r],wy[r],wz[r]]))
                        coorDict.update({'binary':True})
                        coorDict.update({'title':title})
                        coorDict.update({'numatoms':natoms})
                        coorDict.update({'xt':xt})
                        coorDict.update({'coords':np.asarray(coords)})
                        coorDict.update({'ww':ww})
            except(FileNotFoundError, IOError):
                pass
        else:
            try:
                with open(fname, 'r') as fin:
                    title = []
                    natoms = 0
                    countatoms = 0
                    atomlist = []
                    for line in fin:
                        # If the line is a title line
                        # add it to title list
                        if(line.startswith('*') or 
                           line.startswith('TITLE')):
                            title.append(line)
                        # else start readinf COOR format
                        else:
                            if natoms<1:
                                try:
                                    natoms = int(line.split()[0])
                                except(ValueError):
                                    return coorDict
                            else:
                                if countatoms<natoms:
                                    atomlist.append(line.split())
                                countatoms+=1
                    if coorDict is None:
                        coorDict = {}
                    coorDict.update({'binary':False})
                    coorDict.update({'title':title})
                    coorDict.update({'numatoms':natoms})
                    coorDict.update({'atomlist':atomlist})
                    coords=np.asarray(atomlist)
                    if len(coords[0,:])>6:
                        coorDict.update({'coords':np.asarray([[[float(x) for x in v] for v in coords[:,4:7]]])})
                    else:
                        coorDict.update({'coords':np.asarray([[[float(x) for x in v] for v in coords[:,3:6]]])})
            except(FileNotFoundError, IOError):
                pass
    return coorDict

def get_dir_base_extension(file_name):
    """ Splits directory, file base and file extensions

        Returns: directory without leading '/', 
                 file base name, and file extension without '.'
    """
    file_base, file_extension_with_dot = os.path.splitext(os.path.basename(file_name))
    file_extension = file_extension_with_dot.split(".")[-1]
    file_dir = os.path.dirname(file_name)
    return file_dir, file_base, file_extension
            
def pmdConvertTopoDict(topoStor, topoPMD):
    """ Function to convert ParmEd topology info
        from CharmmPsfFile and CharmmParameterSet
        to MDDataAccess

    """
    topo=None
    if(isinstance(topoPMD, pmd.charmm.CharmmParameterSet) or
       isinstance(topoPMD, pmd.charmm.CharmmPsfFile) or
       isinstance(topoPMD, pmd.tinker.tinkerfiles.DynFile) or
       isinstance(topoPMD, pmd.tinker.tinkerfiles.XyzFile) or
       (topoPMD is None and topoStor.charmmcoor is not None) or 
       topoStor.charmmcoortopo is not None):
        topo = topoPMD

        def getatomall(atom):
            atmid=0
            atmname = atom[0]
            atmtyp = atom[1]
            if atom[10]:
                atmmass = float(atom[10])
            else:
                atmmass = None
            if atom[3]:
                atmresid = atom[3]
                atmres = atom[2]
            else:
                atmresid = None
                atmres = None
            if atom[4]:
                atmsegid = atom[4]
            else:
                atmsegid = None
            if atom[11]:
                atmchrg = float(atom[11])
            else:
                atmchrg = None
            if atom[12]:
                atmrad = float(atom[12])
            else:
                atmrad = None
            if atom[9]:
                atmbfac = float(atom[9])
            else:
                atmbfac = None
            atm_unique = ''
            if atmname:
                atm_unique = atm_unique + atmname 
            if atmtyp:
                atm_unique = atm_unique + "-" + atmtyp
            if atom[16]:
                aeps = float(atom[16])
            else:
                aeps = None
            if atom[17]:
                arm = float(atom[17])
            else:
                arm = None
            if atom[18]:
                aeps14 = float(atom[18])
            else:
                aeps14 = None
            if atom[19]:
                arm14 = float(atom[19])
            else:
                arm14 = None
            if atom[20]:
                anb = [[x,float(y)] for x,y in atom[20]]
            else:
                anb = None
            #if atmres:
            #    atm_unique = atm_unique + "-" + atmres
            #if atmsegid:
            #    atm_unique = atm_unique + "-" + atmsegid
            return [atmid,atm_unique,atmname,
                    atmtyp,atmres,atmresid,atmsegid,
                    atmmass,atmchrg,atmrad,atmbfac,
                    aeps,arm,aeps14,arm14,anb,None]

        def checkatomsdiff(atom1,atom2,checklist):
            index=0
            target=len(checklist)
            for t in checklist:
                if atom1[t] is not atom2[t]:
                    index += 1
            if index==target:
                return True
            else:
                return False

        def atom_respairs(atom1,atom2):
            return [get_atomres(atom1),get_atomres(atom2)]

        def atompairs(atom1,atom2):
            return [getatom(atom1),getatom(atom2)]

        def atompairids(atom1,atom2):
            return [atom1[0],atom2[0]]
     
        def getX(x, attr, attr2=None):
            try:
                if attr2 is not None:
                    xa = getattr(x, attr)
                    return getX(xa, attr2)
                else:
                    return getattr(x, attr)
            except AttributeError:
                return None

        structure=None
        segmentList = None
        residueList = None
        atom_element_list = []
        nbepsList = None
        nbrList = None
        nbeps14List = None
        nbr14List = None
        nbfixList = None
        # If topo is a PsfFile, the parameter files were not supplied.
        # Hence, we can only extract topo info from Psf class in topo.
        if(isinstance(topoPMD, pmd.charmm.CharmmPsfFile) or
           isinstance(topoPMD, pmd.tinker.tinkerfiles.XyzFile)):
            unassignedAtomType = False
            try:
                structure = [[getX(x,'name'), getX(x,'type'), getX(x,'residue','name'), 
                              getX(x,'residue','_idx'), getX(x,'residue','segid'), 
                              getX(x,'residue','chain'), getX(x,'altloc'), 
                              getX(x,'irotat'), getX(x,'occupancy'), getX(x,'bfactor'), 
                              getX(x,'mass'), getX(x,'_charge'), getX(x,'solvent_radius'), 
                              getX(x,'atomic_number'), getX(x,'tree'), getX(x,'atom_type','name'), 
                              getX(x,'atom_type','epsilon'), getX(x,'atom_type','rmin'), 
                              getX(x,'atom_type','epsilon_14'), getX(x,'atom_type','rmin_14'), 
                              getX(x,'atom_type','nbfix')] for x in topoPMD.atoms]
            except AttributeError:
                unassignedAtomType = True
                structure = [[getX(x,'name'), getX(x,'type'), getX(x,'residue.name'), 
                              getX(x,'residue._idx'), getX(x,'residue.segid'), 
                              getX(x,'residue.chain'), getX(x,'altloc'), 
                              getX(x,'irotat'), getX(x,'occupancy'), getX(x,'bfactor'), getX(x,'mass'), 
                              getX(x,'._charge'), getX(x,'solvent_radius'), getX(x,'atomic_number'), 
                              getX(x,'tree'), '', None, None, None, None, 
                              None] for x in topoPMD.atoms]

            chainList = [a[5] for a in structure]
            segmentList = [a[4] for a in structure]
            residList = [a[3] for a in structure]
            residueList = [a[2] for a in structure]
            nbepsList = [a[16] for a in structure]
            nbrList = [a[17] for a in structure]
            nbeps14List = [a[18] for a in structure]
            nbr14List = [a[19] for a in structure]
            nbfixList = [a[20] for a in structure]
            atomList = structure
            atomAllList = [getatomall(a) for a in structure]
            atom_name_list = [a[0] for a in structure]
            #for ai, atom in enumerate(structure):
            #    if atom[2] in ELEMENTS_MASS_TABLE.keys():
            #        element = atom[2]
            #    else:
            #        element = get_element_name(atom[2])
            #        if element is None:
            #            element = atom[1]
            #    atom_element_list.append(element) 
        # If topo is ParSet, the parameter file or files might be supplied.
        # If Psf file is also supplied with the parameters, we can extract 
        # full information of the topology with parameter values.
        # However, if there is no Psf file and user only supplied COOR CARD
        # or CRD file, than we can only extract atom definitions through 
        # these two files and we try to fill in the blanks of the topology
        # puzzle in CHARMM with parameter set class in ParmEd.
        elif(isinstance(topoPMD, pmd.charmm.CharmmParameterSet) and
            (topoStor.charmmcoor is not None or 
             topoStor.charmmcoortopo is not None)):
            coortopo=False
            #ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
            if topoStor.charmmcoortopo is not None:
                if len(topoStor.charmmcoortopo[0])>8:
                    structure = [[at[3], at[2], at[2],
                                  at[1], at[7], at[8],
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None] for at in topoStor.charmmcoortopo]
                else:
                    structure = [[at[3], at[2], at[2],
                                  at[1], None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None] for at in topoStor.charmmcoortopo]
                coortopo=True
            elif topoStor.charmmcoor is not None:
                if 'atomlist' in topoStor.charmmcoor:
                    topostruct = topoStor.charmmcoor['atomlist']
                    if len(topostruct[0])>8:
                        structure = [[at[3], at[2], at[2],
                                      at[1], at[7], at[8],
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None] for at in topostruct]
                    else:
                        structure = [[at[3], at[2], at[2],
                                      at[1], None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None] for at in topostruct]
                    coortopo=True
            if coortopo is True:
                param_types = [[topoPMD.atom_types[at].name, topoPMD.atom_types[at].epsilon,
                                topoPMD.atom_types[at].rmin, topoPMD.atom_types[at].epsilon_14,
                                topoPMD.atom_types[at].rmin_14, topoPMD.atom_types[at].nbfix] for at in topoPMD.atom_types]
                nbond_params = {}
                for par in param_types:
                    nbond_params.update({par[0]:par[1:6]})
                structstor = structure[:]
                structure = []
                for ai, at in enumerate(structstor):
                    structure.append(at)
                    structure[ai].extend([None])
                    if at[2] in nbond_params:
                        structure[ai].extend(nbond_params[at[2]][0:4])
                        if nbond_params[at[2]][4]:
                            structure[ai].extend(nbond_params[at[2]][4])
                        else:
                            structure[ai].extend([None])
                chainList = [a[5] for a in structure]
                segmentList = [a[4] for a in structure]
                residList = [a[3] for a in structure]
                residueList = [a[2] for a in structure]
                nbepsList = [a[16] for a in structure]
                nbrList = [a[17] for a in structure]
                nbeps14List = [a[18] for a in structure]
                nbr14List = [a[19] for a in structure]
                nbfixList = [a[20] for a in structure]
                atomList = structure
                atomAllList = [getatomall(a) for a in structure]
                atom_name_list = [a[0] for a in structure]
                #for ai, atom in enumerate(structure):
                #    if atom[2] in ELEMENTS_MASS_TABLE.keys():
                #        element = atom[2]
                #    else:
                #        element = get_element_name(atom[2])
                #        if element is None:
                #            element = atom[1]
                #    atom_element_list.append(element) 

        # There is also one possibility that we do not have 
        # Psf or parameter set but only COOR info. If this is 
        # the case and the info is not from a binary CRD file,        
        # we can only extract names of residues, atom kinds 
        # and types from the coordinate card info and 
        elif(topoStor.charmmcoor is not None or 
             topoStor.charmmcoortopo is not None):
            coortopo=False
            #ATOMNO RESNO   RES  TYPE  X     Y     Z   SEGID RESID Weighting
            if topoStor.charmmcoortopo is not None:
                if len(topoStor.charmmcoortopo[0])>8:
                    structure = [[at[3], at[2], at[2],
                                  at[1], at[7], at[8],
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None] for at in topoStor.charmmcoortopo]
                else:
                    structure = [[at[3], at[2], at[2],
                                  at[1], None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None,
                                  None, None, None] for at in topoStor.charmmcoortopo]
                coortopo=True
            elif topoStor.charmmcoor is not None:
                if 'atomlist' in topoStor.charmmcoor:
                    topostruct = topoStor.charmmcoor['atomlist']
                    if len(topostruct[0])>8:
                        structure = [[at[3], at[2], at[2],
                                      at[1], at[7], at[8],
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None] for at in topostruct]
                    else:
                        structure = [[at[3], at[2], at[2],
                                      at[1], None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None,
                                      None, None, None] for at in topostruct]
                    coortopo=True
            if coortopo is True:
                chainList = [a[5] for a in structure]
                segmentList = [a[4] for a in structure]
                residList = [a[3] for a in structure]
                residueList = [a[2] for a in structure]
                atomList = topo.structure
                atomAllList = [getatomall(a) for a in structure]
                atom_name_list = [a[0] for a in structure]
                #for ai, atom in enumerate(structure):
                #    if atom[2] in ELEMENTS_MASS_TABLE.keys():
                #        element = atom[2]
                #    else:
                #        element = get_element_name(atom[2])
                #        if element is None:
                #            element = atom[1]
                #    atom_element_list.append(element) 

        if structure is None:
            return topoStor

        count=0
        for i in range(len(atomAllList)):
            atomAllList[i][0]=count
            count+=1

        types_list = list(set([a[1] for a in atomAllList]))
        #atom_type_list = [a[3] for a in atomList]

        atom_all_list = np.asarray(atomAllList)
        atom_names = atom_all_list[:,2]
        atom_types = atom_all_list[:,3]
        atom_masses = atom_all_list[:,7]
        atom_charges = atom_all_list[:,8]
        atom_radiuses = atom_all_list[:,9]
        atom_bfactors = atom_all_list[:,10]
        #print("Atom Bonds:",topo.bonds)
        #print("Atom Angles:",topo.angles)
        #print("Atom Dihedrals:",topo.dihedrals)
        #print("Atom Impropers:",topo.impropers)


        atom_type_dict = {}
        atomtypesDict = {}
        atomnameDict = {}
        massDict = {}
        elementDict = {}
        radiusDict = {}
        chargeDict = {}
        bfactorDict = {}
        epsilonDict = {}
        rminDict = {}
        epsilon14Dict = {}
        rmin14Dict = {}
        nbfixDict = {}
        atomlabelList = []
        atomlabelDict = []
        for ielm in range(len(types_list)):
            elm = types_list[ielm]
            atom_type_dict.update({elm : ielm+1})
            typelabelList = []
            for ai, atom in enumerate(atomAllList):
                if atom[2] in ELEMENTS_MASS_TABLE.keys():
                    element = atom[2]
                else:
                    element = get_element_name(atom[2])
                    if element is None:
                        element = atom[1]
                if elm == atom[1]:
                    elementDict.update({elm : element})
                    atomnameDict.update({elm : atom[2]})
                    atomtypesDict.update({elm : atom[3]})
                    massDict.update({elm : atom[7]})
                    chargeDict.update({elm : atom[8]})
                    radiusDict.update({elm : atom[9]})
                    bfactorDict.update({elm : atom[10]})
                    epsilonDict.update({elm : atom[11]})
                    rminDict.update({elm : atom[12]})
                    epsilon14Dict.update({elm : atom[13]})
                    rmin14Dict.update({elm : atom[14]})
                    nbfixDict.update({elm : atom[15]})
                    typelabelList.append(atom[0])
                atomAllList[ai][16]=element
            atomlabelDict.append(typelabelList)

        for ai, atom in enumerate(atomAllList):
            atom_element_list.append(atom[16]) 

        for atom in atomAllList:
            atomlabelList.append([atom[0],atom_type_dict[atom[1]]])
        
        system_name = ''
        if atom_element_list:
            system_name = system_name + ''.join([
                el + str(atom_element_list.count(el)) for el in set(
                    atom_element_list)])
        if segmentList is not None:
            segsysname = '-'.join([str(seg) for seg in list(
                    set(segmentList)) if seg is not None])
            if segsysname:
                system_name = system_name + '-' + segsysname
        if residueList is not None:
            ressysname = '-'.join([str(res) for res in list(
                    set(residueList)) if res is not None])
            if ressysname:
                system_name = system_name + '-' + ressysname

        massList = list(massDict.values())
        atom_type_list = list(atomtypesDict.values())
        name_list = list(atomnameDict.values())
        elementList = list(elementDict.values())
        radiusList = list(radiusDict.values())
        chargesList = list(chargeDict.values())
        bfactorList = list(bfactorDict.values())

        if(isinstance(topoPMD, pmd.charmm.CharmmPsfFile) or
           isinstance(topoPMD, pmd.tinker.tinkerfiles.XyzFile)):
            if getattr(topoPMD, 'bonds'):
                if topoPMD.bonds[0].type is not None:
                    for bi, bondtype in enumerate(topoPMD.bonds[0].type.list):
                        topoPMD.bonds[0].type.list[bi]._idx=bi
                    topo_bond_list = np.asarray([
                        [x.atom1._idx, x.atom2._idx, x.type._idx, 
                         x.type.k, x.type.req] for x in topoPMD.bonds
                        ])
                else:
                    topo_bond_list = np.asarray([
                        [x.atom1._idx, x.atom2._idx] for x in topoPMD.bonds
                        ])
                topbList = topo_bond_list[:,0:2]
            else:
                topbList = []
        elif(isinstance(topoPMD, pmd.charmm.CharmmParameterSet) and
            (topoStor.charmmcoor is not None or 
             topoStor.charmmcoortopo is not None)):
            if topoPMD.bond_types is not None:
                bi=0
                search_types = list(set(atom_type_list))
                search_names = search_types[:]
                search_names.extend(list(set(atom_name_list)))
                search_names.extend(list(set(elementList)))
                search_names = list(set(search_names))
                search_bonds = []
                for bondtype in topoPMD.bond_types:
                    for typ in search_names:
                        if typ in bondtype[0] and bondtype[1] in search_types:
                            topoPMD.bond_types[bondtype]._idx=bi
                            search_bonds.append(bondtype)
                            bi+=1
                        elif typ in bondtype[1] and bondtype[0] in search_types:
                            topoPMD.bond_types[bondtype]._idx=bi
                            search_bonds.append(bondtype)
                            bi+=1
                bond_found_list = []
                bond_found_names = []
                # Full bond search, the worst case here.
                for btyp in search_bonds:
                    for ati, at1 in enumerate(atomAllList[0:len(atomAllList)]):
                        for at2 in atomAllList[ati:len(atomAllList)]:
                            if at1[5] != at2[5]: # We assume bonds are in same residue since
                                continue         # COOR file only includes residue ids for each atom
                            if at1[0] == at2[0]:
                                continue
                            if((at1[2] in btyp[0] or at1[3] in btyp[0] or at1[16] in btyp[0]) and 
                               (at2[2] in btyp[1] or at2[3] in btyp[1] or at2[16] in btyp[1])):
                                if(str(at1[0])+'-'+str(at2[0]) in bond_found_names or
                                   str(at2[0])+'-'+str(at1[0]) in bond_found_names):
                                       continue
                                bond_found_list.append([at1[0], at2[0], topoPMD.bond_types[btyp]._idx, 
                                    topoPMD.bond_types[btyp].k,topoPMD.bond_types[btyp].req])
                                bond_found_names.append(str(at1[0])+'-'+str(at2[0]))
                                bond_found_names.append(str(at2[0])+'-'+str(at1[0]))
                            if((at1[2] in btyp[1] or at1[3] in btyp[1] or at1[16] in btyp[1]) and 
                               (at2[2] in btyp[0] or at2[3] in btyp[0] or at2[16] in btyp[0])):
                                if(str(at1[0])+'-'+str(at2[0]) in bond_found_names or
                                   str(at2[0])+'-'+str(at1[0]) in bond_found_names):
                                       continue
                                bond_found_list.append([at1[0], at2[0], topoPMD.bond_types[btyp]._idx, 
                                    topoPMD.bond_types[btyp].k,topoPMD.bond_types[btyp].req])
                                bond_found_names.append(str(at1[0])+'-'+str(at2[0]))
                                bond_found_names.append(str(at2[0])+'-'+str(at1[0]))

                topo_bond_list = np.asarray(bond_found_list)
                if len(topo_bond_list)>0:
                    topbList = topo_bond_list[:,0:2]
                else:
                    topbList = []
            else:
                topbList = []
            #if topoPMD.angle_types is not None:
            #    for angletype in topoPMD.angle_types:
            #        for typ in set(atom_type_list):
            #            if typ in angletype:
            #                print("PMDTOPO: angles:",angletype, topoPMD.angle_types[angletype].__dict__)
            #if topoPMD.dihedral_types is not None:
            #    for dihtype in topoPMD.dihedral_types:
            #        for typ in set(atom_type_list):
            #            if typ in dihtype:
            #                print("PMDTOPO: dihedrals:",dihtype, topoPMD.dihedral_types[dihtype].__dict__)
            #if topoPMD.improper_types is not None:
            #    for imptype in topoPMD.improper_types:
            #        for typ in set(atom_type_list):
            #            if typ in imptype:
            #                print("PMDTOPO: impropers:",imptype, topoPMD.improper_types[imptype].__dict__)
        else:
            topbList = []

        topbNames = []
        interDict = {}
        interTypeDict = {}
        topbCont=False
        if len(topbList)>0:
            topbCont=True

        if topbCont:
            for pair in topbList:
                topbNames.append(atomAllList[int(pair[0])][1] + '-' + atomAllList[int(pair[1])][1])
            
            topb=list(set(topbNames))
            interNum = 0
            for tb in topb:
                bondList = []
                typeList = []
                noninter = True
                bc = 0 
                for bb in topbList:
                    b=[int(x) for x in bb]
                    topt=atomAllList[b[0]][1] + '-' + atomAllList[b[1]][1]
                    if topt == tb:
                        reslist=[atomAllList[b[0]][4],atomAllList[b[1]][4]]
                        atmlist=[atomAllList[b[0]],atomAllList[b[1]]]
                        atmidlist=[b[0],b[1]]
                        bondList.append(atmidlist)
                        interDict.update({interNum : bondList})
                        if noninter:
                            noninter = False
                            typeList.extend(list([
                                atom_type_dict[atmlist[0][1]],
                                atom_type_dict[atmlist[1][1]]
                                ]))
                            interTypeDict.update({interNum : typeList})
                interNum += 1

#        for ielm in range(len(atom_type_list)-1):
#            for jelm in range(ielm+1, len(atom_type_list)):
#                aelm = atom_type_list[ielm]
#                belm = atom_type_list[jelm]
#                bondList = []
#                typeList = []
#                bondid = 0
#                noninter = True
#                for key in topdk:
#                    topt = topd[key]
#                    for b in topt:
#                        reslist=atom_respairs(b.atoms[0],b.atoms[1])
#                        atmlist=atompairs(b.atoms[0],b.atoms[1])
#                        atmidlist=atompairids(b.atoms[0],b.atoms[1])
#                        if((aelm == str(atmlist[0][2]) and belm == str(atmlist[1][2])) or 
#                           (aelm == str(atmlist[1][2]) and belm == str(atmlist[0][2]))):
#                            bondList.append(atmidlist)
#                            interDict.update({interNum : bondList})
#                            if noninter:
#                                noninter = False
#                                typeList.extend(list([
#                                    atom_type_dict[aelm],
#                                    atom_type_dict[belm]
#                                ]))
#                                interTypeDict.update({interNum : typeList})
#                                interNum += 1
#                        bondid += 1

        #atomIndex = np.arange(len(residueList))
        #atom_to_residue = np.zeros((len(residueList), 2), dtype=int)
        #atom_to_residue[:,0] = atomIndex+1
        #atom_to_residue[:,1] = np.array(residueList)+1
       

        topoStor.topoDict.update({"system_name" : system_name})
        topoStor.topoDict.update({"atom_name_list" : atom_name_list})
        topoStor.topoDict.update({"atom_label_list" : atomlabelList})
        topoStor.topoDict.update({"name_list" : name_list})
        topoStor.topoDict.update({"types_list" : types_list})
        topoStor.topoDict.update({"atom_type_list" : atom_type_list})
        topoStor.topoDict.update({"atom_mass_list" : massList})
        topoStor.topoDict.update({"atom_element_list" : atom_element_list})
        topoStor.topoDict.update({"element_list" : elementList})
        topoStor.topoDict.update({"atom_radius_list" : radiusList})
        topoStor.topoDict.update({"atom_charge_list" : chargesList})
        topoStor.topoDict.update({"interactions_dict" : interDict})
        topoStor.topoDict.update({"interactions_type_dict" : interTypeDict})
        topoStor.topoDict.update({"mol_list" : residueList})
        #topoStor.topoDict.update({"atom_to_mol_list" : atom_to_residue})
        topoStor.topoDict.update({"residue_list" : residueList})
    return topoStor

def gtoConvertTopoDict(topoStor, topoGTO):
    """ Function to convert GromosTopology topology info
        to MDDataAccess

    """
    topo=None
    if isinstance(topoGTO, gto.GromosTopology):
        topo = topoGTO

        def getatomall(atom):
            atmid=atom._id
            atmname = atom.name
            atmtyp = atom.atom_type.name
            atmsegid = None
            if atom.mass:
                atmmass = float(atom.mass)
            else:
                atmmass = None
            if atom.res_id:
                atmresid = int(atom.res_id)
                atmres = atom.residue
            else:
                atmresid = None
                atmres = None
            if atom.charge:
                atmchrg = float(atom.charge)
            else:
                atmchrg = None
            if atom.radius:
                atmrad = float(atom.radius)
            else:
                atmrad = None
            if atom.bfactor:
                atmbfac = float(atom.bfactor)
            else:
                atmbfac = None
            atm_unique = atom.atom_type.unique_name
            return [atmid,atm_unique,atmname,atmtyp,atmres,atmresid,atmsegid,atmmass,atmchrg,atmrad,atmbfac]

        def checkatomsdiff(atom1,atom2,checklist):
            index=0
            target=len(checklist)
            for t in checklist:
                if atom1[t] is not atom2[t]:
                    index += 1
            if index==target:
                return True
            else:
                return False

        def atom_respairs(atom1,atom2):
            return [get_atomres(atom1),get_atomres(atom2)]

        def atompairs(atom1,atom2):
            return [getatom(atom1),getatom(atom2)]

        def atompairids(atom1,atom2):
            return [atom1[0],atom2[0]]

        residList = [i for i, a in enumerate(topo._residues)]
        residueList = topo.residues
        atomList = topo.atoms
        atomAllList = [getatomall(a) for a in topo.atoms]

        #types_list = [a.name for a in topo.atomtypes]
        types_list = list(set([a.unique_name for a in topo.atomtypes]))
        atom_name_list = [a.name for a in topo.atoms]

        atom_element_list = []
        for atom in topo.atoms:
            if atom.name in ELEMENTS_MASS_TABLE.keys():
                element = atom.name
            else:
                element = get_element_name(atom.name)
                if element is None:
                    element = atom.name
            atom_element_list.append(element) 

        system_name = ''
        if atom_element_list:
            system_name = system_name + ''.join([
                el + str(atom_element_list.count(el)) for el in set(
                    atom_element_list)])
        if residueList is not None:
            system_name = system_name + '_' + '-'.join([
                str(res) for res in list(
                    set(residueList)) if res is not None])

        atom_all_list = np.asarray(atomAllList)
        atom_names = atom_all_list[:,2]
        atom_types = atom_all_list[:,3]
        atom_masses = atom_all_list[:,7]
        atom_charges = atom_all_list[:,8]
        atom_radiuses = atom_all_list[:,9]
        atom_bfactors = atom_all_list[:,10]

        atom_type_dict = {}
        atomtypesDict = {}
        atomnameDict = {}
        massDict = {}
        elementDict = {}
        radiusDict = {}
        chargeDict = {}
        bfactorDict = {}
        atomlabelList = []
        atomlabelDict = []
        for ielm in range(len(types_list)):
            elm = types_list[ielm]
            atom_type_dict.update({elm : ielm+1})
            typelabelList = []
            for atom in atomAllList:
                if elm == atom[1]:
                    if atom[2] in ELEMENTS_MASS_TABLE.keys():
                        element = atom[2]
                    else:
                        element = get_element_name(atom[2])
                        if element is None:
                            element = atom[2]
                    elementDict.update({elm : element})
                    atomnameDict.update({elm : atom[2]})
                    atomtypesDict.update({elm : atom[3]})
                    massDict.update({elm : atom[7]})
                    chargeDict.update({elm : atom[8]})
                    radiusDict.update({elm : atom[9]})
                    bfactorDict.update({elm : atom[10]})
                    typelabelList.append(atom[0])
            atomlabelDict.append(typelabelList)

        for atom in atomAllList:
            atomlabelList.append([atom[0],atom_type_dict[atom[1]]])
        
        massList = list(massDict.values())
        atom_type_list = list(atomtypesDict.values())
        name_list = list(atomnameDict.values())
        elementList = list(elementDict.values())
        radiusList = list(radiusDict.values())
        chargesList = list(chargeDict.values())
        bfactorList = list(bfactorDict.values())

        interNum = 0
        interDict = {}
        interTypeDict = {}
        interTypeKind = {}
        interTypeParams = {}
        for b in topo.bondtypes:
            bondList = []
            typeList = []
            noninter = True
            bc = 0 
            for i1, a1 in enumerate(topo.atomtypes):
                for i2, a2 in enumerate(topo.atomtypes):
                    if i2<i1:
                        continue
                    if((b.atom1.name == a1.name and b.atom2.name == a2.name) or
                       (b.atom1.name == a2.name and b.atom2.name == a1.name)):
                        reslist=[a1.res_id, a2.res_id]
                        bondList.append([i1,i2])
                        interDict.update({interNum : bondList})
                        if noninter:
                            noninter = False
                            typeList.extend(list([
                                atom_type_dict[a1.unique_name],
                                atom_type_dict[a2.unique_name]
                                ]))
                            interTypeDict.update({interNum : typeList})
                            if b.constraint:
                                interTypeKind.update({interNum : "bond-constraint"})
                            else:
                                if b.involve_h:
                                    interTypeKind.update({interNum : "H-bond"})
                                else:
                                    interTypeKind.update({interNum : "bond"})
                            interTypeParams.update({interNum : {
                                "quartic_force":b.quartic_force,
                                "harmonic_force":b.harmonic_force,
                                "bond_length":b.eq_length,
                                }})
            interNum += 1

#        for b in topo.ljparameters:
#            bondList = []
#            typeList = []
#            noninter = True
#            bc = 0 
#            for i1, a1 in enumerate(topo.atomtypes):
#                for i2, a2 in enumerate(topo.atomtypes):
#                    if((b.type1_name == a1.name and b.type2_name == a2.name) or
#                       (b.type1_name == a2.name and b.type2_name == a1.name)):
#                        reslist=[a1.res_id, a2.res_id]
#                        bondList.append([i1,i2])
#                        interDict.update({interNum : bondList})
#                        if noninter:
#                            noninter = False
#                            typeList.extend(list([
#                                atom_type_dict[a1.unique_name],
#                                atom_type_dict[a2.unique_name]
#                                ]))
#                            interTypeDict.update({interNum : typeList})
#                            interTypeKind.update({interNum : "VDW-LJ-nonbonded"})
#                            interTypeParams.update({interNum : {
#                                "c12":b.lj_r12,
#                                "c6":b.lj_r6,
#                                "c12_lj14":b.lj14_r12,
#                                "c6_lj14":b.lj14_r6
#                                }})
#            interNum += 1

        topoStor.topoDict.update({"system_name" : system_name})
        topoStor.topoDict.update({"atom_name_list" : atom_name_list})
        topoStor.topoDict.update({"atom_label_list" : atomlabelList})
        topoStor.topoDict.update({"name_list" : name_list})
        topoStor.topoDict.update({"types_list" : types_list})
        topoStor.topoDict.update({"atom_type_list" : atom_type_list})
        topoStor.topoDict.update({"atom_mass_list" : massList})
        topoStor.topoDict.update({"atom_element_list" : atom_element_list})
        topoStor.topoDict.update({"element_list" : elementList})
        topoStor.topoDict.update({"atom_radius_list" : radiusList})
        topoStor.topoDict.update({"atom_charge_list" : chargesList})
        topoStor.topoDict.update({"interactions_dict" : interDict})
        topoStor.topoDict.update({"interactions_type_dict" : interTypeDict})
        topoStor.topoDict.update({"interactions_kind_dict" : interTypeKind})
        topoStor.topoDict.update({"interactions_params_dict" : interTypeParams})
        topoStor.topoDict.update({"mol_list" : residueList})
        #topoStor.topoDict.update({"atom_to_mol_list" : atom_to_residue})
        topoStor.topoDict.update({"residue_list" : residueList})
    return topoStor

def pymConvertTopoDict(topoStor, topoPYM):
    """ Function to convert Pymolfile topology info
        to MDDataAccess

        Pymolfile stores structure data in numpy with following fields
            (name, type, resname, resid, segid, chain, altloc, insertion, 
             occupancy, bfactor, mass, charge, radius, atomicnumber)
        Ex.: ['N' 'NH3' 'ASP' '48' 'PRO1' 'P' ' ' '' 
              '0.0' '0.0' '14.007' '-0.3' '0.0' '0']
    """
    topo=None
    if isinstance(topoPYM, pym.OpenMolfile):
        topo = topoPYM.topology

        def getatomall(atom):
            atmid=0
            atmname = atom[0]
            atmtyp = atom[1]
            if atom[10]:
                atmmass = float(atom[10])
            else:
                atmmass = None
            if atom[3]:
                atmresid = atom[3]
                atmres = atom[2]
            else:
                atmresid = None
                atmres = None
            if atom[4]:
                atmsegid = atom[4]
            else:
                atmsegid = None
            if atom[11]:
                atmchrg = float(atom[11])
            else:
                atmchrg = None
            if atom[12]:
                atmrad = float(atom[12])
            else:
                atmrad = None
            if atom[9]:
                atmbfac = float(atom[9])
            else:
                atmbfac = None
            atm_unique = ''
            if atmname:
                atm_unique = atm_unique + atmname 
            if atmtyp:
                atm_unique = atm_unique + "-" + atmtyp
            #if atmres:
            #    atm_unique = atm_unique + "-" + atmres
            #if atmsegid:
            #    atm_unique = atm_unique + "-" + atmsegid
            return [atmid,atm_unique,atmname,atmtyp,atmres,atmresid,atmsegid,atmmass,atmchrg,atmrad,atmbfac]

        def checkatomsdiff(atom1,atom2,checklist):
            index=0
            target=len(checklist)
            for t in checklist:
                if atom1[t] is not atom2[t]:
                    index += 1
            if index==target:
                return True
            else:
                return False

        def atom_respairs(atom1,atom2):
            return [get_atomres(atom1),get_atomres(atom2)]

        def atompairs(atom1,atom2):
            return [getatom(atom1),getatom(atom2)]

        def atompairids(atom1,atom2):
            return [atom1[0],atom2[0]]

        #chainList = [a[5] if a[5] else '' for a in topo.structure]
        #segmentList = [a[4] if a[4] else '' for a in topo.structure]
        #residList = [a[3] if a[3] else None for a in topo.structure]
        #residueList = [a[2] if a[2] else '' for a in topo.structure]
        chainList = [a[5] for a in topo.structure]
        segmentList = [a[4] for a in topo.structure]
        residList = [a[3] for a in topo.structure]
        residueList = [a[2] for a in topo.structure]
        atomList = topo.structure
        atomAllList = [getatomall(a) for a in topo.structure]
        count=0
        for i in range(len(atomAllList)):
            atomAllList[i][0]=count
            count+=1

        types_list = list(set([a[1] for a in atomAllList]))
        atom_name_list = [a[0] for a in topo.structure]
        #atom_type_list = [a[3] for a in atomList]


        atom_element_list = []
        for atom in topo.structure:
            if atom[1] in ELEMENTS_MASS_TABLE.keys():
                element = atom[1]
            else:
                element = get_element_name(atom[1])
                if element is None:
                    element = atom[1]
            atom_element_list.append(element) 

        system_name = ''
        if atom_element_list:
            system_name = system_name + ''.join([
                el + str(atom_element_list.count(el)) for el in set(
                    atom_element_list)])
        if segmentList is not None:
            system_name = system_name + '-' + '-'.join([
                str(seg) for seg in list(
                    set(segmentList)) if seg is not None])
        if residueList is not None:
            system_name = system_name + '-' + '-'.join([
                str(res) for res in list(
                    set(residueList)) if res is not None])

        atom_all_list = np.asarray(atomAllList)
        atom_names = atom_all_list[:,2]
        atom_types = atom_all_list[:,3]
        atom_masses = atom_all_list[:,7]
        atom_charges = atom_all_list[:,8]
        atom_radiuses = atom_all_list[:,9]
        atom_bfactors = atom_all_list[:,10]
        #print("Atom Bonds:",topo.bonds)
        #print("Atom Angles:",topo.angles)
        #print("Atom Dihedrals:",topo.dihedrals)
        #print("Atom Impropers:",topo.impropers)


        atom_type_dict = {}
        atomtypesDict = {}
        atomnameDict = {}
        massDict = {}
        elementDict = {}
        radiusDict = {}
        chargeDict = {}
        bfactorDict = {}
        atomlabelList = []
        atomlabelDict = []
        for ielm in range(len(types_list)):
            elm = types_list[ielm]
            atom_type_dict.update({elm : ielm+1})
            typelabelList = []
            for atom in atomAllList:
                if elm == atom[1]:
                    if atom[2] in ELEMENTS_MASS_TABLE.keys():
                        element = atom[2]
                    else:
                        element = get_element_name(atom[2])
                        if element is None:
                            element = atom[1]
                    elementDict.update({elm : element})
                    atomnameDict.update({elm : atom[2]})
                    atomtypesDict.update({elm : atom[3]})
                    massDict.update({elm : atom[7]})
                    chargeDict.update({elm : atom[8]})
                    radiusDict.update({elm : atom[9]})
                    bfactorDict.update({elm : atom[10]})
                    typelabelList.append(atom[0])
            atomlabelDict.append(typelabelList)

        for atom in atomAllList:
            atomlabelList.append([atom[0],atom_type_dict[atom[1]]])
        
        massList = list(massDict.values())
        atom_type_list = list(atomtypesDict.values())
        name_list = list(atomnameDict.values())
        elementList = list(elementDict.values())
        radiusList = list(radiusDict.values())
        chargesList = list(chargeDict.values())
        bfactorList = list(bfactorDict.values())

        topbList = np.column_stack((topo.bonds["from"],topo.bonds["to"]))
        topbNames = []
        for pair in topbList:
            topbNames.append(atomAllList[pair[0]-1][1] + '-' + atomAllList[pair[1]-1][1])

        topb=list(set(topbNames))

        interNum = 0
        interDict = {}
        interTypeDict = {}
        for tb in topb:
            bondList = []
            typeList = []
            noninter = True
            bc = 0 
            for b in topbList:
                topt=atomAllList[b[0]-1][1] + '-' + atomAllList[b[1]-1][1]
                if topt == tb:
                    reslist=[atomAllList[b[0]-1][4],atomAllList[b[1]-1][4]]
                    atmlist=[atomAllList[b[0]-1],atomAllList[b[1]-1]]
                    atmidlist=[b[0]-1,b[1]-1]
                    bondList.append(atmidlist)
                    interDict.update({interNum : bondList})
                    if noninter:
                        noninter = False
                        typeList.extend(list([
                            atom_type_dict[atmlist[0][1]],
                            atom_type_dict[atmlist[1][1]]
                            ]))
                        interTypeDict.update({interNum : typeList})
            interNum += 1

#        for ielm in range(len(atom_type_list)-1):
#            for jelm in range(ielm+1, len(atom_type_list)):
#                aelm = atom_type_list[ielm]
#                belm = atom_type_list[jelm]
#                bondList = []
#                typeList = []
#                bondid = 0
#                noninter = True
#                for key in topdk:
#                    topt = topd[key]
#                    for b in topt:
#                        reslist=atom_respairs(b.atoms[0],b.atoms[1])
#                        atmlist=atompairs(b.atoms[0],b.atoms[1])
#                        atmidlist=atompairids(b.atoms[0],b.atoms[1])
#                        if((aelm == str(atmlist[0][2]) and belm == str(atmlist[1][2])) or 
#                           (aelm == str(atmlist[1][2]) and belm == str(atmlist[0][2]))):
#                            bondList.append(atmidlist)
#                            interDict.update({interNum : bondList})
#                            if noninter:
#                                noninter = False
#                                typeList.extend(list([
#                                    atom_type_dict[aelm],
#                                    atom_type_dict[belm]
#                                ]))
#                                interTypeDict.update({interNum : typeList})
#                                interNum += 1
#                        bondid += 1

        #atomIndex = np.arange(len(residueList))
        #atom_to_residue = np.zeros((len(residueList), 2), dtype=int)
        #atom_to_residue[:,0] = atomIndex+1
        #atom_to_residue[:,1] = np.array(residueList)+1
       

        topoStor.topoDict.update({"system_name" : system_name})
        topoStor.topoDict.update({"atom_name_list" : atom_name_list})
        topoStor.topoDict.update({"atom_label_list" : atomlabelList})
        topoStor.topoDict.update({"name_list" : name_list})
        topoStor.topoDict.update({"types_list" : types_list})
        topoStor.topoDict.update({"atom_type_list" : atom_type_list})
        topoStor.topoDict.update({"atom_mass_list" : massList})
        topoStor.topoDict.update({"atom_element_list" : atom_element_list})
        topoStor.topoDict.update({"element_list" : elementList})
        topoStor.topoDict.update({"atom_radius_list" : radiusList})
        topoStor.topoDict.update({"atom_charge_list" : chargesList})
        topoStor.topoDict.update({"interactions_dict" : interDict})
        topoStor.topoDict.update({"interactions_type_dict" : interTypeDict})
        topoStor.topoDict.update({"mol_list" : residueList})
        #topoStor.topoDict.update({"atom_to_mol_list" : atom_to_residue})
        topoStor.topoDict.update({"residue_list" : residueList})
    return topoStor

def mdtConvertTopoDict(topoStor, topoMDT):
    if(isinstance(topoMDT, mdt_Topology) or 
       isinstance(topoMDT, mdt_Trajectory)):
        topo = topoMDT

        topologyTable, topologyBonds = topo.to_dataframe()
        topologyDict = topologyTable.to_dict(orient='list')
        types_list=list(set(topologyDict["name"]))
        atom_type_dict = {}
        massesDict = {}
        elementDict = {}
        radiusDict = {}
        for ielm in range(len(types_list)):
            elm = types_list[ielm]
            atom_type_dict.update({elm : ielm+1})
            for atom in topo.atoms:
                if elm == atom.name:
                    massesDict.update({atom.name : atom.element.mass})
                    elementDict.update({atom.name : atom.element.symbol})
                    radiusDict.update({atom.name : atom.element.radius})

        massesList = list(massesDict.values())
        try:
            atom_type_list = list(topologyDict["name"].values())
        except AttributeError:
            atom_type_list = topologyDict["name"]
        atom_element_list = [atom.element.symbol for atom in topo.atoms]
        elementList = list(elementDict.values())
        radiusList = list(radiusDict.values())

        interNum = 0
        interDict = {}
        interTypeDict = {}
        for ielm in range(len(types_list)-1):
            for jelm in range(ielm+1, len(types_list)):
                aelm = types_list[ielm]
                belm = types_list[jelm]
                bondList = []
                typeList = []
                bondid = 0
                noninter = True
                for bond in topo.bonds:
                    molname1, molatom1 = str(bond[0]).split('-')
                    molname2, molatom2 = str(bond[1]).split('-')
                    if((aelm == str(molatom1) and belm == str(molatom2)) or 
                        (aelm == str(molatom2) and belm == str(molatom1))):
                        bondList.append(list(topologyBonds[bondid]))
                        interDict.update({interNum : bondList})
                        if noninter:
                            noninter = False
                            typeList.extend(list([
                                atom_type_dict[aelm],
                                atom_type_dict[belm]
                                ]))
                            interTypeDict.update({interNum : typeList})
                            interNum += 1
                    bondid += 1

        residueList = topologyDict["resSeq"]
        atomIndex = np.arange(len(residueList))
        atom_to_residue = np.zeros((len(residueList), 2), dtype=int)
        atom_to_residue[:,0] = atomIndex+1
        atom_to_residue[:,1] = np.array(residueList)+1
        
        system_name = ''
        if atom_element_list:
            system_name = system_name + ''.join([
                el + str(atom_element_list.count(el)) for el in set(
                    atom_element_list)])
        if residueList is not None:
            system_name = system_name + '-' + '-'.join([
                str(res) for res in list(
                    set(topologyDict["resName"])) if res is not None])

        topoStor.topoDict.update({"system_name" : system_name})
        topoStor.topoDict.update({"types_list" : types_list})
        topoStor.topoDict.update({"atom_type_list" : atom_type_list})
        topoStor.topoDict.update({"atom_mass_list" : massesList})
        topoStor.topoDict.update({"atom_element_list" : atom_element_list})
        topoStor.topoDict.update({"element_list" : elementList})
        topoStor.topoDict.update({"atom_radius_list" : radiusList})
        topoStor.topoDict.update({"atom_charge_list" : []})
        topoStor.topoDict.update({"interactions_dict" : interDict})
        topoStor.topoDict.update({"interactions_type_dict" : interTypeDict})
        topoStor.topoDict.update({"mol_list" : residueList})
        topoStor.topoDict.update({"atom_to_mol_list" : atom_to_residue})
        topoStor.topoDict.update({"residue_list" : residueList})
    return topoStor

def mdaConvertTopoDict(topoStor, topoMDA):
    if(isinstance(topoMDA, mda_u.Universe) or 
       isinstance(topoMDA, mda_u.Topology)):
        topo = topoMDA

        def getseg(seg):
            regex = re.compile(r"^seg_[0-9]+_b")
            if regex.findall(seg.segid):
                segname = bytes(re.sub(r"^seg_[0-9]+_b","",seg.segid).replace("'",""), "utf-8").decode("utf-8")
            else:
                segname = seg.segid
            return [seg.ix,segname]
        
        def getres(res):
            regex = re.compile(r"^b'")
            if isinstance(res.resname, bytes):
                resname = res.resname.decode("utf-8")
            else:
                resname= res.resname
            if regex.findall(resname):
                resname = bytes(re.sub(r"^b","",resname).replace("'",""), "utf-8").decode("utf-8")
            return [res.ix,resname]

        def get_atomseg(atom):
            return [atom.segment.ix,atom.segment.segname.decode('utf-8')]
        
        def get_atomres(atom):
            return [atom.resname.decode('utf-8'),atom.resid]

        def getatom(atom):
            atmid = atom.ix
            atmname = atom.name.decode('utf-8')
            atmtyp = atom.type.decode('utf-8')
            atmres = atom.resname.decode('utf-8')
            atmresid = atom.resid
            atmsegid = atom.segid
            atm_unique = atmname + "-" + atmtyp
            return [atmid,atm_unique,atmname,atmtyp,atmres,atmresid,atmsegid]

        def getatomall(atom):
            atmid = atom.ix
            atmname = atom.name.decode('utf-8')
            atmtyp = atom.type.decode('utf-8')
            if hasattr(atom, "mass"):
                atmmass = float(atom.mass)
            else:
                atmmass = None
            if hasattr(atom, "resid"):
                atmresid = atom.resid
                atmres = atom.resname.decode('utf-8')
            else:
                atmresid = None
                atmres = None
            if hasattr(atom, "segid"):
                atmsegid = atom.segid
            else:
                atmsegid = None
            if hasattr(atom, "charge"):
                atmchrg = float(atom.charge)
            else:
                atmchrg = None
            if hasattr(atom, "radius"):
                atmrad = float(atom.radius)
            else:
                atmrad = None
            if hasattr(atom, "bfactor"):
                atmbfac = float(atom.bfactor)
            else:
                atmbfac = None
            atm_unique = atmname + "-" + atmtyp
            return [atmid,atm_unique,atmname,atmtyp,atmres,atmresid,atmsegid,atmmass,atmchrg,atmrad,atmbfac]

        def checkatomsdiff(atom1,atom2,checklist):
            index=0
            target=len(checklist)
            for t in checklist:
                if atom1[t] is not atom2[t]:
                    index += 1
            if index==target:
                return True
            else:
                return False

        def atom_respairs(atom1,atom2):
            return [get_atomres(atom1),get_atomres(atom2)]

        def atompairs(atom1,atom2):
            return [getatom(atom1),getatom(atom2)]

        def atompairids(atom1,atom2):
            return [getatom(atom1)[0],getatom(atom2)[0]]

        segmentList = [getseg(a)[1] for a in topo.segments]
        residueList = [getres(a)[1] for a in topo.residues]
        atomList = [getatomall(a) for a in topo.atoms]
        #atomAllList = [getatomall(a) for a in topo.atoms]

        types_list = list(set([a[1] for a in atomList]))
        atom_name_list = [a[2] for a in atomList]
        #atom_type_list = [a[3] for a in atomList]


        atom_element_list = []
        for atom in atomList:
            try:
                guessed_element = mda.topology.guessers.guess_atom_element(atom[2])
            except (TypeError, ValueError, AttributeError):
                guessed_element = atom[2]
            if guessed_element in ELEMENTS_MASS_TABLE.keys():
                element = guessed_element
            else:
                element = get_element_name(atom[2])
                #if element is None:
                #    element = get_element_name(atom[3])
                if element is None:
                    element = atom[2]
            atom_element_list.append(element) 

        system_name = ''
        if atom_element_list:
            system_name = system_name + ''.join([
                el + str(atom_element_list.count(el)) for el in set(
                    atom_element_list)])
        if segmentList is not None:
            system_name = system_name + '-' + '-'.join([
                str(seg) for seg in list(
                    set(segmentList)) if seg is not None])
        if residueList is not None:
            system_name = system_name + '-' + '-'.join([
                str(res) for res in list(
                    set(residueList)) if res is not None])

        attrlist = dir(topo.atoms)
        atom_names = []
        atom_types = []
        atom_masses = []
        atom_radiuses = []
        atom_bfactors = []
        atom_charges = []
        if "names" in attrlist:
            atom_names = topo.atoms.names
        if "types" in attrlist:
            atom_types = topo.atoms.types
        if "masses" in attrlist:
            atom_masses = topo.atoms.masses
        if "charges" in attrlist:
            atom_charges = topo.atoms.charges
        if "radiuses" in attrlist:
            atom_radiuses = topo.atoms.radiuses
        if "bfactors" in attrlist:
            atom_bfactors = topo.atoms.bfactors
        #print("Atom Bonds:",MDdata.topohandler.atoms.bonds)
        #print("Atom Angles:",MDdata.topohandler.atoms.angles)
        #print("Atom Dihedrals:",MDdata.topohandler.atoms.dihedrals)
        #print("Atom Impropers:",MDdata.topohandler.atoms.impropers)


        atom_type_dict = {}
        atomtypesDict = {}
        atomnameDict = {}
        massDict = {}
        elementDict = {}
        radiusDict = {}
        chargeDict = {}
        bfactorDict = {}
        atomlabelList = []
        atomlabelDict = []
        for ielm in range(len(types_list)):
            elm = types_list[ielm]
            atom_type_dict.update({elm : ielm+1})
            typelabelList = []
            for atom in atomList:
                if elm == atom[1]:
                    try:
                        guessed_element = mda.topology.guessers.guess_atom_element(atom[2])
                    except (TypeError, ValueError, AttributeError):
                        guessed_element = atom[2]
                    if guessed_element in ELEMENTS_MASS_TABLE.keys():
                        element = guessed_element
                    else:
                        element = get_element_name(atom[2])
                        #if element is None:
                        #    element = get_element_name(atom[3])
                        if element is None:
                            element = atom[2]
                    elementDict.update({elm : element})
                    atomnameDict.update({elm : atom[2]})
                    atomtypesDict.update({elm : atom[3]})
                    massDict.update({elm : atom[7]})
                    chargeDict.update({elm : atom[8]})
                    radiusDict.update({elm : atom[9]})
                    bfactorDict.update({elm : atom[10]})
                    typelabelList.append(atom[0])
            atomlabelDict.append(typelabelList)

        for atom in atomList:
            atomlabelList.append([atom[0],atom_type_dict[atom[1]]])
        
        massList = list(massDict.values())
        atom_type_list = list(atomtypesDict.values())
        name_list = list(atomnameDict.values())
        elementList = list(elementDict.values())
        radiusList = list(radiusDict.values())
        chargesList = list(chargeDict.values())
        bfactorList = list(bfactorDict.values())

        topd = topo.bonds.topDict
        topdk = topd.keys()

        interNum = 0
        interDict = {}
        interTypeDict = {}
        for key in topdk:
            bondList = []
            typeList = []
            noninter = True
            topt = topd[key]
            for b in topt:
                reslist=atom_respairs(b.atoms[0],b.atoms[1])
                atmlist=atompairs(b.atoms[0],b.atoms[1])
                atmidlist=atompairids(b.atoms[0],b.atoms[1])
                bondList.append(atmidlist)
                interDict.update({interNum : bondList})
                if noninter:
                    noninter = False
                    typeList.extend(list([
                        atom_type_dict[atmlist[0][1]],
                        atom_type_dict[atmlist[1][1]]
                    ]))
                    interTypeDict.update({interNum : typeList})
            interNum += 1

#        for ielm in range(len(atom_type_list)-1):
#            for jelm in range(ielm+1, len(atom_type_list)):
#                aelm = atom_type_list[ielm]
#                belm = atom_type_list[jelm]
#                bondList = []
#                typeList = []
#                bondid = 0
#                noninter = True
#                for key in topdk:
#                    topt = topd[key]
#                    for b in topt:
#                        reslist=atom_respairs(b.atoms[0],b.atoms[1])
#                        atmlist=atompairs(b.atoms[0],b.atoms[1])
#                        atmidlist=atompairids(b.atoms[0],b.atoms[1])
#                        if((aelm == str(atmlist[0][2]) and belm == str(atmlist[1][2])) or 
#                           (aelm == str(atmlist[1][2]) and belm == str(atmlist[0][2]))):
#                            bondList.append(atmidlist)
#                            interDict.update({interNum : bondList})
#                            if noninter:
#                                noninter = False
#                                typeList.extend(list([
#                                    atom_type_dict[aelm],
#                                    atom_type_dict[belm]
#                                ]))
#                                interTypeDict.update({interNum : typeList})
#                                interNum += 1
#                        bondid += 1

        #atomIndex = np.arange(len(residueList))
        #atom_to_residue = np.zeros((len(residueList), 2), dtype=int)
        #atom_to_residue[:,0] = atomIndex+1
        #atom_to_residue[:,1] = np.array(residueList)+1
       

        topoStor.topoDict.update({"system_name" : system_name})
        topoStor.topoDict.update({"atom_name_list" : atom_name_list})
        topoStor.topoDict.update({"atom_label_list" : atomlabelList})
        topoStor.topoDict.update({"name_list" : name_list})
        topoStor.topoDict.update({"types_list" : types_list})
        topoStor.topoDict.update({"atom_type_list" : atom_type_list})
        topoStor.topoDict.update({"atom_mass_list" : massList})
        topoStor.topoDict.update({"atom_element_list" : atom_element_list})
        topoStor.topoDict.update({"element_list" : elementList})
        topoStor.topoDict.update({"atom_radius_list" : radiusList})
        topoStor.topoDict.update({"atom_charge_list" : chargesList})
        topoStor.topoDict.update({"interactions_dict" : interDict})
        topoStor.topoDict.update({"interactions_type_dict" : interTypeDict})
        topoStor.topoDict.update({"mol_list" : residueList})
        #topoStor.topoDict.update({"atom_to_mol_list" : atom_to_residue})
        topoStor.topoDict.update({"residue_list" : residueList})
    return topoStor

def aseConvertTopoDict(topoStor, topoASE):
    pass

class MDDataConverter(object):
    def __init__(self):
        pass

    def topology_to_dictionary(self):
        """ This function generates self.topologyDict dictionary
            and converts/stores all the topology information 
            according to the meta info definitions
        """
        newDict = {}
        return newDict

    def topology_num_topo_mol(topodict, itemdict):
        """ Function to generate data for number_of_topology_molecules
        """
        if "residue_list" in topodict:
            residueList = topodict["residue_list"]
            return False, len(residueList), itemdict
        else:
            return False, None, itemdict

    def topology_system_name(topodict, itemdict):
        """ Function to generate data for system_name
        """
        if "filename" in topodict:
            filename = topodict["filename"]
            fdir, fbase, fext = get_dir_base_extension(filename)
            fdir = fdir.split('/')[-1]
        else:
            fdir = ''
            fbase = ''
        if "system_name" in topodict:
            system_name = topodict["system_name"]
        else:
            system_name = ''
        system_name = system_name + '-' + fdir + '-' + fbase
        return False, system_name, itemdict

    def topology_atom_to_mol(topodict, itemdict):
        """ Function to generate data for atom_to_molecule 
        """
        if "residue_list" in topodict:
            atom_to_residue = topodict["residue_list"]
            return False, atom_to_residue, itemdict
        else:
            return False, None, itemdict

class MDDataBase(dict):
    """ Stores user data in MDDataAccess
        classes.
    """
    def __init__(self, *args, **kwargs):
        super(MDDataBase, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if k in self:
                        self[k] = v 
        if kwargs:
            for k, v in kwargs.items():
                if k in self:
                    self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(MDDataBase, self).__setitem__(key, value)
        self.__dict__.update({key: value})

class MDDataCommon(MDDataBase):
    """ Stores the base data of MDDataAccess
        classes.

        time is (1 x n_step) numpy array.
        n_ variables are in numbers.
        iread is expected to be a generator 
            that loops all n_steps of given array.
    """
    def __init__(self,
        step_no = 0,
        iread = None,
        n_steps = None,
        time = None
        ):
        pass

class MDDataTopology(MDDataCommon):
    """ Stores the topology data in a conveniate 
        way that individual components of topology 
        data can be accessed as numpy arrays or 
        dictionaries.
        
        n_ values are numbers
        _types values are numbers in numpy array
        _labels values are a list of strings

        if not stated, all variables are in numpy arrays.
    """
    def __init__(self,
        system_name = "",
        _topology = {},
        atoms = None,
        bonds = None,
        dihedrals = [],
        impropers = [],
        residues = [],
        n_atoms = None,
        _numAtoms = None,
        numAtoms = None,
        n_bonds = None,
        n_dihedrals = None,
        n_impropers = None,
        n_residues = None,
        atom_types = None,
        atom_labels = None,
        bond_types = None,
        bond_labels = None,
        dihedral_types = None,
        dihedral_labels = None,
        improper_types = None,
        improper_labels = None,
        residue_types = None,
        residue_labels = None,
        atom_in_mol = None,
        charmmcoor = None,
        charmmpsf = None,
        charmmcoortopo = None,
        topoDict = {}
        ):
        pass

class MDDataTrajectory(MDDataCommon):
    """ Stores the trajectory data as well as 
        velocities and forces, if they exist.
    """
    def __init__(self,
        title = None,
        natoms = None,
        positions = None,
        forces = None,
        time = None,
        velocities = None,
        unitcell_vectors = None,
        unitcell_lengths = None,
        unitcell_angles = None,
        trajDict = None,
        frame_no = None
        ):
        pass

class MDDataThermo(MDDataCommon):
    """ Stores the thermostat values (energies) 
        data as well as, if it exists.

        _tensor variables are Nx3x3 numpy arrays 
            that hold 3x3 tensors of all N steps (n_steps).
    """
    def __init__(self,
        nsteps = None,
        step_no = None,
        time = None,
        ljsr = None,
        coulombsr = None,
        pe = None,
        ke = None,
        te = None,
        press = None,
        temp = None,
        bond = None,
        ub = None,
        proper = None,
        improper = None,
        cmapdih = None,
        lj = None,
        coulomb = None,
        conserveden = None,
        constrrmsd = None,
        virial_tensor = None,
        pressure_tensor = None,
        thermoDict = None
        ):
        pass

class MDDataAccess(object):
    """Data access class for various different atomic coordinate files.

    Reading topology and trajectory files using one interface 
    for different MD data analyses packages.

    In current version, the following packages are supported:
        mdtraj       : MDtraj 1.8.0 (9 Nov 2016) and 1.9.0 (3 Sep 2017)
        mdanalysis   : MDAnalysis 0.16.2 (27 Jun 2017) and 0.17.1 (4 Sep 2017)
        ase          : ASE 3.14.1 (28 Jun 2017)
        usersupplied : self.UserSuppliedInterface
    
    If user does not supply a package name or an ordering list for the 
    interfaces using the names at first column above, the packages will 
    be used for accessing data in given files with the order above.

    A user defined function can also be supplied to the high-level 
    data access in MDDataAccess. (See self.user_formats).

    Additional Interfaces:
        A common interface to access thermodynamical quantities for GROMACS
        is added through panedr package. To initiate the file handler use 
        thermoFileHandler function and to access the values at each step 
        use thermo_iread which returns a dictionary with all properties.
        For a full list of properties see docs of thermo_iread.
        
        Interface is supplied by the following package:
          Panedr 0.2 (15 Jan 2016)

    Returns topology and/or trajectory data as numpy arrays if it exists.
    """
    def __init__(self):
        try:
            if self.initialized:
                pass
        except (AttributeError, ValueError, TypeError):
            self.initialize()
            self.initialized = True
            self.user_formats = {
                "myformat":    {"My Fancy Molecule Format", self.custom_iread}
            }

    def initialize(self):
        self.natoms = None
        self.topology = None # Main storage for topology data
        self.trajectory = None # Main storage for trajectory data
        self.inputcoords = None
        self.outputcoords = None
        self.trajtype = None
        self.forcefield = None # Main storage for force field parameters
        self.thermostats = None # Main storage for thermodynamical quantities 
                                   #  and properties such as energies, temperatures
        self.forcefieldhandler = None # Main storage for force field parameters
        self.access_ui = None
        self.interfaceorder = None
        self.interfacematch = None
        self.UserSuppliedInterface = None
        self.set_defaults()
        self.init_topo()
        self.init_traj()
        self.init_incoord()
        self.init_outcoord()
        self.init_thermo()

    def set_defaults(self):
        self.trajchunk = 100 # The chunk size for reading trajectory file (only supports MDtraj for now)
        self.readfirst = False

    def init_topo(self):
        self.topofile = None # File name of the topology file with path if needed
        self.topoformat = None # Format type of topology file
        self.topoplugin = None # Format type of topology file
        self.topocode = None # To explicitly define the parsing library (pymolfile, mdtraj, ASE ...) for topology files
        self.topohandler = None # The object parsing the topology file
        self.topo_n_atoms = None # Number of atoms at topology file
        self.topostream = None # List of lines in a text stream of topology
        self.topocharmmpsf = None # List of lines in a text stream of topology
        self.topology = None

    def init_traj(self):
        self.trajformat = None # Format type of trajectory file
        self.trajhandler = None # The object parsing the trajectory file
        self.trajplugin = None # The object parsing the trajectory file
        self.trajiter = None # The object parsing the trajectory file
        self.trajcode = None # To explicitly define the parsing library (pymolfile, mdtraj, ASE ...) for trajectory files
        self.trajfile = None # File name of the trajectory file with path if needed
        self.trajstream = None # List of lines in a text stream of input trajectory 
        self.trajectory = None

    def init_incoord(self):
        self.incoord_natoms = None
        self.incoordformat = None # Format type of trajectory file
        self.incoordhandler = None # The object parsing the trajectory file
        self.incoordplugin = None # The object parsing the trajectory file
        self.incoorditer = None # The object parsing the trajectory file
        self.incoordcode = None # To explicitly define the parsing library (pymolfile, mdtraj, ASE ...) for trajectory files
        self.incoordfile = None # File name of the trajectory file with path if needed
        self.incoordstream = None # List of lines in a text stream of input trajectory 
        self.inputpositions = None

    def init_outcoord(self):
        self.outcoord_natoms = None
        self.outcoordformat = None # Format type of output file
        self.outcoordhandler = None # The object parsing the output file
        self.outcoordplugin = None # The plugin parsing the output file
        self.outcoorditer = None # The iteration at the output file
        self.outcoordcode = None # To explicitly define the parsing library (pymolfile, mdtraj, ASE ...) for the files
        self.outcoordfile = None # File name of the file with path if needed
        self.outcoordstream = None # List of lines in a text stream of input trajectory 
        self.outputpositions = None

    def init_thermo(self):
        self.thermofile = None # File name of the thermostat file with path if needed
        self.thermoformat = None # Format type of topology file
        self.thermocode = None # To explicitly define the parsing library (Panedr or ??) for topology files
        self.thermohandler = None # The object parsing the thermostat file (Ex.: ener.edr of Gromacs)
        self.thermostats = None

    def reset_all(self):
        self.forcefield = None
        self.set_defaults()
        self.init_topo()
        self.init_traj()
        self.init_incoord()
        self.init_outcoord()
        self.init_thermo()

    def set_topo(self, filename, fileformat=None):
        self.topofile = filename
        self.topoformat = fileformat

    def set_traj(self, filename, fileformat=None):
        self.trajfile = filename
        self.trajformat = fileformat

    def set_incoord(self, filename, fileformat=None):
        self.incoordfile = filename
        self.incoordformat = fileformat

    def set_outcoord(self, filename, fileformat=None):
        self.outcoordfile = filename
        self.outcoordformat = fileformat

    def set_thermo(self, filename, fileformat=None):
        self.thermofile = filename
        self.thermoformat = fileformat

    def load(self, interface=None):
        """Loads the file handles for trajectory and/or topology
        """
        
        #if self.topohandler is None:
        #    if self.topofile:
        self.check_topology_format_support()

        if self.trajhandler is None:
            self.check_trajectory_format_support("traj")

        return self.trajhandler
    
    def load_topology(self):
        """Loads the file handles for topology only
        """
        
        #if self.topohandler is None:
        #    if self.topofile:
        #        self.check_topology_format_support()
        #else:
        if self.topofile is None:
            self.topofile = ''
        self.check_topology_format_support()

        return self.topohandler
    
    def load_incoord(self):
        """Loads the file handles for coordinates only
        """
        
        if self.incoordhandler is None:
        #    if self.incoordfile:
            self.check_trajectory_format_support("input")

        return self.incoordhandler
    
    def load_outcoord(self):
        """Loads the file handles for coordinates only
        """
        
        if self.outcoordhandler is None:
            if self.outcoordfile:
                self.check_trajectory_format_support("output")

        return self.outcoordhandler
    
    def load_thermo(self):
        """Loads the file handles for topology only
        """
        
        if self.thermohandler is None:
            if self.thermofile:
                self.thermohandler = panedr.edr_to_df(self.thermofile, verbose=False)

        if self.thermohandler is not None:
            self.thermocode = "panedr"
            self.thermoformat = ".edr"

        return self.thermohandler

    def get_dir_base_extension(self, file_name):
        """ Splits directory, file base and file extensions

            Returns: directory without leading '/', 
            file base name, and file extension without '.'
        """
        file_base, file_extension_with_dot = os.path.splitext(os.path.basename(file_name))
        file_extension = file_extension_with_dot.split(".")[-1]
        file_dir = os.path.dirname(file_name)
        return file_dir, file_base, file_extension
    
    def get_file_format(self, filename, givenformat):
        """Returns format of a file from extension
        """
        file_format = None
        if givenformat is not None:
            # MDTraj expects that file name extensions start with dot
            if '.' not in givenformat:
                file_format = "." + str(givenformat).lower()
            else:
                file_format = str(givenformat).lower()
        else:
            filedir, filebase, fileext  = self.get_dir_base_extension(filename)
            if '.' not in fileext:
                file_format = "." + fileext.lower()
            else:
                file_format = fileext.lower()
        return file_format

    def is_class_of_module(self, obj, mod, mods=None):
        """ Check if the obj1 is a class of obj2 module
            or its modules.

            Return: True if it is, else False
        """
        isinlist = False
        if getattr(obj, '__module__', None) == mod.__name__:
            isinlist = True
        else:
            isinlist = False
        if (mods is not None and isinlist is False):
            for submod in mods:
                if getattr(obj, '__module__', None) == submod.__name__:
                    isinlist = True
                    break
        if isinlist:
            return True
        else:
            return False

    def check_topology_format_support(self):
        """Check if the given format is supported.
        """
        topofilename=None
        if self.topofile:
            topofilename = os.path.basename(self.topofile)
        if self.topoformat is None and topofilename is not None:
            fileloadformat = self.get_file_format(topofilename, self.topoformat)
            self.topoformat = fileloadformat
        else:
            fileloadformat = self.topoformat

        usedefault=True
        # Use the given order to check topology
        if self.interfaceorder:
            for interface in self.interfaceorder:
                if "gromosread" in interface:
                    zipfile = None
                    ziptype = None
                    filetopoformat, zipfile, ziptype = get_zipType(self.topofile)
                    ftopoform = re.sub('[.]', '', filetopoformat)
                    if('GROMOSTOP' == ftopoform.upper() or 
                       'GROMOSCNF' == ftopoform.upper(),
                       'CNF' == ftopoform.upper(),
                       'TOP' == ftopoform.upper()):
                        self.topohandler = self.load_gromos_topology(ftopoform, base_topo=self.topohandler)
                    if self.topohandler:
                        usedefault=False
                        self.topocode = "gromosread"
                        break
                if "charmmcoor" in interface:
                    topohandler_check = None
                    charmmcoor_dict = None
                    filetopoformat = re.sub('[.]', '', fileloadformat)
                    if('CHARMMCOR' == filetopoformat.upper() or 
                       'CHARMMCOOR' == filetopoformat.upper() or 
                       #'COR' == file_format.upper() or 
                       #'CRD' == file_format.upper() or 
                       'COOR' == filetopoformat.upper()):
                        if self.topofile:
                            charmmcoor_dict = charmm_coor_reader(self.topofile)
                    if('CHARMMSTRCOR' == filetopoformat.upper() or 
                       'CHARMMSTRCRD' == filetopoformat.upper() or
                       'CHARMMSTREAM' == filetopoformat.upper()):
                        if self.topostream is not None:
                            if 'cor' in self.topostream:
                                charmmcoor_dict = charmm_coor_reader(self.topofile, ftextmem=self.topostream['cor'])
                        else:
                            charmmcoor_dict = charmm_coor_reader(chkfile, ftextmem=None)
                    if charmmcoor_dict:
                        if 'binary' in charmmcoor_dict:
                            if charmmcoor_dict['binary'] is False:
                                topohandler_check = charmmcoor_dict
                    if topohandler_check:
                        if self.topohandler is None:
                            if self.topology is None:
                                self.topology = MDDataTopology()
                                self.topology.charmmcoortopo = self.charmm_coor_toporead(topohandler_check)
                        else:
                            self.topo_n_atoms = topohandler_check['numatoms']
                            if self.topology is None:
                                self.topology = MDDataTopology()
                            self.topology.charmmcoor = topohandler_check
                        if self.topohandler:
                            usedefault=False
                            self.topocode = "charmmcoor"
                            break
                elif "parmed" in interface:
                    filetopoformat = re.sub('[.]', '', fileloadformat)
                    self.topohandler = self.load_parmed_topology(filetopoformat, base_topo=self.topohandler)
                    if self.topohandler:
                        usedefault=False
                        self.topocode = "parmed"
                        break
                elif "pymolfile" in interface:
                    if self.topohandler is None:
                        self.topohandler = self.load_pymolfile_topology(fileloadformat)
                    if self.topohandler:
                        usedefault=False
                        self.topocode = "pymolfile"
                        break
                elif "mdtraj" in interface:
                    if self.topohandler is None:
                        self.topohandler = self.load_mdtraj_topology(fileloadformat)
                    if self.topohandler:
                        usedefault=False
                        self.topocode = "mdtraj"
                        break
                elif "mdanalysis" in interface:
                    if self.topohandler is None:
                        self.topohandler = self.load_mdanalysis_topology(fileloadformat)
                    if self.topohandler:
                        usedefault=False
                        self.topocode = "mdanalysis"
                        break
                elif "ase" in interface:
                    if self.topohandler is None:
                        if self.topofile:
                            self.topohandler = self.load_ase_support(self.topofile, file_format=fileloadformat)
                    if self.topohandler:
                        usedefault=False
                        self.topocode = "ase"
                        break
                elif self.UserSuppliedInterface:
                    if isinstance(self.UserSuppliedInterface, MDDataAccess.UserSuppliedInterface):
                        if self.UserSuppliedInterface.name in interface:
                            if self.topohandler is None:
                                self.topohandler = self.UserSuppliedInterface.topology_support(self.topofile, file_format=fileloadformat)
                            if self.topohandler:
                                usedefault=False
                                self.topocode = self.UserSuppliedInterface.name
                                break

        # If no given order or the items in list does not match with the supported pacakges
        # use the default order in heroistic mode.
        if usedefault:
            if self.topohandler is None:
                # Nothing to lose to be heroistic here.
                self.topohandler = self.load_pymolfile_topology(fileloadformat)
            if self.topohandler is None:
                self.topocode = "pymolfile"

            if self.topohandler is None:
                self.topohandler = self.load_mdtraj_topology(fileloadformat)
                if self.topohandler:
                    self.topocode = "mdtraj"

            # If MDTraj does not have support for the format 
            # or can not load the topology, use MDAnalysis and ASE.
            if self.topohandler is None:
                self.topohandler = self.load_mdanalysis_topology(fileloadformat)
                if self.topohandler:
                    self.topocode = "mdanalysis"

            # Fall back to check ASE support
            if self.topohandler is None:
                ase_support = False
                ase_support = self.get_ase_format_support(fileloadformat)
                # May still have chance that ASE can recognize the 
                # format with its filetype checking function
                if ase_support is None:
                    if self.topofile:
                        ase_support = ase_io.formats.filetype(self.topofile)
                        self.topohandler = self.load_ase_support(self.topofile, file_format=ase_support)
                        self.topocode = "ase"
                else:
                    if self.topofile:
                        self.topohandler = self.load_ase_support(self.topofile, file_format=fileloadformat)
                        self.topocode = "ase"

        # If no success after all attempts return False
        if self.topohandler is None:
            return False
        else:
            return True

    def check_trajectory_format_support(self, filetype):
        """Check if the given format is supported.
        """
        if "input" in filetype:
            chkfile = self.incoordfile
            chkformat = self.incoordformat
            self.trajtype = "input"
        elif "output" in filetype:
            chkfile = self.outcoordfile
            chkformat = self.outcoordformat
            self.trajtype = "output"
        else:
            chkfile = self.trajfile
            chkformat = self.trajformat
            self.trajtype = "traj"

        numatoms=None
        chkfilename = None
        if chkfile:
            chkfilename = os.path.basename(chkfile)
        if chkformat is None:
            fileloadformat = self.get_file_format(chkfile, chkformat)
            chkformat = fileloadformat
            if "input" in filetype:
                self.incoordformat = fileloadformat
            elif "output" in filetype:
                self.outcoordformat = fileloadformat
            else:
                self.trajformat = fileloadformat
        else:
            fileloadformat = chkformat

        usedefault=True
        # Use the given order to check topology
        if self.interfaceorder:
            for interface in self.interfaceorder:
                if "gromosread" in interface:
                    trajhandler_check = None
                    zipfile = None
                    ziptype = None
                    filetrajformat, zipfile, ziptype = get_zipType(chkfile)
                    ftrajform = re.sub('[.]', '', filetrajformat)
                    if('GROMOSCNF' == ftrajform.upper() or 
                       'GROMOSTRC' == ftrajform.upper(),
                       'GROMOSTRV' == ftrajform.upper(),
                       'CNF' == ftrajform.upper(),
                       'TRC' == ftrajform.upper(),
                       'TRV' == ftrajform.upper()):
                        trajhandler_check = gto.GromosTrajectory(trc=chkfile, trc_format=ftrajform)
                    if trajhandler_check:
                        trajhandler = None
                        trajhandler = self.gromos_iread(traj_handler=trajhandler_check)
                        if trajhandler:
                            if self.interfacematch:
                                if interface in self.topocode:
                                    usedefault=False
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "mdtraj"
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "mdtraj"
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "mdtraj"
                                    break
                                else:
                                    trajhandler = None
                            else:
                                usedefault=False
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "gromosread"
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "gromosread"
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "gromosread"
                                break
                if "charmmcoor" in interface:
                    trajhandler_check = None
                    charmmcoor_dict = None
                    filetrajformat = re.sub('[.]', '', fileloadformat)
                    if('CHARMMCOOR' == filetrajformat.upper() or 
                       'CHARMMCOR' == filetrajformat.upper() or 
                       #'COR' == file_format.upper() or 
                       #'CRD' == file_format.upper() or 
                       #'COORBIN' == file_format.upper() or 
                       'COOR' == filetrajformat.upper()):
                        charmmcoor_dict = charmm_coor_reader(chkfile)
                    if('CHARMMSTRCOR' == filetrajformat.upper() or 
                       'CHARMMSTRCRD' == filetrajformat.upper() or
                       'CHARMMSTREAM' == filetrajformat.upper()):
                        if('CHARMMSTR' in filetrajformat.upper() and 
                           self.topostream is not None and 
                           (self.trajstream is None and 
                            self.incoordstream is None and 
                            self.outcoordstream is None )):
                            if 'cor' in self.topostream:
                                charmmcoor_dict = charmm_coor_reader(chkfile, ftextmem=self.topostream['cor'])
                        elif('CHARMMSTR' in filetrajformat.upper() and self.trajstream is not None):
                            if 'cor' in self.trajstream:
                                charmmcoor_dict = charmm_coor_reader(chkfile, ftextmem=self.trajstream['cor'])
                        elif('CHARMMSTR' in filetrajformat.upper() and self.incoordstream is not None):
                            if 'cor' in self.incoordstream:
                                charmmcoor_dict = charmm_coor_reader(chkfile, ftextmem=self.incoordstream['cor'])
                        elif('CHARMMSTR' in filetrajformat.upper() and self.outcoordstream is not None):
                            if 'cor' in self.outcoordstream:
                                charmmcoor_dict = charmm_coor_reader(chkfile, ftextmem=self.outcoordstream['cor'])
                        else:
                            charmmcoor_dict = charmm_coor_reader(chkfile, ftextmem=None)
                        if charmmcoor_dict:
                            trajhandler_check = charmmcoor_dict
                        if trajhandler_check:
                            trajhandler = None
                            trajhandler = self.charmmcoor_iread(coorDict=trajhandler_check)
                            if trajhandler:
                                if self.interfacematch:
                                    if interface in self.topocode:
                                        usedefault=False
                                        if "input" in filetype:
                                            self.incoordhandler = trajhandler
                                            self.incoordcode = "charmmcoor"
                                            self.incoordplugin = charmmcoor_dict
                                            self.incoord_natoms = numatoms
                                        elif "output" in filetype:
                                            self.outcoordhandler = trajhandler
                                            self.outcoordcode = "charmmcoor"
                                            self.outcoordplugin = charmmcoor_dict
                                            self.outcoord_natoms = numatoms
                                        else:
                                            self.trajhandler = trajhandler
                                            self.trajcode = "charmmcoor"
                                            self.trajplugin = charmmcoor_dict
                                            self.natoms = numatoms
                                        break
                                    else:
                                        trajhandler = None
                                else:
                                    usedefault=False
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "charmmcoor"
                                        self.incoordplugin = charmmcoor_dict
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "charmmcoor"
                                        self.outcoordplugin = charmmcoor_dict
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "charmmcoor"
                                        self.trajplugin = charmmcoor_dict
                                        self.natoms = numatoms
                                    break
                if "pymolfile" in interface:
                    filetrajformat = re.sub('[.]', '', chkformat)
                    trajhandler_check = None
                    molfile_traj = None
                    if self.topohandler is not None and chkfile is not None:
                        numatoms = self.get_natoms_from_topo(self.topocode)
                        if isinstance(self.topohandler, pym.OpenMolfile):
                            molfile_traj = pym.OpenMolfile(chkfile, file_format=filetrajformat, topology=self.topohandler, silent=False)
                        elif numatoms is not None:
                            if numatoms > 0:
                                molfile_traj = pym.OpenMolfile(chkfile, file_format=filetrajformat, natoms=numatoms, silent=False)
                        if molfile_traj is not None:
                            if molfile_traj.trajectory is not None:
                                trajhandler_check = molfile_traj
                        if trajhandler_check:
                            trajhandler = None
                            trajhandler = self.pymolfile_iread(trajhandler_check)
                            if trajhandler:
                                if self.interfacematch:
                                    if interface in self.topocode:
                                        usedefault=False
                                        if "input" in filetype:
                                            self.incoordhandler = trajhandler
                                            self.incoordcode = "pymolfile"
                                            self.incoordplugin = molfile_traj
                                            self.incoord_natoms = numatoms
                                        elif "output" in filetype:
                                            self.outcoordhandler = trajhandler
                                            self.outcoordcode = "pymolfile"
                                            self.outcoordplugin = molfile_traj
                                            self.outcoord_natoms = numatoms
                                        else:
                                            self.trajhandler = trajhandler
                                            self.trajcode = "pymolfile"
                                            self.trajplugin = molfile_traj
                                            self.natoms = numatoms
                                        break
                                    else:
                                        trajhandler = None
                                else:
                                    usedefault=False
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "pymolfile"
                                        self.incoordplugin = molfile_traj
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "pymolfile"
                                        self.outcoordplugin = molfile_traj
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "pymolfile"
                                        self.trajplugin = molfile_traj
                                        self.natoms = numatoms
                                    break
                if "mdtraj" in interface:
                    numatoms = None
                    if self.topohandler is not None:
                        numatoms = self.get_natoms_from_topo(self.topocode)
                    trajhandler_check = None
                    try:
                        trajhandler_check = mdt_FormatRegistry.fileobjects[fileloadformat]
                    except KeyError:
                        pass
                    else:
                        if trajhandler_check:
                            trajhandler = None
                            trajhandler = self.mdtraj_iread(mdtraj_handler=trajhandler_check)
                            if trajhandler:
                                if self.interfacematch:
                                    if interface in self.topocode:
                                        usedefault=False
                                        if "input" in filetype:
                                            self.incoordhandler = trajhandler
                                            self.incoordcode = "mdtraj"
                                            if self.incoord_natoms is None:
                                                self.incoord_natoms = numatoms
                                        elif "output" in filetype:
                                            self.outcoordhandler = trajhandler
                                            self.outcoordcode = "mdtraj"
                                            if self.outcoord_natoms is None:
                                                self.outcoord_natoms = numatoms
                                        else:
                                            self.trajhandler = trajhandler
                                            self.trajcode = "mdtraj"
                                            if self.natoms is None:
                                                self.natoms = numatoms
                                        break
                                    else:
                                        trajhandler = None
                                else:
                                    usedefault=False
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "mdtraj"
                                        if self.incoord_natoms is None:
                                            self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "mdtraj"
                                        if self.outcoord_natoms is None:
                                            self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "mdtraj"
                                        if self.natoms is None:
                                            self.natoms = numatoms
                                    break
                if "mdanalysis" in interface:
                    mdanalysis_format = re.sub('[.]', '', fileloadformat)
                    mdanalysis_format = mdanalysis_format.upper()
                    if self.topohandler is not None and chkfile is not None:
                        # if the topology handler is a MDAnalysis universe, 
                        # we may try replacing the trajectory data in it.
                        if isinstance(self.topohandler, mda_u.Universe):
                            try:
                                trajhandler = None
                                self.topohandler.load_new(chkfile, format=mdanalysis_format)
                                trajhandler = self.mdanalysis_iread(mdanalysis_handler=self.topohandler.trajectory)
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "mdanalysis"
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "mdanalysis"
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "mdanalysis"
                                    self.natoms = numatoms
                                usedefault=False
                                break
                            except (AttributeError, IOError, OSError, ValueError, TypeError):
                                try: 
                                    universe = mda_u.Universe(self.topohandler, self.trajfile, format=mdanalysis_format)
                                    trajhandler = None
                                    if isinstance(universe, mda_u.Universe):
                                        self.topohandler = universe
                                        trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe.trajectory)
                                        if "input" in filetype:
                                            self.incoordhandler = trajhandler
                                            self.incoordcode = "mdanalysis"
                                            self.incoord_natoms = numatoms
                                        elif "output" in filetype:
                                            self.outcoordhandler = trajhandler
                                            self.outcoordcode = "mdanalysis"
                                            self.outcoord_natoms = numatoms
                                        else:
                                            self.trajhandler = trajhandler
                                            self.trajcode = "mdanalysis"
                                            self.natoms = numatoms
                                        usedefault=False
                                        break
                                    elif self.is_class_of_module(universe, 
                                            mda_c.Trajectory, mda_coordinates_modules):
                                        trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe.trajectory)
                                        if "input" in filetype:
                                            self.incoordhandler = trajhandler
                                            self.incoordcode = "mdanalysis"
                                            self.incoord_natoms = numatoms
                                        elif "output" in filetype:
                                            self.outcoordhandler = trajhandler
                                            self.outcoordcode = "mdanalysis"
                                            self.outcoord_natoms = numatoms
                                        else:
                                            self.trajhandler = trajhandler
                                            self.trajcode = "mdanalysis"
                                            self.natoms = numatoms
                                        usedefault=False
                                        break
                                #except (AttributeError, IOError, OSError, ValueError, TypeError):
                                except IOError:
                                    pass
                    # if topology handler is not a MDAnalysis Universe 
                    # or is not initialized, we can try accessing only 
                    # trajectory file. (A stub/converted topology might be supplied 
                    # for MDAnalysis to solve the needed topology data later)
                    # Hang on it: PDB has both topology and trajectory data.
                    # DuckTyper: If we can load only trajectory file, 
                    #            we may access the data.
                    try:
                        universe = mda_u.Universe(self.trajfile, format=mdanalysis_format)
                        if isinstance(universe, mda_u.Universe):
                            if self.interfacematch:
                                if interface in self.topocode:
                                    self.topohandler = universe
                                    trajhandler = None
                                    trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "mdanalysis"
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "mdanalysis"
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "mdanalysis"
                                        self.natoms = numatoms
                                    usedefault=False
                                    break
                                else:
                                    universe = None
                            else:
                                self.topohandler = universe
                                trajhandler = None
                                trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "mdanalysis"
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "mdanalysis"
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "mdanalysis"
                                    self.natoms = numatoms
                                usedefault=False
                                break
                        elif self.is_class_of_module(universe, 
                                mda_c.Trajectory, mda_coordinates_modules):
                            if self.interfacematch:
                                if interface in self.topocode:
                                    trajhandler = None
                                    trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "mdanalysis"
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "mdanalysis"
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "mdanalysis"
                                        self.natoms = numatoms
                                    usedefault=False
                                    break
                                else:
                                    universe = None
                            else:
                                trajhandler = None
                                trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "mdanalysis"
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "mdanalysis"
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "mdanalysis"
                                    self.natoms = numatoms
                                usedefault=False
                                break
                    except (AttributeError, IOError, OSError, ValueError, TypeError):
                        pass
                if "ase" in interface:
                    ase_support = None
                    ase_support = self.get_ase_format_support(fileloadformat)
                    if ase_support is None:
                        ase_support = ase_io.formats.filetype(self.trajfile)
                        trajhandler = None
                        trajhandler = self.ase_iread(self.trajfile, fileformat=ase_support)
                        try:
                            trajtest = next(trajhandler)
                        except StopIteration:
                            trajtest = None
                        trajhandler = None
                        if trajtest is not None:
                            trajhandler = self.ase_iread(self.trajfile, fileformat=ase_support)
                        if trajhandler and trajtest is not None:
                            if self.interfacematch:
                                if interface in self.topocode:
                                    usedefault=False
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "ase"
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "ase"
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "ase"
                                        self.natoms = numatoms
                                    break
                                else:
                                    trajhandler = None
                            else:
                                usedefault=False
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "ase"
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "ase"
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "ase"
                                    self.natoms = numatoms
                                break
                    else:
                        trajhandler = None
                        trajhandler = self.ase_iread(self.trajfile, fileformat=fileloadformat)
                        trajtest = None
                        try:
                            trajtest = next(trajhandler)
                        except StopIteration:
                            trajtest = None
                        trajhandler = None
                        if trajtest is not None:
                            trajhandler = self.ase_iread(self.trajfile, fileformat=ase_support)
                        if trajhandler and trajtest is not None:
                            if self.interfacematch:
                                if interface in self.topocode:
                                    usedefault=False
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "ase"
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "ase"
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "ase"
                                        self.natoms = numatoms
                                    break
                                else:
                                    trajhandler = None
                            else:
                                usedefault=False
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "ase"
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "ase"
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "ase"
                                    self.natoms = numatoms
                                break
                elif hasattr(self, "UserSuppliedInterface"):
                    if self.UserSuppliedInterface is not None:
                        if isinstance(self.UserSuppliedInterface, MDDataAccess.UserSuppliedInterface):
                            if self.UserSuppliedInterface.name in interface:
                                trajhandler = None
                                trajhandler = self.UserSuppliedInterface.trajectory_support(self.trajfile, file_format=fileloadformat)
                                if trajhandler:
                                    if self.interfacematch:
                                        if interface in self.topocode:
                                            usedefault=False
                                            self.trajhandler = trajhandler
                                            self.trajcode = self.UserSuppliedInterface.name
                                            break
                                        else:
                                            trajhandler = None
                                    else:
                                        usedefault=False
                                        self.trajhandler = trajhandler
                                        self.trajcode = self.UserSuppliedInterface.name
                                        break

        if usedefault:
            if "input" in filetype:
                trajhandler = self.incoordhandler 
            elif "output" in filetype:
                trajhandler = self.outcoordhandler 
            else:
                trajhandler = self.trajhandler 
            if trajhandler is None:
                filetrajformat = re.sub('[.]', '', chkformat)
                # First check whether pymolfile has support for the file type
                trajhandler_check = None
                if self.topohandler is not None:
                    numatoms = self.get_natoms_from_topo(topocode=self.topocode)
                    try:
                        pymHasTopo = getattr(pym.OpenMolfile, "topology")
                    except AttributeError:
                        pymHasTopo = None
                    if pymHasTopo:
                        if isinstance(self.topohandler, pym.OpenMolfile.topology):
                            topoIsPymTopo = True
                        else:
                            topoIsPymTopo = False
                    else:
                        topoIsPymTopo = False
                    molfile_traj = None
                    if chkfile is not None:
                        if(isinstance(self.topohandler, pym.OpenMolfile) or topoIsPymTopo):
                            molfile_traj = pym.OpenMolfile(chkfile, file_format=filetrajformat, topology=self.topohandler, silent=False)
                        elif numatoms is not None:
                            if numatoms > 0:
                                molfile_traj = pym.OpenMolfile(chkfile, file_format=filetrajformat, natoms=numatoms, silent=False)
                    if molfile_traj is not None:
                        if molfile_traj.trajectory is not None:
                            trajhandler_check = molfile_traj
                    if trajhandler_check:
                        trajhandler = None
                        trajhandler = self.pymolfile_iread(trajhandler_check)
                        if trajhandler:
                            if self.interfacematch:
                                if "pymolfile" in self.topocode:
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "pymolfile"
                                        self.incoordplugin = molfile_traj
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "pymolfile"
                                        self.outcoordplugin = molfile_traj
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "pymolfile"
                                        self.trajplugin = molfile_traj
                                        self.natoms = numatoms
                                else:
                                    trajhandler = None
                            else:
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "pymolfile"
                                    self.incoordplugin = molfile_traj
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "pymolfile"
                                    self.outcoordplugin = molfile_traj
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "pymolfile"
                                    self.trajplugin = molfile_traj
                                    self.natoms = numatoms
                else:
                    molfile_traj = None
                    if numatoms is not None and chkfile is not None:
                        if numatoms > 0:
                            molfile_traj = pym.OpenMolfile(chkfile, file_format=filetrajformat, natoms=numatoms, silent=False)
                    if molfile_traj is not None:
                        if molfile_traj.trajectory is not None:
                            trajhandler_check = molfile_traj
                    if trajhandler_check:
                        trajhandler = None
                        trajhandler = self.pymolfile_iread(trajhandler_check)
                        if trajhandler:
                            if self.interfacematch:
                                if interface in self.topocode:
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "pymolfile"
                                        self.incoordplugin = molfile_traj
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "pymolfile"
                                        self.outcoordplugin = molfile_traj
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "pymolfile"
                                        self.trajplugin = molfile_traj
                                        self.natoms = numatoms
                                else:
                                    trajhandler = None
                            else:
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "pymolfile"
                                    self.incoordplugin = molfile_traj
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "pymolfile"
                                    self.outcoordplugin = molfile_traj
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "pymolfile"
                                    self.trajplugin = molfile_traj
                                    self.natoms = numatoms

            if "input" in filetype:
                trajhandler = self.incoordhandler 
            elif "output" in filetype:
                trajhandler = self.outcoordhandler 
            else:
                trajhandler = self.trajhandler 
            if trajhandler is None:
                # Second,check whether MDtraj has support for the file type
                # trajhandler_check = mdt_FormatRegistry.loaders[file_format]
                try:
                    trajhandler_check = mdt_FormatRegistry.fileobjects[fileloadformat]
                except KeyError:
                    pass
                else:
                    if trajhandler_check:
                        trajhandler = None
                        trajhandler = self.mdtraj_iread(mdtraj_handler=trajhandler_check)
                        handlerresult = None
                        try:
                            handlerresult = iter(trajhandler)
                        except:
                            trajhandler = None
                        if (trajhandler is not None and 
                            handlerresult is not None):
                            # Now, we know that trajhandler works. We can ignore it at memory 
                            # and set trajhandler again in its original home at self.
                            trajhandler = None
                            if self.interfacematch:
                                if "mdtraj" in self.topocode:
                                    trajhandler = self.mdtraj_iread(mdtraj_handler=trajhandler_check)
                            else:
                                trajhandler = self.mdtraj_iread(mdtraj_handler=trajhandler_check)
                            if "input" in filetype:
                                self.incoordhandler = trajhandler
                                self.incoordcode = "mdtraj"
                                self.incoord_natoms = numatoms
                            elif "output" in filetype:
                                self.outcoordhandler = trajhandler
                                self.outcoordcode = "mdtraj"
                                self.outcoord_natoms = numatoms
                            else:
                                self.trajhandler = trajhandler
                                self.trajcode = "mdtraj"
                                self.natoms = numatoms

            # If MDTraj does not have support for the format 
            # or can not load the trajectory, use MDAnalysis or ASE.
            if "input" in filetype:
                trajhandler = self.incoordhandler 
            elif "output" in filetype:
                trajhandler = self.outcoordhandler 
            else:
                trajhandler = self.trajhandler 
            if trajhandler is None:
                mdanalysis_format = re.sub('[.]', '', fileloadformat)
                if self.topohandler is not None and chkfile is not None:
                    if isinstance(self.topohandler, mda_u.Universe):
                        try:
                            self.topohandler.load_new(chkfile, format=mdanalysis_format)
                            trajhandler = self.mdanalysis_iread(mdanalysis_handler=self.topohandler)
                            if "input" in filetype:
                                self.incoordhandler = trajhandler
                                self.incoordcode = "mdanalysis"
                                self.incoord_natoms = numatoms
                            elif "output" in filetype:
                                self.outcoordhandler = trajhandler
                                self.outcoordcode = "mdanalysis"
                                self.outcoord_natoms = numatoms
                            else:
                                self.trajhandler = trajhandler
                                self.trajcode = "mdanalysis"
                                self.natoms = numatoms
                        except (AttributeError, IOError, OSError, ValueError, TypeError):
                            try: 
                                universe = mda_u.Universe(self.topohandler, self.trajfile, format=mdanalysis_format)
                                if isinstance(universe, mda_u.Universe):
                                    self.topohandler = universe
                                    trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "mdanalysis"
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "mdanalysis"
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "mdanalysis"
                                        self.natoms = numatoms
                                elif self.is_class_of_module(self.universe, 
                                        mda_c.Trajectory, mda_coordinates_modules):
                                    trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                    if "input" in filetype:
                                        self.incoordhandler = trajhandler
                                        self.incoordcode = "mdanalysis"
                                        self.incoord_natoms = numatoms
                                    elif "output" in filetype:
                                        self.outcoordhandler = trajhandler
                                        self.outcoordcode = "mdanalysis"
                                        self.outcoord_natoms = numatoms
                                    else:
                                        self.trajhandler = trajhandler
                                        self.trajcode = "mdanalysis"
                                        self.natoms = numatoms
                            except (AttributeError, IOError, OSError, ValueError, TypeError):
                                pass
                try:
                    universe = mda_u.Universe(self.trajfile, format=mdanalysis_format)
                    if isinstance(universe, mda_u.Universe):
                        trajhandler = None
                        if self.interfacematch:
                            if "mdanalysis" in self.topocode:
                                self.topohandler = universe
                                trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                if "input" in filetype:
                                    self.incoordhandler = trajhandler
                                    self.incoordcode = "mdanalysis"
                                    self.incoord_natoms = numatoms
                                elif "output" in filetype:
                                    self.outcoordhandler = trajhandler
                                    self.outcoordcode = "mdanalysis"
                                    self.outcoord_natoms = numatoms
                                else:
                                    self.trajhandler = trajhandler
                                    self.trajcode = "mdanalysis"
                                    self.natoms = numatoms
                            else:
                                universe = None
                        else:
                            self.topohandler = universe
                            trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                            if "input" in filetype:
                                self.incoordhandler = trajhandler
                                self.incoordcode = "mdanalysis"
                                self.incoord_natoms = numatoms
                            elif "output" in filetype:
                                self.outcoordhandler = trajhandler
                                self.outcoordcode = "mdanalysis"
                                self.outcoord_natoms = numatoms
                            else:
                                self.trajhandler = trajhandler
                                self.trajcode = "mdanalysis"
                                self.natoms = numatoms
                    elif self.is_class_of_module(self.universe, 
                            mda_c.Trajectory, mda_coordinates_modules):
                        trajhandler = None
                        trajcode = None
                        if self.interfacematch:
                            if "mdanalysis" in self.topocode:
                                trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                                trajcode = "mdanalysis"
                            else:
                                universe = None
                        else:
                            trajhandler = self.mdanalysis_iread(mdanalysis_handler=universe)
                            trajcode = "mdanalysis"
                        if "input" in filetype:
                            self.incoordhandler = trajhandler
                            self.incoordcode = trajcode
                            self.incoord_natoms = numatoms
                        elif "output" in filetype:
                            self.outcoordhandler = trajhandler
                            self.outcoordcode = trajcode
                            self.outcoord_natoms = numatoms
                        else:
                            self.trajhandler = trajhandler
                            self.trajcode = trajcode
                            self.natoms = numatoms
                except (AttributeError, IOError, OSError, ValueError, TypeError):
                    pass

            if "input" in filetype:
                trajhandler = self.incoordhandler 
            elif "output" in filetype:
                trajhandler = self.outcoordhandler 
            else:
                trajhandler = self.trajhandler 
            if trajhandler is None:
                ase_support = None
                ase_support = self.get_ase_format_support(fileloadformat)
                trajcode=None
                # May still have chance that ASE can recognize the 
                # format with its filetype checking function
                if ase_support is None:
                    if chkfile is not None:
                        try:
                            ase_support = ase_io.formats.filetype(chkfile)
                        except (FileNotFoundError,IOError):
                            pass
                    trajhandler = None
                    trajtest = None
                    if chkfile is not None:
                        trajhandler = self.ase_iread(chkfile, fileformat=ase_support)
                        try:
                            trajtest = next(trajhandler)
                        except StopIteration:
                            trajtest = None
                        trajhandler = None
                        if trajtest is not None:
                            trajhandler = self.ase_iread(self.trajfile, fileformat=ase_support)
                    if trajhandler and trajtest is not None:
                        trajcode=None
                        if self.interfacematch:
                            if "ase" in self.topocode:
                                trajcode = "ase"
                            else:
                                trajhandler = None
                        else:
                            trajcode = "ase"
                else:
                    trajhandler = None
                    trajtest = None
                    if chkfile is not None:
                        trajhandler = self.ase_iread(chkfile, fileformat=fileloadformat)
                        try:
                            trajtest = next(trajhandler)
                        except StopIteration:
                            trajtest = None
                        trajhandler = None
                        if trajtest is not None:
                            trajhandler = self.ase_iread(self.trajfile, fileformat=ase_support)
                    trajcode=None
                    if trajhandler and trajtest is not None:
                        if self.interfacematch:
                            if "ase" in self.topocode:
                                trajcode = "ase"
                            else:
                                trajhandler = None
                        else:
                            trajcode = "ase"
                if "input" in filetype:
                    self.incoordhandler = trajhandler
                    self.incoordcode = trajcode
                    self.incoord_natoms = numatoms
                elif "output" in filetype:
                    self.outcoordhandler = trajhandler
                    self.outcoordcode = trajcode
                    self.outcoord_natoms = numatoms
                else:
                    self.trajhandler = trajhandler
                    self.trajcode = trajcode
                    self.natoms = numatoms

        if "input" in filetype:
            trajhandler = self.incoordhandler 
        elif "output" in filetype:
            trajhandler = self.outcoordhandler 
        else:
            trajhandler = self.trajhandler 
        if trajhandler is None:
            return False
        else:
            return True

    def iread(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        iterator_object = iter(self.trajhandler)
        try:
            while True:
                self.trajiter = next(iterator_object)
                self.atompositions = self.trajiter
                return self.trajiter
        except StopIteration:
            pass
        finally:
            del iterator_object

    def incoord_iread(self):
        """Returns an iterator that goes through the given file one
        configuration at a time.
        """
        try:
            iterator_object = iter(self.incoordhandler)
            try:
                while True:
                    try:
                        self.incoorditer = next(iterator_object)
                        return self.incoorditer
                    except ValueError:
                        pass
            except StopIteration:
                pass
            finally:
                del iterator_object
        except TypeError:
            pass

    def outcoord_iread(self):
        """Returns an iterator that goes through the given file one
        configuration at a time.
        """
        try:
            iterator_object = iter(self.outcoordhandler)
            try:
                while True:
                    try:
                        self.outcoorditer = next(iterator_object)
                        return self.outcoorditer
                    except ValueError:
                        pass
            except StopIteration:
                pass
            finally:
                del iterator_object
        except TypeError:
            pass

    def topology_iread(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        try:
            iterator_object = iter(self.topohandler)
            try:
                while True:
                    try:
                        self.trajiter = next(iterator_object)
                        return self.trajiter
                    except ValueError:
                        pass
            except StopIteration:
                pass
            finally:
                del iterator_object
        except TypeError:
            pass

    def get_topology(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        return self.topohandler

    def get_unitcell(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        return self.trajhandler.cell_lengths
    
    def load_pymolfile_topology(self, file_format=None):
        """ function to call pymolfile topology reader

        Returns
        -------
        topology : pymolfile Topology in OpenMolfile Class
        """

        topology = None
        universe = None
        if file_format:
            fileformat = re.sub('[.]', '', file_format)
        else:
            fileformat = file_format
        if self.topofile:
            try:
                molfile_topo = pym.OpenMolfile(self.topofile, file_format=fileformat)
                if molfile_topo.topology is not None:
                    topology = molfile_topo
                else:
                    if molfile_topo.trajectory is None:
                        molfile_topo = None
                        del molfile_topo
            except (AttributeError, IOError, OSError, ValueError):
                pass

        return topology

    def load_gromos_topology(self, file_format, base_topo=None):
        """
        Get the topology and parameters from file loaders in parmed

        Returns
        -------
        topology : GromosTopoObjects
        """

        if base_topo is not None:
            topology = base_topo
        else:
            topology = None

        ext = file_format
        top = self.topofile

        if ext in ['gromostop', 'gromostop.gz', 'gromostop.zip',
                   'top', 'top.gz', 'top.zip']:
            if base_topo is None:
                base_topo = gto.GromosTopology()
            try:
                if self.topofile:
                    if isinstance(base_topo, gto.GromosTopology):
                        base_topo.load_top(top, 'gromostop')
                topology = base_topo
            except(ValueError, AttributeError, IOError):
                pass
        elif ext in ['gromoscnf', 'gromoscnf.gz', 'gromoscnf.zip',
                     'cnf', 'cnf.gz', 'cnf.zip']:
            if base_topo is None:
                base_topo = gto.GromosTopology()
            try:
                if self.topofile:
                    if isinstance(base_topo, gto.GromosTopology):
                        base_topo.load_cnf(top, 'gromoscnf')
                topology = base_topo
            except(ValueError, AttributeError, IOError):
                pass

        return topology

    def load_parmed_topology(self, file_format, base_topo=None):
        """
        Get the topology and parameters from file loaders in parmed

        Returns
        -------
        topology : parmed.topologyobjects
        """

        if base_topo is not None:
            topology = base_topo
        else:
            topology = None

        ext = file_format
        top = self.topofile

        if ext in ['mol2', 'mol2.gz', 'mol2.bz2']:
            if self.topofile:
                try:
                    topology = pmd.load_file(top, structure=True)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['psf', 'psf.gz', 'psf.bz2',
                     'pdb', 'pdb.gz', 'pdb.bz2',
                     'cif', 'cif.gz', 'cif.bz2',
                     'sdf', 'sdf.gz', 'sdf.bz2']:
            if self.topofile:
                parmset_topo = None
                if base_topo is not None:
                    if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                        parmset_topo = base_topo
                try:
                    topology = pmd.load_file(top)
                    if parmset_topo is not None:
                        topology.load_parameters(parmset_topo, copy_parameters=True)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['mdl', 'mdl.gz', 'mdl.bz2']:
            if self.topofile:
                try:
                    topology = pmd.amber.AmberFormat(top)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['prmtop', 'prmtop.gz', 'prmtop.bz2',
                     'parm7', 'parm7.gz', 'parm7.bz2']:
            if self.topofile:
                try:
                    topology = pmd.amber.LoadParm(top)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['top', 'gromacstop']:
            if self.topofile:
                try:
                    topology = pmd.gromacs.GromacsTopologyFile(top)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['tinkertop', 'tinkerxyz', 'txyz']:
            if self.topofile:
                try:
                    topology = pmd.tinker.tinkerfiles.XyzFile(top)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['tinkercor', 'tinkerdyn', 'dyn']:
            if self.topofile:
                try:
                    topology = pmd.tinker.tinkerfiles.DynFile(top)
                except(ValueError, AttributeError, IOError):
                    pass
        elif ext in ['rtf', 'charmmtop', 'charmmstrrtf']: # .top is used by Gromacs.
            if base_topo is None:
                base_topo = pmd.charmm.CharmmParameterSet()
            try:
                # top can be str or list of lines
                if ext in ['charmmstrrtf'] and self.topostream is not None:
                    if 'rtf' in self.topostream:
                        if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                            base_topo.read_topology_file(iter(self.topostream['rtf']))
                        elif isinstance(base_topo, pmd.charmm.CharmmPsfFile):
                            self.topocharmmpsf = True
                            parmset_topo = pmd.charmm.CharmmParameterSet()
                            parmset_topo.read_topology_file(iter(self.topostream['rtf']))
                            base_topo.load_parameters(parmset_topo, copy_parameters=True)
                else:
                    if self.topofile:
                        if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                            base_topo.read_topology_file(top)
                        elif isinstance(base_topo, pmd.charmm.CharmmPsfFile):
                            self.topocharmmpsf = True
                            parmset_topo = pmd.charmm.CharmmParameterSet()
                            parmset_topo.read_topology_file(top)
                            base_topo.load_parameters(parmset_topo, copy_parameters=True)
                topology = base_topo
            except(ValueError, AttributeError, IOError):
                pass
        elif ext in ['stream', 'str', 'charmmstream']:
            if base_topo is None:
                base_topo = pmd.charmm.CharmmParameterSet()
            try:
                # top can be str or list of lines
                if ext in ['charmmstream'] and self.topostream is not None:
                    if isinstance(self.topostream,dict):
                        if 'rtf' in self.topostream:
                            if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                                base_topo.read_topology_file(iter(self.topostream['rtf']))
                            elif isinstance(base_topo, pmd.charmm.CharmmPsfFile):
                                self.topocharmmpsf = True
                                parmset_topo = pmd.charmm.CharmmParameterSet()
                                parmset_topo.read_topology_file(iter(self.topostream['rtf']))
                                base_topo.load_parameters(parmset_topo, copy_parameters=True)
                        if 'par' in self.topostream:
                            if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                                base_topo.read_parameter_file(iter(self.topostream['par']))
                            elif isinstance(base_topo, pmd.charmm.CharmmPsfFile):
                                self.topocharmmpsf = True
                                parmset_topo = pmd.charmm.CharmmParameterSet()
                                parmset_topo.read_parameter_file(iter(self.topostream['par']))
                                base_topo.load_parameters(parmset_topo, copy_parameters=True)
                else:
                    if self.topofile:
                        base_topo.read_stream_file(top)
                topology = base_topo
            except(ValueError, AttributeError, IOError):
                pass
        elif ext in ['par', 'prm', 'charmmpar', 'charmmstrpar']:
            if base_topo is None:
                base_topo = pmd.charmm.CharmmParameterSet()
            try:
                # top can be str or list of lines
                if ext in ['charmmstrpar'] and self.topostream is not None:
                    if 'par' in self.topostream:
                        if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                            base_topo.read_parameter_file(iter(self.topostream['par']))
                        elif isinstance(base_topo, pmd.charmm.CharmmPsfFile):
                            self.topocharmmpsf = True
                            parmset_topo = pmd.charmm.CharmmParameterSet()
                            parmset_topo.read_parameter_file(iter(self.topostream['par']))
                            base_topo.load_parameters(parmset_topo, copy_parameters=True)
                else:
                    if self.topofile:
                        if isinstance(base_topo, pmd.charmm.CharmmParameterSet):
                            base_topo.read_parameter_file(top)
                        elif isinstance(base_topo, pmd.charmm.CharmmPsfFile):
                            self.topocharmmpsf = True
                            parmset_topo = pmd.charmm.CharmmParameterSet()
                            parmset_topo.read_parameter_file(top)
                            base_topo.load_parameters(parmset_topo, copy_parameters=True)
                topology = base_topo
            except(ValueError, AttributeError, IOError):
                pass

        return topology

    def load_mdtraj_topology(self, file_format, **kwargs):
        """Borrowed from mdtraj.trajectory to remove 
           IOError and TypeError warnings
    
        Get the topology from an argument of indeterminate type.
        If top is a string, we try loading a pdb, if its a trajectory
        we extract its topology.

        Returns
        -------
        topology : mdtraj.Topology
        """

        topology = None
        wrapkwargs = kwargs.copy()
        wrapkwargs.pop("top", None)
        wrapkwargs.pop("atom_indices", None)
        wrapkwargs.pop("frame", None)

        ext=file_format
        top=self.topofile

        if top:
            if ext in ['.pdb', '.pdb.gz', '.h5','.lh5']:
                _traj = mdt.load_frame(top, 0, **wrapkwargs)
                topology = _traj.topology
            elif ext in ['.prmtop', '.parm7']:
                topology = mdt_load.prmtop(top, **wrapkwargs)
            elif ext in ['.psf']:
                topology = mdt_load.psf(top, **wrapkwargs)
            elif ext in ['.mol2']:
                topology = mdt_load.mol2(top, **wrapkwargs).topology
            elif ext in ['.gro']:
                topology = mdt.core.trajectory.load_gro(top, **wrapkwargs).topology
            elif ext in ['.arc']:
                #topology = mdt_load.arc(top, **wrapkwargs).topology
                topology = mdt.core.trajectory.load_arc(top, **wrapkwargs).topology
            elif ext in ['.hoomdxml']:
                topology = mdt_load.hoomdxml(top, **wrapkwargs).topology
            elif isinstance(top, mdt_Trajectory):
                topology = top.topology
            elif isinstance(top, mdt_Topology):
                topology = top

        return topology

    def load_mdanalysis_topology(self, file_format=None):
        """ Support level for MDAnalysis package

        If filename is given, we try loading the format, if its a trajectory
        we extract its topology.

        Returns
        -------
        topology : MDAnalysis.core.topology.Topology
        """

        topology = None
        universe = None
        if file_format:
            fileformat = re.sub('[.]', '', file_format)
        else:
            fileformat = file_format
        top=self.topofile
        if top:
            try:
                universe = mda_u.Universe(top, topology_format=fileformat)
                if isinstance(universe, mda_u.Universe):
                    topology = universe
                elif isinstance(universe, mda_u.Topology):
                    topology = universe
                elif isinstance(universe._topology, mda_u.Topology):
                    topology = universe._topology
            except (AttributeError, IOError, OSError, ValueError):
                pass

        return topology

    def load_ase_support(self, file_name, file_format=None):
        """Get the topology from a argument of indeterminate type
        If top is a string, we try loading a pdb, if its a trajectory
        we extract its topology.

        Returns
        -------
        topology : generator
        """
        ioformat = None
        ioformat = self.get_ase_format_support(file_format)
        if ioformat is not None:
            handler = None
            try:
                handler = ase_io.read(file_name, format=ioformat)
            except (ModuleNotFoundError, ImportError, AttributeError, 
                IOError, OSError, ValueError):
                pass
            if handler is None:
                try:
                    handler_iread = ase_io.iread(file_name, index=":", format=ioformat)
                    try:
                        handler_iter = iter(handler_iread)
                        handler_iread = None
                        handler = ase_io.iread(file_name, index=":", format=ioformat)
                    except (ModuleNotFoundError, ImportError, AttributeError, 
                        IOError, OSError, ValueError):
                        pass
                except (ModuleNotFoundError, ImportError, AttributeError, 
                    IOError, OSError, ValueError):
                    pass
            return handler
        return None

    def get_ase_format_support(self, file_format):
        """ Inspired from ase.io.formats.get_ioformat

        Get support level of ASE for file formats

        Returns
        -------
        IOformat list
        """

        ioformat = None
        module_handler = None
        if file_format:
            _format = file_format.replace('-', '_')
            module_name = ase_io.formats.format2modulename.get(file_format, _format)
            try:
                module_handler = ase_io.formats.import_module('ase.io.' + module_name)
            except ImportError:
                pass
        if module_handler:
            module_read = getattr(module_handler, 'read_' + _format, None)
            if module_read and not inspect.isgeneratorfunction(module_read):
                read = functools.partial(wrap_read_function, module_read)
                if module_read:
                    ioformat = all_formats[file_format]

        return ioformat

    def ase_iread(self, filename, fileformat):
        """Used to wrap an iterator returned by ase.io.iread so that it returns
        the positions instead of the ase.Atoms object.
        """
        # The newest ASE version found in Github has an iread function.
        # After reading the ASE source code, it seems that the ASE iread does
        # actually read the entire file into memory and the yields the
        # configurations from it. Should be checked at some point.
        handler = None
        try:
            handler = ase_io.iread(filename, index=":", format=fileformat)
        except (AttributeError, IOError, OSError, 
                ValueError, ImportError):
            return handler

        if handler is not None:
            try:
                for value in handler:
                    yield value.get_positions()
            except (AttributeError, IOError, OSError, 
                    ValueError, ImportError):
                return handler
        else:
            return handler
    
    def charmm_coor_toporead(self, coorDict=None):
        """Returns the atom list for CHARMM COOR input
        """
        if coorDict is not None:
            if 'atomlist' in coorDict:
                return coorDict['atomlist']
            else:
                return None
        else:
            return None

    def charmmcoor_iread(self, coorDict=None):
        """Returns the iterator for CHARMM COOR trajectory
        """
        if "input" in self.trajtype:
            traj = self.inputcoords
        elif "output" in self.trajtype:
            traj = self.outputcoords
        else:
            traj = self.trajectory

        if coorDict is not None:
            if traj is None:
                traj = MDDataTrajectory()
            if 'ww' in coorDict:
                traj.charmm = {'ww':coorDict['ww']}
            if "input" in self.trajtype:
                self.inputcoords = traj
            elif "output" in self.trajtype:
                self.outputcoords = traj
            else:
                self.trajectory = traj
            if 'coords' in coorDict:
                for positions in coorDict['coords']:
                    if len(positions)>0:
                        yield np.asarray(positions)
                    else:
                        return None
        return None

    def gromos_iread(self, traj_handler=None):
        """Returns the iterator for gromos trajectory
        """
        if "input" in self.trajtype:
            traj = self.inputcoords
        elif "output" in self.trajtype:
            traj = self.outputcoords
        else:
            traj = self.trajectory

        if traj_handler is not None:
            while True:
                positions = traj_handler.iread()
                if positions is not None:
                    if traj is None:
                        traj = MDDataTrajectory()
                    if traj_handler.velocities is not None:
                        traj.velocities = traj_handler.velocities
                    if traj_handler.unitcell is not None:
                        traj.unitcell_lengths = np.asarray(
                            traj_handler.unitcell[0:3])
                        traj.unitcell_angles = np.asarray(
                            traj_handler.unitcell[3:6])
                        a = traj_handler.unitcell[0]
                        b = traj_handler.unitcell[1]
                        c = traj_handler.unitcell[2]
                        alpha = traj_handler.unitcell[3]
                        beta = traj_handler.unitcell[4]
                        gamma = traj_handler.unitcell[5]
                        lx = a
                        xy = (b * math.cos(gamma*math.pi/180.0))
                        xy = 0.0 if xy<1e-10 else xy
                        xz = (c * math.cos(beta*math.pi/180.0))
                        xz = 0.0 if xz<1e-10 else xz
                        bxy = b*b-xy*xy
                        if bxy>0:
                            ly = math.sqrt(bxy)
                        else:
                            ly = b if b>0 else 1.0
                        yz = (b*c*math.cos(alpha*math.pi/180.0)-xy*xz)/ly
                        yz = 0.0 if yz<1e-10 else yz
                        cxz = c*c-xz*xz-yz*yz
                        if cxz>0:
                            lz = math.sqrt(cxz)
                        else:
                            lz = c if c>0 else 1.0
                        traj.unitcell_vectors = np.zeros((3,3))
                        traj.unitcell_vectors[0][0] = lx if lx != 0 else 1.0
                        traj.unitcell_vectors[0][1] = xy
                        traj.unitcell_vectors[0][2] = xz
                        traj.unitcell_vectors[1][1] = ly if ly != 0 else 1.0
                        traj.unitcell_vectors[1][2] = yz
                        traj.unitcell_vectors[2][2] = lz
                    if "input" in self.trajtype:
                        self.inputcoords = traj
                    elif "output" in self.trajtype:
                        self.outputcoords = traj
                    else:
                        self.trajectory = traj
                    yield traj_handler.positions
                else:
                    return None

    def pymolfile_iread(self, traj_handler=None):
        """Returns the iterator for pymolfile trajectory
        """
        if "input" in self.trajtype:
            traj = self.inputcoords
        elif "output" in self.trajtype:
            traj = self.outputcoords
        else:
            traj = self.trajectory

        if traj_handler.trajectory is not None:
            while True:
                positions = traj_handler.trajectory.iread()
                if positions is not None:
                    if traj is None:
                        traj = MDDataTrajectory()
                    if "has_velocities" in positions:
                        if positions["has_velocities"]>0:
                            traj.velocities = positions["velocities"]
                    if("A" in positions and 
                       "B" in positions and 
                       "C" in positions):
                        traj.unitcell_lengths = np.asarray(
                                [positions["A"],positions["B"],positions["C"]])
                    if("alpha" in positions and
                       "beta" in positions and 
                       "gamma" in positions):
                        traj.unitcell_angles = np.asarray(
                                [positions["alpha"],positions["beta"],positions["gamma"]])
                    if("A" in positions and "B" in positions and 
                       "C" in positions and "alpha" in positions and
                       "beta" in positions and "gamma" in positions):
                        a = positions["A"]
                        b = positions["B"]
                        c = positions["C"]
                        alpha = positions["alpha"]
                        beta = positions["beta"]
                        gamma = positions["gamma"]
                        lx = a
                        xy = (b * math.cos(gamma*math.pi/180.0))
                        xy = 0.0 if xy<1e-10 else xy
                        xz = (c * math.cos(beta*math.pi/180.0))
                        xz = 0.0 if xz<1e-10 else xz
                        bxy = b*b-xy*xy
                        if bxy>0:
                            ly = math.sqrt(bxy)
                        else:
                            ly = b if b>0 else 1.0
                        yz = (b*c*math.cos(alpha*math.pi/180.0)-xy*xz)/ly
                        yz = 0.0 if yz<1e-10 else yz
                        cxz = c*c-xz*xz-yz*yz
                        if cxz>0:
                            lz = math.sqrt(cxz)
                        else:
                            lz = c if c>0 else 1.0
                        traj.unitcell_vectors = np.zeros((3,3))
                        traj.unitcell_vectors[0][0] = lx if lx != 0 else 1.0
                        traj.unitcell_vectors[0][1] = xy
                        traj.unitcell_vectors[0][2] = xz
                        traj.unitcell_vectors[1][1] = ly if ly != 0 else 1.0
                        traj.unitcell_vectors[1][2] = yz
                        traj.unitcell_vectors[2][2] = lz
                    if "input" in self.trajtype:
                        self.inputcoords = traj
                    elif "output" in self.trajtype:
                        self.outputcoords = traj
                    else:
                        self.trajectory = traj
                    yield positions["coords"]
                else:
                    return None

    def mdanalysis_iread(self, mdanalysis_handler=None):
        """Returns the iterator for mdanalysis trajectory
        """
        if mdanalysis_handler is not None:
            for tstep in mdanalysis_handler:
                if tstep.has_positions:
                    data = tstep.positions
                    if isinstance(data, tuple):
                        positions = data[0]
                    else:
                        positions = data
                    if len(positions) == 0:
                        break
                    else:
                        yield positions

    def mdtraj_iread(self, mdtraj_handler=None, **kwargs):
        """Generator function that is used to read an atomic configuration file (MD
        trajectory, geometry optimization, static snapshot) from a file one frame
        at a time. Only the xyz positions are returned from the file, and no unit
        conversion is done, so you have to be careful with units.

        By using a generator pattern we can avoid loading the entire trajectory
        file into memory. This function will instead load a chunk of the file into
        memory (with MDTraj you can decide the chunk size, with ASE it seems to
        always be one frame), and serve individual files from that chunk. Once the
        frames in one chunk are iterated, the chunk will be garbage collected and
        memory is freed.

        Args:
            filename: String for the file path.
            file_format: String for the file format. If not given the format is
                automatically detected from the extension.

        Yields:
            numpy array containing the atomic positions in one frame.
        """

        # Try to open the file with MDTraj first. With a brief inspection it seems
        # that MDTraj is better performance wise, because it can iteratively load a
        # "chunk" of frames, and still serve the individual frames one by one. ASE
        # on the other hand will iteratively read frames one by one (unnecessary
        # IO).

        # Must use the low level MDTraj API to open files without topology.
        # mdtraj_supported_format = FormatRegistry.loaders[self.fileformat]
        if "input" in self.trajtype:
            traj = self.inputcoords
            trajformat = self.incoordformat
            trajfile = self.incoordfile
            chunk = 1
        elif "output" in self.trajtype:
            traj = self.outputcoords
            trajformat = self.outcoordformat
            trajfile = self.outcoordfile
            chunk = 1
        else:
            traj = self.trajectory
            trajformat = self.trajformat
            trajfile = self.trajfile
            chunk = self.trajchunk

        if mdtraj_handler is not None:
            if chunk == 0:
                try:
                    loader = mdt_FormatRegistry.loaders[trajformat]
                except KeyError:
                    return
                # If chunk was 0 then we want to avoid filetype-specific code
                # in case of undefined behavior in various file parsers.
                # TODO: this will first apply stride, then skip!
                if trajformat not in mdt_TOPOLOGY_EXTS:
                    topkwargs = kwargs.copy()
                    topkwargs['top'] = self.topohandler
                    # standard_names is a valid keyword argument only for files containing topologies
                    topkwargs.pop('standard_names', None)
                    positions = loader(trajfile, **topkwargs)
                else:
                    positions = loader(trajfile)
                for pos in positions:
                    yield pos
            elif trajformat in ('.pdb', '.pdb.gz'):
                # the PDBTrajectoryFile class doesn't follow the standard API. Fixing it here
                try:
                    loader = mdt_FormatRegistry.loaders[trajformat]
                except KeyError:
                    return
                t = loader(trajfile)
                for i in range(0, len(t), chunk):
                    positions = t[i:i+chunk]
                    for pos in positions:
                        yield pos
            else:
                if self.topohandler is None:
                    n_atoms_set=None
                    if trajformat in ('.crd', '.mdcrd'):
                        return
                else:
                    n_atoms_set = None
                    n_atoms_set = self.get_natoms_from_topo(self.topocode)
                    self.mddata_topology = None
                    if self.topology is not None:
                        self.mddata_topology = self.topology
                        if isinstance(self.topohandler, mdt_Topology):
                            self.topology = self.topohandler
                        else:
                            self.topology = mdt_Topology()
                    else:
                        self.topology = mdt_Topology()
                    self.topology._numAtoms = n_atoms_set
                try:
                    with (lambda x: mdtraj_handler(x, n_atoms=n_atoms_set)
                          if trajformat in ('.crd', '.mdcrd')
                          else mdtraj_handler(trajfile, mode="r"))(trajfile) as f:
                        empty = False
                        while not empty:
                            if trajformat not in mdt_TOPOLOGY_EXTS:
                                data = f.read_as_traj(self.topology, n_frames=chunk)
                            else:
                                data = f.read_as_traj(n_frames=chunk)
                            for pos in data.xyz:
                                if self.trajectory is not None:
                                    self.trajectory.unitcell_vectors = data.unitcell_vectors
                                    self.trajectory.unitcell_lengths = data.unitcell_lengths
                                    self.trajectory.unitcell_angles = data.unitcell_angles
                                    self.trajectory.time = data.time
                                yield pos
                            empty = True
                    if self.mddata_topology is not None:
                        self.topology = self.mddata_topology
                except IOError:
                    self.topology = self.mddata_topology
                    return

    def custom_iread(self, file_handle, format):
        """
        """
        pass

    def topologyFileHandler(self, filename, fileformatlist=None, interfacelist=None):
        self.init_topo()
        self.set_topo(filename)
        if (self.interfaceorder is None and 
            interfacelist is not None):
            self.interfaceorder = interfacelist
        self.load_topology()
        if self.topohandler is None:
            for fileformat in fileformatlist:
                self.init_topo()
                self.interfaceorder = interfacelist
                self.set_topo(filename, fileformat)
                self.load_topology()
                if self.topohandler is not None:
                    return fileformat

    def trajectoryFileHandler(self, filename, fileformatlist=None, interfacelist=None):
        self.init_traj()
        self.set_traj(filename)
        self.interfaceorder = interfacelist
        traj_loaded = self.load()
        self.atompositions = self.iread()
        if self.atompositions is None:
            for fileformat in fileformatlist:
                self.init_traj()
                self.set_traj(filename, fileformat)
                traj_loaded = self.load()
                self.atompositions = self.iread()
                if self.atompositions is not None:
                    return fileFormat

    def incoordFileHandler(self, filename, fileformatlist=None, interfacelist=None):
        self.init_incoord()
        self.set_incoord(filename)
        self.interfaceorder = interfacelist
        incoord_loaded = self.load_incoord()
        self.inputpositions = self.incoord_iread()
        if self.inputpositions is None:
            for fileformat in fileformatlist:
                self.init_incoord()
                self.set_incoord(filename, fileformat)
                incoord_loaded = self.load_incoord()
                self.inputpositions = self.incoord_iread()
                if self.inputpositions is not None:
                    return fileFormat

    def outcoordFileHandler(self, filename, fileformatlist=None, interfacelist=None):
        self.init_outcoord()
        self.set_outcoord(filename)
        self.interfaceorder = interfacelist
        outcoord_loaded = self.load_outcoord()
        self.outputpositions = self.outcoord_iread()
        if self.outputpositions is None:
            for fileformat in fileformatlist:
                self.init_outcoord()
                self.set_outcoord(filename, fileformat)
                outcoord_loaded = self.load_outcoord()
                self.outputpositions = self.outcoord_iread()
                if self.outputpositions is not None:
                    return fileFormat

    def thermoFileHandler(self, filename, fileformatlist=None, interfacelist=None):
        self.init_thermo()
        self.set_thermo(filename)
        if (self.interfaceorder is None and 
            interfacelist is not None):
            self.interfaceorder = interfacelist
        self.load_thermo()
        if self.thermohandler is None:
            for fileformat in fileformatlist:
                self.init_thermo()
                self.interfaceorder = interfacelist
                self.set_thermo(filename, fileformat)
                self.load_thermo()
                if self.thermohandler is not None:
                    return fileformat

    def initializeFileHandlers(self, parser_ui):
        # Files will be loaded using their extensions initially.
        # If this fails, the fileFormat lists will be used in loading process.
        topoformat = None
        trajformat = None
        self.atompositions = None
        self.topohandler = None
        self.thermohandler = None
        self.forcefieldhandler = None
        for fileItem in parser_ui.fileDict:
            fileformatlist = None
            interfacelist = None
            if (parser_ui.fileDict[fileItem].fileSupplied and
                parser_ui.fileDict[fileItem].activeInfo):
                filename = parser_ui.fileDict[fileItem].fileName
                fileformatlist = parser_ui.fileDict[fileItem].fileFormat
                interfacelist = parser_ui.fileDict[fileItem].fileInterface
                if 'topology' in parser_ui.fileDict[fileItem].infoPurpose:
                    topoformat = self.topologyFileHandler(filename, fileformatlist, interfacelist)
                #if 'forcefield' in parser_ui.fileDict[fileItem].infoPurpose:
                #    forcefieldformat = self.forcefieldFileHandler(filename, fileformatlist, interfacelist)
                if 'inputcoordinates' in parser_ui.fileDict[fileItem].infoPurpose:
                    incoordformat = self.incoordFileHandler(filename, fileformatlist, interfacelist)
                if 'outputcoordinates' in parser_ui.fileDict[fileItem].infoPurpose:
                    outcoordformat = self.outcoordFileHandler(filename, fileformatlist, interfacelist)
                if 'trajectory' in parser_ui.fileDict[fileItem].infoPurpose:
                    trajformat = self.trajectoryFileHandler(filename, fileformatlist, interfacelist)
                if 'thermostats' in parser_ui.fileDict[fileItem].infoPurpose:
                    thermoformat = self.thermoFileHandler(filename, fileformatlist, interfacelist)
        if self.topohandler is not None and self.topocode is not None:
            parser_ui.topology = self.set_TopologyData()
        if self.inputpositions is not None:
            parser_ui.inputcoords = self.set_InputData()
        if self.outputpositions is not None:
            parser_ui.outputcoords = self.set_OutputData()
        if self.atompositions is not None:
            parser_ui.trajectory = self.set_TrajectoryData()
        if self.thermohandler is not None:
            parser_ui.thermostats = self.set_ThermoData()
            parser_ui.thermoDict = parser_ui.thermostats.thermoDict
        #if self.forcefieldhandler is not None:
        #    #self.set_ForceFieldData()
        #    pass

    def initializeTopologyFileHandlers(self, parser_ui):
        # Files will be loaded using their extensions initially.
        # If this fails, the fileFormat lists will be used in loading process.
        self.init_topo()
        topoformat = None
        self.topohandler = None
        for fileItem in parser_ui.fileDict:
            fileformatlist = None
            interfacelist = None
            if (parser_ui.fileDict[fileItem].fileSupplied and
                parser_ui.fileDict[fileItem].activeInfo):
                filename = parser_ui.fileDict[fileItem].fileName
                fileformatlist = parser_ui.fileDict[fileItem].fileFormat
                interfacelist = parser_ui.fileDict[fileItem].fileInterface
                if 'topology' in parser_ui.fileDict[fileItem].infoPurpose:
                    if interfacelist is not None:
                        self.interfaceorder = interfacelist
                    for fileformat in fileformatlist:
                        self.set_topo(filename, fileformat)
                        if 'strDict' in parser_ui.fileDict[fileItem]:
                            if parser_ui.fileDict[fileItem].strDict:
                                if self.topostream is None:
                                    self.topostream=parser_ui.fileDict[fileItem].strDict
                                else:
                                    self.topostream.update(parser_ui.fileDict[fileItem].strDict)
                        self.load_topology()
                        if self.topohandler is not None:
                            break
        if self.topohandler is not None and self.topocode is not None:
            parser_ui.topology = self.set_TopologyData()

    def initializeTrajectoryFileHandlers(self, parser_ui):
        # Files will be loaded using their extensions initially.
        # If this fails, the fileFormat lists will be used in loading process.
        self.init_traj()
        self.trajectory = None
        trajformat = None
        self.atompositions = None
        for fileItem in parser_ui.fileDict:
            fileformatlist = None
            interfacelist = None
            if (parser_ui.fileDict[fileItem].fileSupplied and
                parser_ui.fileDict[fileItem].activeInfo):
                filename = parser_ui.fileDict[fileItem].fileName
                fileformatlist = parser_ui.fileDict[fileItem].fileFormat
                interfacelist = parser_ui.fileDict[fileItem].fileInterface
                if 'trajectory' in parser_ui.fileDict[fileItem].infoPurpose:
                    if interfacelist is not None:
                        self.interfaceorder = interfacelist
                    for fileformat in fileformatlist:
                        if self.atompositions is not None:
                            continue
                        self.set_traj(filename, fileformat)
                        if 'strDict' in parser_ui.fileDict[fileItem]:
                            if parser_ui.fileDict[fileItem].strDict:
                                self.trajstream.update(parser_ui.fileDict[fileItem].strDict)
                        traj_loaded = self.load()
                        try: 
                            self.atompositions = self.iread()
                        except TypeError:
                            pass
        if self.atompositions is not None:
            parser_ui.trajectory = self.set_TrajectoryData()

    def initializeOutputCoordinateFileHandlers(self, parser_ui):
        # Files will be loaded using their extensions initially.
        # If this fails, the fileFormat lists will be used in loading process.
        self.init_outcoord()
        trajformat = None
        self.outputpositions = None
        for fileItem in parser_ui.fileDict:
            fileformatlist = None
            interfacelist = None
            if (parser_ui.fileDict[fileItem].fileSupplied and
                parser_ui.fileDict[fileItem].activeInfo):
                filename = parser_ui.fileDict[fileItem].fileName
                fileformatlist = parser_ui.fileDict[fileItem].fileFormat
                interfacelist = parser_ui.fileDict[fileItem].fileInterface
                if 'outputcoordinates' in parser_ui.fileDict[fileItem].infoPurpose:
                    if interfacelist is not None:
                        self.interfaceorder = interfacelist
                    for fileformat in fileformatlist:
                        if self.outputpositions is not None:
                            continue
                        self.set_outcoord(filename, fileformat)
                        if 'strDict' in parser_ui.fileDict[fileItem]:
                            if parser_ui.fileDict[fileItem].strDict:
                                self.incoordstream.update(parser_ui.fileDict[fileItem].strDict)
                        outcoord_loaded = self.load_outcoord()
                        self.outputpositions = self.outcoord_iread()
        if self.outputpositions is not None:
            parser_ui.outputcoords = self.set_OutputData()

    def initializeInputCoordinateFileHandlers(self, parser_ui):
        # Files will be loaded using their extensions initially.
        # If this fails, the fileFormat lists will be used in loading process.
        self.init_incoord()
        trajformat = None
        self.inputpositions = None
        for fileItem in parser_ui.fileDict:
            fileformatlist = None
            interfacelist = None
            if (parser_ui.fileDict[fileItem].fileSupplied and
                parser_ui.fileDict[fileItem].activeInfo):
                filename = parser_ui.fileDict[fileItem].fileName
                fileformatlist = parser_ui.fileDict[fileItem].fileFormat
                interfacelist = parser_ui.fileDict[fileItem].fileInterface
                if 'inputcoordinates' in parser_ui.fileDict[fileItem].infoPurpose:
                    if interfacelist is not None:
                        self.interfaceorder = interfacelist
                    for fileformat in fileformatlist:
                        if self.inputpositions is not None:
                            continue
                        self.set_incoord(filename, fileformat)
                        if 'strDict' in parser_ui.fileDict[fileItem]:
                            if parser_ui.fileDict[fileItem].strDict:
                                if self.incoordstream is None:
                                    self.incoordstream=parser_ui.fileDict[fileItem].strDict
                                else:
                                    self.incoordstream.update(parser_ui.fileDict[fileItem].strDict)
                        incoord_loaded = self.load_incoord()
                        self.inputpositions = self.incoord_iread()
        if self.inputpositions is not None:
            parser_ui.inputcoords = self.set_InputData()

    def set_TopologyData(self):
        if self.topology is None:
            self.topology = MDDataTopology()
        if self.topocharmmpsf:
            self.topology.charmmpsf = True
        self.topology.filename = self.topofile
        self.topology.topoDict = {}
        if("gromosread" in self.topocode):
            self.topology = gtoConvertTopoDict(self.topology, self.topohandler)
        if("parmed" in self.topocode or 
           "charmmcoor" in self.topocode):
            self.topology = pmdConvertTopoDict(self.topology, self.topohandler)
        if((self.trajcode is not None or 
            self.incoordcode is not None or 
            self.outcoordcode is not None) and 
            self.topology is None):
            if "charmmcoor" in self.trajcode:
                self.topology = pmdConvertTopoDict(self.topology, self.trajhandler)
            if "charmmcoor" in self.incoordcode:
                self.topology = pmdConvertTopoDict(self.topology, self.incoordhandler)
            if "charmmcoor" in self.outcoordcode:
                self.topology = pmdConvertTopoDict(self.topology, self.outcoordhandler)
        if "pymolfile" in self.topocode:
            self.topology = pymConvertTopoDict(self.topology, self.topohandler)
        if(self.trajcode and not self.topology):
            if ("pymolfile" in self.topocode and 
                "pymolfile" in self.trajcode):
                self.topology = pymConvertTopoDict(self.topology, self.trajhandler.Topology)
        if "mdtraj" in self.topocode:
            self.topology = mdtConvertTopoDict(self.topology, self.topohandler)
        if(self.trajcode and not self.topology):
            if ("mdtraj" in self.topocode and 
                "mdtraj" in self.trajcode):
                self.topology = mdtConvertTopoDict(self.topology, self.trajhandler.get_topology())
        if "mdanalysis" in self.topocode:
            self.topology = mdaConvertTopoDict(self.topology, self.topohandler)
        if(self.topology.n_atoms is None and
           self.topology.topo_n_atoms is not None):
            self.topology.n_atoms = self.topo_n_atoms
        return self.topology

    def trajgen(self):
        if self.trajectory.frame_no == -1:
            self.trajectory.natoms = self.natoms
            self.trajectory.frame_no += 1
            return self.atompositions
        else:
            self.atompositions = self.iread()
            self.trajectory.natoms = self.natoms
            if self.atompositions is not None:
                self.trajectory.frame_no += 1
            return self.atompositions

    def incoordgen(self):
        if self.inputcoords.frame_no == -1:
            self.inputcoords.natoms = self.incoord_natoms
            return self.inputpositions
        else:
            self.inputpositions = self.incoord_iread()
            self.inputcoords.natoms = self.incoord_natoms
            if self.inputpositions is not None:
                self.inputcoords.frame_no += 1
            return self.inputpositions

    def outcoordgen(self):
        if self.outputcoords.frame_no == -1:
            self.outputcoords.natoms = self.outcoord_natoms
            return self.outputpositions
        else:
            self.outputpositions = self.outcoord_iread()
            self.outputcoords.natoms = self.outcoord_natoms
            if self.outputpositions is not None:
                self.outputcoords.frame_no += 1
            return self.outputpositions

    def set_TrajectoryData(self):
        if self.trajectory is None:
            self.trajectory = MDDataTrajectory()
        self.trajectory.filename = self.trajfile
        self.trajectory.nsteps = self.set_nsteps()
        self.trajectory.frame_no = -1
        self.trajectory.trajDict = None
        self.trajectory.positions = self.trajgen
        return self.trajectory

    def set_InputData(self):
        if self.inputcoords is None:
            self.inputcoords = MDDataTrajectory()
        self.inputcoords.frame_no = -1
        self.inputcoords.nsteps = 1
        self.inputcoords.filename = self.incoordfile
        self.inputcoords.positions = self.incoordgen
        return self.inputcoords

    def set_OutputData(self):
        if self.outputcoords is None:
            self.outputcoords = MDDataTrajectory()
        self.outputcoords.frame_no = -1
        self.outputcoords.nsteps = 1
        self.outputcoords.filename = self.outcoordfile
        self.outputcoords.positions = self.outcoordgen
        return self.outputcoords

    def set_ThermoData(self):
        try:
            self.thermohandler.keys()
            self.thermostats = MDDataThermo()
            self.thermostats.thermoDict = {} 
        except (AttributeError, IOError, OSError, ValueError, TypeError):
            return

        setflag = 0
        nframes = 1
        self.thermostats.nsteps = nframes
        for dfkey in self.thermohandler.keys():
            if setflag == 0:
                setflag = 1
                nframes = int(len(self.thermohandler[dfkey]))
                self.thermostats.nsteps = nframes
            if 'Time' in dfkey:
                self.thermostats.time = np.asarray(self.thermohandler[dfkey])
            if 'LJ (SR)' in dfkey:
                self.thermostats.ljsr = np.asarray(self.thermohandler[dfkey])
            if 'Coulomb (SR)' in dfkey:
                self.thermostats.coulombsr = np.asarray(self.thermohandler[dfkey])
            if 'Potential' in dfkey:
                self.thermostats.pe = np.asarray(self.thermohandler[dfkey])
            if 'Kinetic' in dfkey:
                self.thermostats.ke = np.asarray(self.thermohandler[dfkey])
            if 'Total Energy' in dfkey:
                self.thermostats.te = np.asarray(self.thermohandler[dfkey])
            if 'Pressure' in dfkey:
                self.thermostats.press = np.asarray(self.thermohandler[dfkey])
            if 'Temperature' in dfkey:
                self.thermostats.temp = np.asarray(self.thermohandler[dfkey])
            if 'Bond' in dfkey:
                self.thermostats.bond = np.asarray(self.thermohandler[dfkey])
            if 'U-B' in dfkey:
                self.thermostats.ub = np.asarray(self.thermohandler[dfkey])
            if 'Proper Dih' in dfkey:
                self.thermostats.proper = np.asarray(self.thermohandler[dfkey])
            if 'Improper Dih' in dfkey:
                self.thermostats.improper = np.asarray(self.thermohandler[dfkey])
            if 'CMAP Dih' in dfkey:
                self.thermostats.cmapdih = np.asarray(self.thermohandler[dfkey])
            if 'LJ-14' in dfkey:
                self.thermostats.lj = np.asarray(self.thermohandler[dfkey])
            if 'Coulomb-14' in dfkey:
                self.thermostats.coulomb = np.asarray(self.thermohandler[dfkey])
            if 'Conserved En' in dfkey:
                self.thermostats.conserveden = np.asarray(self.thermohandler[dfkey])
            if 'Constr. rmsd' in dfkey:
                self.thermostats.constrrmsd = np.asarray(self.thermohandler[dfkey])

        if 'Vir-' in self.thermohandler.keys():
            self.thermostats.virial_tensor = np.zeros(nframes,3,3)
            for step in range(nframes):
                self.thermostats.virial_tensor[step] = np.asarray([
                    df[u'Vir-XX'][time[step]], df[u'Vir-XY'][time[step]], df[u'Vir-XZ'][time[step]],
                    df[u'Vir-YX'][time[step]], df[u'Vir-YY'][time[step]], df[u'Vir-YZ'][time[step]],
                    df[u'Vir-ZX'][time[step]], df[u'Vir-ZY'][time[step]], df[u'Vir-ZZ'][time[step]]
                    ]).reshape((3, 3))
        if 'Pres-' in self.thermohandler.keys():
            self.thermostats.pressure_tensor = np.zeros(nframes,3,3)
            for step in range(nframes):
                self.thermostats.pressure_tensor[step] = np.asarray([
                    df[u'Pres-XX'][time[step]], df[u'Pres-XY'][time[step]], df[u'Pres-XZ'][time[step]],
                    df[u'Pres-YX'][time[step]], df[u'Pres-YY'][time[step]], df[u'Pres-YZ'][time[step]],
                    df[u'Pres-ZX'][time[step]], df[u'Pres-ZY'][time[step]], df[u'Pres-ZZ'][time[step]]
                    ]).reshape((3, 3))

        # pandas DataFrame columns to keys
        dfkeys = self.thermohandler.columns.get_values()
        dfid = 0
        for dfkey, val in self.thermohandler.items():
            self.thermostats.thermoDict.update({dfkeys[dfid]:np.asarray(val)})
            dfid += 1

        def thermogen(self):
            for step in range(self.thermostats.nsteps):
                self.thermostats.step_no = step
                yield step

        self.thermostats.iread = thermogen(self)
        return self.thermostats

    def get_natoms_from_topo(self, topocode=None):
        """Read the first configuration of the coordinate file to extract the
           number of atoms in it.
        """
        if(topocode is None and 
            self.topocode is not None):
            topocode = self.topocode
        n_atoms=None
        if self.topohandler is not None and topocode is not None:
            if "gromosread" in topocode:
                n_atoms = len(self.topohandler.atoms)
            elif "pymolfile" in topocode:
                n_atoms = self.topohandler.topology.natoms
            elif "mdtraj" in topocode:
                n_atoms = self.topohandler.topology.natoms
            elif "parmed" in topocode:
                if isinstance(self.topohandler, pmd.charmm.CharmmParameterSet):
                    n_atoms = len([i for i in self.topohandler.atom_types])
                elif isinstance(self.topohandler, pmd.charmm.CharmmPsfFile):
                    self.topocharmmpsf = True
                    n_atoms = len(self.topohandler.atoms)
                elif self.topology.n_atoms is not None:
                    n_atoms = self.topology.n_atoms
            elif "charmmcoor" in topocode:
                n_atoms = int(self.topology.charmmcoor['numatoms'])
            elif "mdanalysis" in topocode:
                n_atoms = len(self.topohandler.atoms)
            elif "ase" in topocode:
                n_atoms = len(self.topohandler.get_positions())
        return n_atoms

    def get_natoms_from_traj(self, trajcode=None):
        """Read the first configuration of the coordinate file to extract the
           number of atoms in it.
        """
        if(trajcode is None and 
            self.trajcode is not None):
            trajcode = self.trajcode
        n_atoms=None
        if self.trajhandler is not None:
            if "gromosread" in trajcode:
                n_atoms = self.trajhandler.natoms
            elif "pymolfile" in trajcode:
                n_atoms = self.trajhandler.natoms
            elif "mdtraj" in trajcode:
                n_atoms = self.trajiter.xyz.shape[1]
            elif "mdanalysis" in trajcode:
                n_atoms = len(self.trajiter)
            elif "ase" in trajcode:
                n_atoms = len(self.trajiter)
        return n_atoms

    def set_nsteps(self):
        #if "pymolfile" in self.trajcode:
        #    return len(self.trajhandler.nsteps)
        #if "mdtraj" in self.trajcode:
        #    return len(self.trajhandler.trajectory)
        if "mdanalysis" in self.trajcode:
            return len(self.trajhandler.trajectory)
        if "ase" in self.trajcode:
            return len(self.trajhandler.trajectory)
        return None

if __name__ == "__main__":
    topo_file = sys.argv[1]
    topo_form = sys.argv[2]
    traj_file = sys.argv[3]
    traj_form = sys.argv[4]
    MDdata = MDDataAccess()
    MDdata.trajfile = traj_file
    MDdata.trajformat = traj_form
    MDdata.topofile = topo_file
    MDdata.topoformat = topo_form
    if 'None' in topo_format:
        pass
    else:
        MDdata.topoformat = topo_format
    if 'None' in traj_format:
        pass
    else:
        MDdata.trajformat = traj_format
    MDdata.interfaceorder = ["mdanalysis", "pymolfile", "mdtraj", "ase"]
    #MDdata.interfaceorder = ["pymolfile", "mdanalysis", "mdtraj", "ase"]
    #MDdata.interfaceorder = ["charmmcoor","parmed","pymolfile", "mdanalysis", "mdtraj"]
    MDdata.interfacematch = False
    traj_iterator = MDdata.load()
    MDdata.trajchunk=300
    MDaccess = False
    print("TopoFile:" + MDdata.topofile)
    if MDdata.topocode:
        print("TopoCode:" + MDdata.topocode)
        print("Topo Handler:" + str(MDdata.topohandler))
    atoms = MDdata.topohandler
    print("TrajFile:" + MDdata.trajfile)
    if MDdata.trajcode:
        print("TrajCode:" + MDdata.trajcode)
        print("Traj Handler:" + str(MDdata.trajhandler))
        print(type(MDdata.trajhandler))
        #print(MDdata.topohandler.trajectory)

    atom = MDdata.iread()
    if atom is not None:
        readaccess = True
    print("ReadAccess:" + str(readaccess))

    i=0
    if readaccess:
        while True:
            if atom is not None:
                print(atom)
                i += 1
            else:
                break
            atom = MDdata.iread()
    print(i)

