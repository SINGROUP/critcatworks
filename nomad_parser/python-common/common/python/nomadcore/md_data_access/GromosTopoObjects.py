import numpy as np
import sys
import os
import re
import zipfile
import gzip

def get_fileExtensions(fname):
    if fname is None:
        fname = ''
    if fname != '':
        outList=[]
        if '.' in fname:
            base, ext = os.path.splitext(fname)
            if base.startswith('.'):
                base = base[1:]
            if ext.startswith('.'):
                ext = ext[1:]
            if '.' in base:
                baseList = get_fileExtensions(base)
            else:
                baseList = [base]
            if '.' in ext:
                extList = get_fileExtensions(ext)
            else:
                extList = [ext]
            for bitem in baseList:
                outList.append(bitem)
            for eitem in extList:
                outList.append(eitem)
        return outList
    else:
        return None

def get_zipType(tfile):
    basename = os.path.basename(tfile)
    exts=get_fileExtensions(basename)
    if 'gz' in exts[-1]:
        ttype = exts[-2]
        ztype = 'gz'
        zfile = basename.replace('.gz', '')
    elif 'zip' in exts[-1]:
        ttype = exts[-2]
        ztype = 'zip'
        zfile = basename.replace('.zip', '')
    else:
        ttype = exts[-1]
        ztype = None
        zfile = tfile
    return ttype, zfile, ztype

def _changes(Func):
    def newFunc(self, *args, **kwargs):
        self.need_update = True
        return Func(self, *args, **kwargs)
    return newFunc

class IndexedList(list):
    """ This list tracks the index of items 
        and assign them to its siblings _id value.
    """
    def __init__(self, *args):
        self.size = 0
        self.need_update = False
        return list.__init__(self, *args)

    def update_index(self):
        if self.need_update is True:
            for i, item in enumerate(self):
                if hasattr(item, '_id'):
                    self[i]._id=i

    @_changes
    def append(self, item):
        if hasattr(item, '_id'):
            item._id = self.size
        self.size += 1
        if hasattr(item, 'list'):
            item.list = self
        self.update_index()
        return list.append(self, item)

    @_changes
    def __delitem__(self, index):
        try:
            ind_list = range(*index.indices(len(self)))
        except AttributeError:
            ind_list = [index]

        for i in ind_list:
            if hasattr(self[i], '_id'):
                self[i]._id=-1
            if hasattr(self[i], 'list'):
                self[i].list = None

        theList = list.__delitem__(self, index)
        self.update_index()
        self.size = len(self)
        return theList

    @_changes
    def __delslice__(self, first, last):
        self.__delitem__(slice(first, last))
        self.update_index()
        self.size = len(self)

    @_changes
    def pop(self, lid=-1):
        item = list.pop(self, lid)
        if hasattr(item, '_id'):
            item._id = -1
            item.list = None
        self.update_index()
        self.size = len(self)
        return item

    @_changes
    def remove(self, item):
        try:
            list.remove(self, item)
            if hasattr(item, '_id'):
                item._id = -1
                item.list = None
            self.update_index()
            self.size = len(self)
        except AttributeError:
            pass

    @_changes
    def __setitem__(self, index, value):
        if hasattr(self[index], '_id'):
            self[index]._id=-1
        if hasattr(self[index], 'list'):
            self[index].list = None

        list.__setitem__(self, index, value)
        if hasattr(self[index], '_id'):
            self[index]._id=index
        if hasattr(self[index], 'list'):
            self[index].list = self

    #append = _changes(list.append)
    extend = _changes(list.extend)
    insert = _changes(list.insert)
    #__setitem__ = _changes(list.__setitem__)
    __iadd__ = _changes(list.__iadd__)
    __imul__ = _changes(list.__imul__)

    def __add__(self, item):
        return IndexedList(list.__add__(self, item))

    def __mul__(self, num):
        return IndexedList(list.__mul__(self, num))

class GromosTopoObject(dict):
    """ Stores data in TopoObject classes.
    """
    def __init__(self, *args, **kwargs):
        super(GromosTopoObject, self).__init__(*args, **kwargs)
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
        super(GromosTopoObject, self).__setitem__(key, value)
        self.__dict__.update({key: value})

class GromosAtomType(GromosTopoObject):

    def __init__(self, list=None, _id=-1, name='', _type='', other_type='', 
                 atom_number=0, atomic_number=-1, type_id=None, mass=0.0, charge=None, 
                 bfactor=0.0, altloc='', occupancy=0.0, solute_atom=False, 
                 solvent_atom=False, residue='', res_id=None, charge_group=None, 
                 exclude_atoms=None, lj14_atoms=None):
        self.list = list
        self._id = _id
        self.name = name.strip()
        self._type = _type.strip()
        self.type_id = type_id
        self.other_type = other_type.strip()
        self.atom_number = atom_number
        self.atomic_number = atomic_number
        self.mass = mass 
        self.charge = charge 
        self.bfactor = bfactor 
        self.altloc = altloc 
        self.occupancy = occupancy
        self.solute_atom = solute_atom 
        self.solvent_atom = solvent_atom 
        self.residue = residue
        self.res_id = res_id
        self.charge_group = charge_group
        self.exclude_atoms = exclude_atoms
        self.lj14_atoms = lj14_atoms
        self.nonbonded_list = None
        self.unique_name = ''
        if self.name != '' and self.name is not None:
            self.unique_name = self.unique_name + self.name
        if self._type != '' and self._type is not None:
            self.unique_name = self.unique_name + '_' + self._type
        if self.res_id is not None and self.residue != '':
            self.unique_name = self.unique_name + '_' + str(self.residue) + str(self.res_id)
        if self.other_type != '' and self.other_type is not None:
            self.unique_name = self.unique_name + '_' + self.other_type
        if self.unique_name == '':
            self.unique_name = None

    def __repr__(self):
        repname = '<AtomType:%s id:%d' % (self.name, self._id)
        if self.residue is not None:
            repname = repname + ', Residue:%s' % (self.residue)
        return repname + '>'

class GromosAtom(GromosTopoObject):

    def __init__(self, list=None, _id=-1, name='', atom_type=None, 
                 atom_number=None, mass=None, charge=None, bfactor=None, 
                 altloc=None, occupancy=None, solute_atom=None, 
                 solvent_atom=None, residue=None, res_id=None, charge_group=None, 
                 exclude_atoms=None, lj14_atoms=None, solvent_mol=None, 
                 mol_id=None):
        self.list = list
        self._id = _id
        self.atom_type = atom_type
        if name:
            self.name = name.strip()
        self.unique_name = ''
        if self.name != '' and self.name is not None:
            self.unique_name = self.unique_name + self.name
        if self.atom_type.name != '' and self.atom_type is not None:
            self.unique_name = self.unique_name + '_' + self.atom_type.name
        if self.atom_type._type != '' and self.atom_type is not None:
            self.unique_name = self.unique_name + '_' + self.atom_type._type
        if self.unique_name == '':
            self.unique_name = None
        if self.atom_type is not None:
            if self.mass is None: 
                self.mass = self.atom_type.mass
            if self.charge is None:
                self.charge = self.atom_type.charge
            if self.bfactor is None:
                self.bfactor = self.atom_type.bfactor
            if self.altloc is None: 
                self.altloc = self.atom_type.altloc
            if self.occupancy is None: 
                self.occupancy = self.atom_type.occupancy
            if self.solute_atom is None: 
                self.solute_atom = self.atom_type.solute_atom 
            if self.solvent_atom is None: 
                self.solvent_atom = self.atom_type.solvent_atom
            if self.residue is None: 
                self.residue = self.atom_type.residue
            if self.res_id is None: 
                self.res_id = self.atom_type.res_id
            if self.charge_group is None: 
                self.charge_group = self.atom_type.charge_group
            if self.exclude_atoms is None: 
                self.exclude_atoms = self.atom_type.exclude_atoms
            if self.lj14_atoms is None: 
                self.lj14_atoms = self.atom_type.lj14_atoms
            if self.solvent_mol is None: 
                self.solvent_mol = self.atom_type.solvent_mol

    def __repr__(self):
        repname = '<Atom:%d, %s' % (self._id, self.name)
        if self.atom_type is not None:
            repname = repname + ', Type: %s' % (self.atom_type._type)
        if self.residue:
            repname = repname + ', Residue: %s' % (self.residue)
        return repname + '>'

class GromosBondParam(GromosTopoObject):

    def __init__(self, list=None, _id=-1, name='', 
                 quartic_force=None, harmonic_force=None,
                 eq_length=None, solvent_constraint=False):
        self.list = list
        self.atoms_list = {}
        self._id = _id
        self.name = name.strip()
        self.quartic_force = quartic_force
        self.harmonic_force = harmonic_force
        self.eq_length = eq_length
        self.solvent_constraint = solvent_constraint
        self.used = False

    def __repr__(self):
        repname = '<BondParam:%d' % (self._id)
        if self.solvent_constraint:
            repname = repname + ' Solvent Constraint'
        if self.quartic_force is not None:
            repname = repname + ', CB:%f' % (self.quartic_force)
        if self.harmonic_force is not None:
            repname = repname + ', CHB:%f' % (self.harmonic_force)
        if self.eq_length is not None:
            repname = repname + ', B0:%f' % (self.eq_length)
        return repname + '>'

class GromosBondType(GromosTopoObject):

    def __init__(self, list=None, _id=-1, atom1=None, atom2=None, 
                 bond_param=None, involve_h=False):
        self.list = list
        self._id = _id
        self.atom1 = atom1
        self.atom2 = atom2
        self.atoms_list = {}
        self.bond_param = bond_param
        self.involve_h = involve_h
        self.constraint = False
        if self.bond_param is not None:
            if self.bond_param.solvent_constraint is True:
                self.constraint = True
        if self.bond_param.atoms_list is not None:
            bname = self.atom1.unique_name + ',' + self.atom2.unique_name
            if bname in self.bond_param.atoms_list:
                bname_list = self.bond_param.atoms_list[bname]
                inlist = False
                for bitem in bname_list:
                    if((self.atom1._id == bitem[0] and
                        self.atom2._id == bitem[1]) or
                       (self.atom1._id == bitem[1] and
                        self.atom2._id == bitem[0])):
                           inlist = True
                if inlist is False:
                    bname_list.append([self.atom1._id, self.atom2._id])
                    self.bond_param.atoms_list.update({bname:bname_list})
            else:
                self.bond_param.atoms_list.update({bname:[
                    [self.atom1._id, self.atom2._id]
                    ]})

    def __repr__(self):
        repname = '<BondType:%d' % (self._id)
        if self.constraint:
            repname = repname + ' Solvent Constraint'
        repname = repname + ', Param:%d AtomTypes[%d-%d]' % (
                self.bond_param._id,
                self.atom1._id, self.atom2._id)
        if self.involve_h:
            repname = repname + ' H-Bond'
        return repname + '>'

class GromosBond(GromosTopoObject):

    def __init__(self, list=None, _id=-1, atom1=None, atom2=None, 
                 bond_type=None, involve_h=False):
        self.list = list
        self._id = _id
        self.atom1 = atom1
        self.atom2 = atom2
        self.bond_type = bond_type
        self.involve_h = involve_h
        self.constraint = False
        if self.bond_type is not None:
            if self.bond_type.constraint is True:
                self.constraint = True
        if self.bond_type.atoms_list is not None:
            bname = self.atom1.unique_name + ',' + self.atom2.unique_name
            if bname in self.bond_type.atoms_list:
                bname_list = self.bond_type.atoms_list[bname]
                inlist = False
                for bitem in bname_list:
                    if((self.atom1._id == bitem[0] and
                        self.atom2._id == bitem[1]) or
                       (self.atom1._id == bitem[1] and
                        self.atom2._id == bitem[0])):
                           inlist = True
                if inlist is False:
                    bname_list.append([self.atom1._id, self.atom2._id])
                    self.bond_type.atoms_list.update({bname:bname_list})
            else:
                self.bond_type.atoms_list.update({bname:[
                    [self.atom1._id, self.atom2._id]
                    ]})

    def __repr__(self):
        repname = '<Bond:%d' % (self._id)
        if self.constraint:
            repname = repname + ' Solvent Constraint'
        repname = repname + ', Type:%d Atoms[%d-%d]' % (
                self.bond_type._id,
                self.atom1._id, self.atom2._id)
        return repname + '>'

class GromosNonbondedParam(GromosTopoObject):

    def __init__(self, list=None, _id=-1, name='', 
                 type1=None, type2=None,
                 type1_name='', type2_name='',
                 lj_r12=None, lj_r6=None,
                 lj14_r12=None, lj14_r6=None):
        self.list = list
        self._id = _id
        self.name = name.strip()
        self.type1 = type1
        self.type2 = type2
        self.type1_name = type1_name
        self.type2_name = type2_name
        #if self.type1 is not None:
        #    self.type1_name = self._atomtypes[self.type1-1]
        #if self.type2 is not None:
        #    self.type2_name = self._atomtypes[self.type2-1]
        self.unique_name = ''
        if self.type1_name != '' and self.type1 is not None:
            self.unique_name = self.type1_name
        if self.type2_name != '' and self.type2 is not None:
            self.unique_name = self.unique_name + '_' + self.type2_name
        self.lj_r12 = lj_r12
        self.lj_r6 = lj_r6
        self.lj14_r12 = lj14_r12
        self.lj14_r6 = lj14_r6

    def __repr__(self):
        repname = '<NonbondedParam:%d' % (self._id)
        if self.unique_name is not None or self.unique_name != '':
            repname = repname + ', [%s]' % (self.unique_name)
        if self.lj_r12 is not None:
            repname = repname + ', C12:%f' % (self.lj_r12)
        if self.lj_r6 is not None:
            repname = repname + ', C6:%f' % (self.lj_r6)
        if self.lj14_r12 is not None:
            repname = repname + ', CS12:%f' % (self.lj14_r12)
        if self.lj14_r6 is not None:
            repname = repname + ', CS6:%f' % (self.lj14_r6)
        return repname + '>'

class GromosNonbondedType(GromosTopoObject):

    def __init__(self, list=None, _id=-1, atom1=None, atom2=None, 
                 nonbond_param=None, solvent_mol=False):
        self.list = list
        self._id = _id
        self.atom1 = atom1
        self.atom2 = atom2
        self.nonbond_param = nonbond_param
        self.involve_h = involve_h

    def __repr__(self):
        repname = '<NonbondedType:%d Param:%d AtomTypes[%d-%d]' % (
                self._id, self.nonbond_param._id,
                self.atom1._id, self.atom2._id)
        return repname + '>'


class GromosAngleParam(GromosTopoObject):

    def __init__(self, list=None, _id=-1, name='', 
                 harmonic_cosine_force=None, harmonic_force=None,
                 eq_angle=None):
        self.list = list
        self.atoms_list = {}
        self._id = _id
        self.name = name.strip()
        self.harmonic_cosine_force = harmonic_cosine_force
        self.harmonic_force = harmonic_force
        self.eq_angle = eq_angle
        self.used = False

    def __repr__(self):
        repname = '<AngleParam:%d' % (self._id)
        if self.harmonic_cosine_force is not None:
            repname = repname + ', CT:%f' % (self.harmonic_cosine_force)
        if self.harmonic_force is not None:
            repname = repname + ', CHT:%f' % (self.harmonic_force)
        if self.eq_angle is not None:
            repname = repname + ', T0:%f' % (self.eq_angle)
        return repname + '>'

class GromosAngleType(GromosTopoObject):

    def __init__(self, list=None, _id=-1, atom1=None, atom2=None, 
                 atom3=None, angle_param=None, involve_h=False):
        self.list = list
        self._id = _id
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.angle_param = angle_param
        self.involve_h = involve_h
        angname = self.atom1.unique_name 
        angname = angname + ',' + self.atom2.unique_name 
        angname = angname + ',' + self.atom3.unique_name 
        angname_list = []
        if self.angle_param.atoms_list is not None:
            if angname in self.angle_param.atoms_list:
                angname_list = self.angle_param.atoms_list[angname]
        inlist = False
        for angitem in angname_list:
            if((self.atom1._id == angitem[0] and
                self.atom2._id == angitem[1] and
                self.atom3._id == angitem[2]) or
               (self.atom1._id == angitem[2] and
                self.atom2._id == angitem[1] and
                self.atom3._id == angitem[0])):
                   inlist = True
        if inlist is False:
            angname_list.append([
                self.atom1._id, 
                self.atom2._id,
                self.atom3._id
                ])
            self.angle_param.atoms_list.update(
                    {angname:angname_list})
        else:
            self.angle_param.atoms_list.update({angname:[
                [self.atom1._id, 
                 self.atom2._id,
                 self.atom3._id]
                ]})

    def __repr__(self):
        repname = '<AngleType:%d Param:%d AtomTypes[%d-%d-%d]' % (
                self._id, self.angle_param._id, self.atom1._id, 
                self.atom2._id, self.atom3._id)
        return repname + '>'

class GromosDihedralParam(GromosTopoObject):

    def __init__(self, list=None, _id=-1, name='', improper=False,
                 force_const=None, eq_phase_angle=None,
                 multiplicity=None):
        self.list = list
        self.dihedral_atoms_list = {}
        self.improper_atoms_list = {}
        self._id = _id
        self.name = name.strip()
        self.improper = improper
        self.force_const = force_const
        self.eq_phase_angle = eq_phase_angle
        self.multiplicity = multiplicity
        self.used = False

    def __repr__(self):
        repname = '<DihedralParam:%d' % (self._id)
        if self.improper is True:
            repname = repname + ' Improper'
        if self.force_const is not None:
            repname = repname + ', C:%f' % (self.force_const)
        if self.eq_phase_angle is not None:
            repname = repname + ', Q0/PD:%f' % (self.eq_phase_angle)
        return repname + '>'

class GromosDihedralType(GromosTopoObject):

    def __init__(self, list=None, _id=-1, atom1=None, atom2=None, 
                 atom3=None, atom4=None, dihedral_param=None, 
                 involve_h=False):
        self.list = list
        self._id = _id
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        if dihedral_param.improper:
            self.improper = True
        else:
            self.improper = False
        self.dihedral_param = dihedral_param
        self.involve_h = involve_h
        dname = self.atom1.unique_name 
        dname = dname + ',' + self.atom2.unique_name 
        dname = dname + ',' + self.atom3.unique_name 
        dname = dname + ',' + self.atom4.unique_name 
        dname_list=[]
        if self.improper:
            if self.dihedral_param.improper_atoms_list is not None:
                if dname in self.dihedral_param.improper_atoms_list:
                    dname_list = self.dihedral_param.improper_atoms_list[dname]
        else:
            if self.dihedral_param.dihedral_atoms_list is not None:
                if dname in self.dihedral_param.dihedral_atoms_list:
                    dname_list = self.dihedral_param.dihedral_atoms_list[dname]
        inlist = False
        for ditem in dname_list:
            if((self.atom1._id == ditem[0] and
                self.atom2._id == ditem[1] and
                self.atom3._id == ditem[2] and
                self.atom4._id == ditem[3]) or
               (self.atom1._id == ditem[3] and
                self.atom2._id == ditem[2] and
                self.atom3._id == ditem[1] and
                self.atom4._id == ditem[0])):
                   inlist = True
        if inlist is False:
            dname_list.append([
                self.atom1._id, 
                self.atom2._id,
                self.atom3._id,
                self.atom4._id
                ])
            if self.improper:
                self.dihedral_param.improper_atoms_list.update(
                        {dname:dname_list})
            else:
                self.dihedral_param.dihedral_atoms_list.update(
                        {dname:dname_list})
        else:
            if self.improper:
                self.dihedral_param.improper_atoms_list.update({dname:[
                    [self.atom1._id, 
                     self.atom2._id,
                     self.atom3._id,
                     self.atom4._id]
                    ]})
            else:
                self.dihedral_param.dihedral_atoms_list.update({dname:[
                    [self.atom1._id, 
                     self.atom2._id,
                     self.atom3._id,
                     self.atom4._id]
                    ]})

    def __repr__(self):
        repname = '<DihedralType:%d' % (self._id) 
        if self.improper is True:
            repname = repname + ' Improper'
        if self.dihedral_param is not None:
            repname = repname + ' Param:%d' % (self.dihedral_param._id)
        if self.improper_param is not None:
            repname = repname + ' Param:%d' % (self.improper_param._id)
        repname = repname + ' AtomTypes[%d-%d-%d-%d]' % (self.atom1._id, 
                self.atom2._id, self.atom3._id, self.atom4._id)
        return repname + '>'

class GromosTopology(GromosTopoObject):

    def __init__(self, top=None, cnf=None,
                 top_format=None, cnf_format=None): 
        self.reset_top()
        self.reset_cnf()
        self.load(topofile=top, topotype=top_format,
                  conffile=cnf, conftype=cnf_format)

    def reset_top(self):
        self.atomtypes = None
        self._atomtypes = None
        self.bondtypes = None
        self.bondparams = None
        self._bondparams = None
        self.angletypes = None
        self._angleparams = None
        self.angleparams = None
        self._dihedralparams = None
        self.dihedralparams = None
        self.improperparams = None
        self._improper_param_ids = None
        self._dihedral_param_ids = None
        self.impropertypes = None
        self.dihedraltypes = None
        self.crossterms = None
        self._ljparameters = None
        self.ljparameters = None
        self.exception_nonbonds = None
        self.residues = None
        self._residues = None
        self.solutemols = None
        self.molecules = None
        self.pressure_groups = None
        self.temperature_groups = None
        self.solventatoms = None
        self.solvent_constraints = None
        self.solvent = False
        self.physical_constants = None
        self.version = None
        self._id = None
        self.title = None

    def reset_cnf(self):
        self.conf_title = None
        self.nbounds = None
        self.time = None
        self.step = None
        self.atoms = None
        self.bonds = None
        self.angles = None
        self.dihedrals = None
        self.impropers = None
        self.positions = None
        self.velocities = None
        self.a = None # box dimensions
        self.b = None
        self.c = None
        self.alpha = None # box angles
        self.beta = None
        self.gamma = None
        self.phi = None # rotation angles
        self.theta= None
        self.psi = None
        self.X = None
        self.Y = None
        self.Z = None

    def load(self, topofile=None, topotype=None,
             conffile=None, conftype=None, 
             clear_top=False, clear_cnf=True):
        if topofile is not None:
            self.load_top(topofile, topotype, clear=clear_top)
        if conffile is not None:
            self.load_cnf(conffile, conftype, clear=clear_cnf)

    def load_top(self, topofile, topotype, clear=False):
        if clear:
            self.reset_top()
        readOK=False
        ztype=None
        zfile=None
        if topotype is None:
            topotype, zfile, ztype = get_zipType(topofile)
        if('gromostop' in topotype or
           'top' in topotype):
            typtop, zfile, ztype = get_zipType(topofile)
            if topofile is not None:
                readOK=True
        if readOK:
            if ztype is not None:
                if 'zip' in ztype:
                    with zipfile.ZipFile(topofile) as zin:
                        with zin.open(zfile, 'r') as fin:
                            self.read_top(fin, decode=True)
                elif 'gz' in ztype:
                    with gzip.open(topofile, 'rt') as fin:
                        self.read_top(fin)
            else:
                with open(topofile, 'r') as fin:
                    self.read_top(fin)
        if self._bondparams:
            if self.bondparams is None:
                self.bondparams = []
            for bparam in self._bondparams:
                if bparam.atoms_list:
                    self.bondparams.append(bparam)
        if self._angleparams:
            if self.angleparams is None:
                self.angleparams = []
            for angparam in self._angleparams:
                if angparam.atoms_list:
                    self.angleparams.append(angparam)
        if self._dihedralparams:
            if self.improperparams is None:
                self.improperparams = []
            if self.dihedralparams is None:
                self.dihedralparams = []
            for dparam in self._dihedralparams:
                if dparam.improper:
                    if dparam.improper_atoms_list:
                        self.improperparams.append(dparam)
                else:
                    if dparam.dihedral_atoms_list:
                        self.dihedralparams.append(dparam)
        if self._ljparameters:
            if self.ljparameters is None:
                self.ljparameters = []
            typenames=[atype._type for atype in self.atomtypes]
            for pari, parm in enumerate(self._ljparameters):
                if(parm.type1_name in typenames and
                   parm.type2_name in typenames):
                    self.ljparameters.append(self._ljparameters[pari])

    def load_cnf(self, conffile, conftype, clear=True):
        if clear:
            self.reset_cnf()
        readOK=False
        ztype=None
        zfile=None
        if conftype is None:
            conftype, zfile, ztype = get_zipType(conffile)
        if('gromoscnf' in conftype or
           'gromos96' in conftype or
           'cnf' in conftype or
           'g96' in conftype):
            typcnf, zfile, ztype = get_zipType(conffile)
            if conffile is not None:
                readOK=True
        if readOK:
            if ztype is not None:
                if 'zip' in ztype:
                    with zipfile.ZipFile(conffile) as zin:
                        with zin.open(zfile, 'r') as fin:
                            self.read_cnf(fin, decode=True)
                elif 'gz' in ztype:
                    with gzip.open(conffile, 'rt') as fin:
                        self.read_cnf(fin)
            else:
                with open(conffile, 'r') as fin:
                    self.read_cnf(fin)
            if self._residues is not None:
                if self.residues is None:
                    self.residues = {}
                for res in self._residues:
                    reslist=[]
                    for ai, at in enumerate(self.atoms):
                        if res in at.residue:
                            reslist.append(self.atoms[ai])
                    self.residues.update({res:reslist})

    def read_top(self, fin, decode=False):
        emptyLine = re.compile(r"^\s*$")
        physc=0
        atn=0
        atnum=0
        atset=False
        resn=0
        resnum=0
        resset=False
        satn=0
        satnum=0
        satnumset=False
        satset=False
        satnn=0
        solvn=0
        solvnum=0
        solvnumset=False
        solvset=False
        solvnn=0
        atnm = None
        mres = None
        panm = None
        iac = None
        mass = None
        cg = None
        cgc = None
        ine = -1
        inelist = []
        ine14 = -1
        ine14list = []
        btn=0
        btnum=0
        btnumset=False
        btset=False
        btnn=0
        cb = None
        cbh = None
        b0 = None
        #--
        bhn=0
        bhnum=0
        bhnumset=False
        bhset=False
        bhnn=0
        ibh = None
        jbh = None
        icbh = None
        #--
        bn=0
        bnum=0
        bnumset=False
        bset=False
        bnn=0
        ib = None
        jb = None
        ic = None
        jc = None
        icb = None
        #--
        angtn=0
        angtnum=0
        angtnumset=False
        angtset=False
        angtnn=0
        ct = None
        cth = None
        t0 = None
        #--
        anghn=0
        anghnum=0
        anghnumset=False
        anghset=False
        anghnn=0
        ith = None
        jth = None
        kth = None
        icth = None
        #--
        angn=0
        angnum=0
        angnumset=False
        angset=False
        angnn=0
        it = None
        jt = None
        kt = None
        ict = None
        #--
        imptn=0
        imptnum=0
        imptnumset=False
        imptset=False
        imptnn=0
        dihtn=0
        dihtnum=0
        dihtnumset=False
        dihtset=False
        dihtnn=0
        cq = None
        q0 = None
        np = None
        #--
        imphn=0
        imphnum=0
        imphnumset=False
        imphset=False
        imphnn=0
        dihhn=0
        dihhnum=0
        dihhnumset=False
        dihhset=False
        dihhnn=0
        iqh = None
        jqh = None
        kqh = None
        lqh = None
        icqh = None
        #--
        impn=0
        impnum=0
        impnumset=False
        impset=False
        impnn=0
        dihn=0
        dihnum=0
        dihnumset=False
        dihset=False
        dihnn=0
        iq = None
        jq = None
        kq = None
        lq = None
        icq = None
        #--
        ljn=0
        ljnum=0
        ljnumset=False
        ljset=False
        ljnn=0
        ljen=0
        ljenum=0
        ljenumset=False
        ljeset=False
        ljenn=0
        ilj = None
        jlj = None
        c12 = None
        c6 = None
        cs12 = None
        cs6 = None
        #--
        smn=0
        smnum=0
        smnumset=False
        smset=False
        smnn=0
        #--
        scn=0
        scnum=0
        scnumset=False
        scset=False
        scnn=0
        #--
        tgn=0
        tgnum=0
        tgnumset=False
        tgset=False
        tgnn=0
        #--
        pgn=0
        pgnum=0
        pgnumset=False
        pgset=False
        pgnn=0
        for line in fin:
            if decode:
                line = line.decode('utf-8')
            if emptyLine.findall(line):
                continue
            cmdLine = ' '.join([x for x in line.strip().split() if x])
            cmdLineUp = cmdLine.upper()
            if cmdLineUp.startswith('#'):
                continue
            if cmdLineUp.startswith('END'):
                section = None
            elif cmdLineUp.startswith('TITLE'):
                section = 'TITLE'
                continue
            elif cmdLineUp.startswith('PHYSICALCONSTANTS'):
                section = 'PHYSICALCONSTANTS'
                continue
            elif cmdLineUp.startswith('TOPVERSION'):
                section = 'TOPVERSION'
                continue
            elif cmdLineUp.startswith('ATOMTYPENAME'):
                section = 'ATOMTYPENAME'
                continue
            elif cmdLineUp.startswith('RESNAME'):
                section = 'RESNAME'
                continue
            elif cmdLineUp.startswith('SOLUTEATOM'):
                section = 'SOLUTEATOM'
                continue
            elif cmdLineUp.startswith('BONDSTRETCHTYPE'):
                section = 'BONDSTRETCHTYPE'
                continue
            elif cmdLineUp.startswith('BONDH'):
                section = 'BONDH'
                continue
            elif cmdLineUp.split()[0] == 'BOND':
                section = 'BOND'
                continue
            elif cmdLineUp.startswith('BONDANGLEBENDTYPE'):
                section = 'BONDANGLEBENDTYPE'
                continue
            elif cmdLineUp.startswith('BONDANGLEH'):
                section = 'BONDANGLEH'
                continue
            elif cmdLineUp.split()[0] == 'BONDANGLE':
                section = 'BONDANGLE'
                continue
            elif cmdLineUp.startswith('IMPDIHEDRALTYPE'):
                section = 'IMPDIHEDRALTYPE'
                continue
            elif cmdLineUp.startswith('IMPDIHEDRALH'):
                section = 'IMPDIHEDRALH'
                continue
            elif cmdLineUp.split()[0] == 'IMPDIHEDRAL':
                section = 'IMPDIHEDRAL'
                continue
            elif cmdLineUp.startswith('TORSDIHEDRALTYPE'):
                section = 'TORSDIHEDRALTYPE'
                continue
            elif cmdLineUp.startswith('DIHEDRALH'):
                section = 'DIHEDRALH'
                continue
            elif cmdLineUp.split()[0] == 'DIHEDRAL':
                section = 'DIHEDRAL'
                continue
            elif cmdLineUp.startswith('LJPARAMETERS'):
                section = 'LJPARAMETERS'
                continue
            elif cmdLineUp.startswith('SOLUTEMOLECULES'):
                section = 'SOLUTEMOLECULES'
                continue
            elif cmdLineUp.startswith('TEMPERATUREGROUPS'):
                section = 'TEMPERATUREGROUPS'
                continue
            elif cmdLineUp.startswith('PRESSUREGROUPS'):
                section = 'PRESSUREGROUPS'
                continue
            elif cmdLineUp.startswith('LJEXCEPTIONS'):
                section = 'LJEXCEPTIONS'
                continue
            elif cmdLineUp.startswith('SOLVENTATOM'):
                section = 'SOLVENTATOM'
                continue
            elif cmdLineUp.startswith('SOLVENTCONSTR'):
                section = 'SOLVENTCONSTR'
                continue
            if section == 'TITLE':
                if self.title is None:
                    self.title = ''
                self.title = self.title + ' ' + cmdLine
            if section == 'PHYSICALCONSTANTS':
                pok=False
                if self.physical_constants is None:
                    self.physical_constants = {}
                val = cmdLine.split()[0]
                try:
                    val = float(val)
                except(TypeError, ValueError):
                    pass
                if physc == 0:
                    pok=True
                    self.physical_constants.update({'FPEPSI' : val})
                elif physc == 1:
                    pok=True
                    self.physical_constants.update({'HBAR' : val})
                elif physc == 2:
                    pok=True
                    self.physical_constants.update({'SPDL-nm_ps' : val})
                elif physc == 3:
                    pok=True
                    self.physical_constants.update({'BOLTZ' : val})
                if pok:
                    physc += 1
            if section == 'TOPVERSION':
                self.version = cmdLine.split()[0]
            if section == 'ATOMTYPENAME':
                atok=False
                val = cmdLine.split()[0]
                if atn == 0:
                    atok=True
                    try:
                        atnum = int(val)
                        atset = True
                    except(TypeError, ValueError):
                        atset = False
                else:
                    atok=True
                    if self._atomtypes is None:
                        self._atomtypes = []
                    self._atomtypes.append(val)
                    if atset is False:
                        atnum = len(self._atomtypes)
                if atok:
                    atn += 1
            if section == 'RESNAME':
                resok=False
                val = cmdLine.split()[0]
                if resn == 0:
                    resok=True
                    try:
                        resnum = int(val)
                        resset = True
                    except(TypeError, ValueError):
                        resset = False
                else:
                    resok=True
                    if self._residues is None:
                        self._residues = []
                    self._residues.append(val)
                    if resset is False:
                        resnum = len(self._residues)
                if resok:
                    resn += 1
            if section == 'SOLUTEATOM':
                satok=False
                if satn == 0:
                    satok=True
                    try:
                        val = cmdLine.split()[0]
                        satnum = int(val)
                        satnumset = True
                    except(TypeError, ValueError):
                        satnumset = False
                else:
                    if satset is False:
                        for item in cmdLine.split():
                            if satnn == 0:
                                atnm = int(item)
                            elif satnn == 1:
                                mres = int(item)
                            elif satnn == 2:
                                panm = item
                            elif satnn == 3:
                                iac = int(item)
                            elif satnn == 4:
                                mass = float(item)
                            elif satnn == 5:
                                cg = float(item)
                            elif satnn == 6:
                                cgc = int(item)
                            elif satnn == 7:
                                ine = int(item)
                            elif satnn>7 and ine>0:
                                if satnn>7 and satnn<=7+ine:
                                    inelist.append(int(item))
                                elif satnn>7+ine and ine14<0:
                                    ine14 = int(item)
                                    if ine14<1:
                                        satset=True
                                elif((satnn>8+ine and satnn<=8+ine+ine14) and 
                                      ine14>0):
                                    ine14list.append(int(item))
                                    satset = True
                            elif satnn>7 and ine<1 and ine14<0:
                                ine14 = int(item)
                                if ine14<1:
                                    satset=True
                            elif((satnn>8 and satnn<=8+ine14) and 
                                  ine<1 and ine14>0):
                                ine14list.append(int(item))
                                satset = True
                            satnn += 1
                    if satset is True:
                        if self.atomtypes is None:
                            self.atomtypes = IndexedList()
                        self.atomtypes.append(GromosAtomType(name=panm, 
                            _type=self._atomtypes[iac-1], type_id=iac, 
                            atom_number=atnm, mass=mass, charge=cg, 
                            solute_atom=True, charge_group=cgc, 
                            res_id=mres, residue=self._residues[mres-1], 
                            exclude_atoms=inelist, lj14_atoms=ine14list))
                        satnn=0
                        atnm = None
                        mres = None
                        panm = None
                        iac = None
                        mass = None
                        cg = None
                        cgc = None
                        ine = -1
                        inelist = []
                        ine14 = -1
                        ine14list = []
                        satset = False
                    satok=True
                if satok:
                    satn += 1
            if section == 'BONDSTRETCHTYPE':
                btok=False
                if btn == 0:
                    btok=True
                    try:
                        val = cmdLine.split()[0]
                        btnum = int(val)
                        btnumset = True
                    except(TypeError, ValueError):
                        btnumset = False
                else:
                    if btset is False:
                        for item in cmdLine.split():
                            if btnn < 1:
                                cb = float(item)
                            elif btnn < 2:
                                chb = float(item)
                            elif btnn < 3:
                                b0 = float(item)
                                btset = True
                            btnn += 1
                    if btset is True:
                        if self._bondparams is None:
                            self._bondparams = IndexedList()
                        self._bondparams.append(GromosBondParam(
                            quartic_force=cb, harmonic_force=chb, 
                            eq_length=b0))
                        btnn=0
                        cb = None
                        cbh = None
                        b0 = None
                        btset = False
                    btok=True
                if btok:
                    btn += 1
            if section == 'BONDH':
                bhok=False
                if bhn == 0:
                    bhok=True
                    try:
                        val = cmdLine.split()[0]
                        bhnum = int(val)
                        bhnumset = True
                    except(TypeError, ValueError):
                        bhnumset = False
                else:
                    if bhset is False:
                        for item in cmdLine.split():
                            if bhnn < 1:
                                ibh = int(item)
                            elif bhnn < 2:
                                jbh = int(item)
                            elif bhnn < 3:
                                icbh = int(item)
                                bhset = True
                            bhnn += 1
                    if bhset is True:
                        if self.bondtypes is None:
                            self.bondtypes = IndexedList()
                        self._bondparams[icbh-1].used=True
                        self.bondtypes.append(GromosBondType(
                            atom1=self.atomtypes[ibh-1], 
                            atom2=self.atomtypes[jbh-1], 
                            bond_param=self._bondparams[icbh-1], 
                            involve_h=True
                            ))
                        bhnn=0
                        ibh = None
                        jbh = None
                        icbh = None
                        bhset = False
                    bhok=True
                if bhok:
                    bhn += 1
            if section == 'BOND':
                bok=False
                if bn == 0:
                    bok=True
                    try:
                        val = cmdLine.split()[0]
                        bnum = int(val)
                        bnumset = True
                    except(TypeError, ValueError):
                        bnumset = False
                else:
                    if bset is False:
                        for item in cmdLine.split():
                            if bnn < 1:
                                ib = int(item)
                            elif bnn < 2:
                                jb = int(item)
                            elif bnn < 3:
                                icb = int(item)
                                bset = True
                            bnn += 1
                    if bset is True:
                        if self.bondtypes is None:
                            self.bondtypes = IndexedList()
                        self._bondparams[icb-1].used=True
                        self.bondtypes.append(GromosBondType(
                            atom1=self.atomtypes[ib-1], 
                            atom2=self.atomtypes[jb-1], 
                            bond_param=self._bondparams[icb-1]
                            ))
                        bnn=0
                        ib = None
                        jb = None
                        icb = None
                        bset = False
                    bok=True
                if bok:
                    bn += 1
            if section == 'BONDANGLEBENDTYPE':
                angtok=False
                if angtn == 0:
                    angtok=True
                    try:
                        val = cmdLine.split()[0]
                        angtnum = int(val)
                        angtnumset = True
                    except(TypeError, ValueError):
                        angtnumset = False
                else:
                    if angtset is False:
                        for item in cmdLine.split():
                            if angtnn < 1:
                                ct = float(item)
                            elif angtnn < 2:
                                cht = float(item)
                            elif angtnn < 3:
                                t0 = float(item)
                                angtset = True
                            angtnn += 1
                    if angtset is True:
                        if self._angleparams is None:
                            self._angleparams = IndexedList()
                        self._angleparams.append(GromosAngleParam(
                            harmonic_cosine_force=ct, 
                            harmonic_force=cht, eq_angle=t0))
                        angtnn=0
                        ct = None
                        cht = None
                        t0 = None
                        angtset = False
                    angtok=True
                if angtok:
                    angtn += 1
            if section == 'BONDANGLEH':
                anghok=False
                if anghn == 0:
                    anghok=True
                    try:
                        val = cmdLine.split()[0]
                        anghnum = int(val)
                        anghnumset = True
                    except(TypeError, ValueError):
                        anghnumset = False
                else:
                    if anghset is False:
                        for item in cmdLine.split():
                            if anghnn < 1:
                                ith = int(item)
                            elif anghnn < 2:
                                jth = int(item)
                            elif anghnn < 3:
                                kth = int(item)
                            elif anghnn < 4:
                                icth = int(item)
                                anghset = True
                            anghnn += 1
                    if anghset is True:
                        if self.angletypes is None:
                            self.angletypes = IndexedList()
                        self._angleparams[icth-1].used=True
                        self.angletypes.append(GromosAngleType(
                            atom1=self.atomtypes[ith-1], 
                            atom2=self.atomtypes[jth-1], 
                            atom3=self.atomtypes[kth-1], 
                            angle_param=self._angleparams[icth-1],
                            involve_h=True
                            ))
                        anghnn=0
                        ith = None
                        jth = None
                        kth = None
                        icth = None
                        anghset = False
                    anghok=True
                if anghok:
                    anghn += 1
            if section == 'BONDANGLE':
                angok=False
                if angn == 0:
                    angok=True
                    try:
                        val = cmdLine.split()[0]
                        angnum = int(val)
                        angnumset = True
                    except(TypeError, ValueError):
                        angnumset = False
                else:
                    if angset is False:
                        for item in cmdLine.split():
                            if angnn < 1:
                                it = int(item)
                            elif angnn < 2:
                                jt = int(item)
                            elif angnn < 3:
                                kt = int(item)
                            elif angnn < 4:
                                ict = int(item)
                                angset = True
                            angnn += 1
                    if angset is True:
                        if self.angletypes is None:
                            self.angletypes = IndexedList()
                        self._angleparams[ict-1].used=True
                        self.angletypes.append(GromosAngleType(
                            atom1=self.atomtypes[it-1], 
                            atom2=self.atomtypes[jt-1], 
                            atom3=self.atomtypes[kt-1], 
                            angle_param=self._angleparams[ict-1]
                            ))
                        angnn=0
                        it = None
                        jt = None
                        kt = None
                        ict = None
                        angset = False
                    angok=True
                if angok:
                    angn += 1
            if section == 'IMPDIHEDRALTYPE':
                imptok=False
                if imptn == 0:
                    imptok=True
                    try:
                        val = cmdLine.split()[0]
                        imptnum = int(val)
                        imptnumset = True
                    except(TypeError, ValueError):
                        imptnumset = False
                else:
                    if imptset is False:
                        for item in cmdLine.split():
                            if imptnn < 1:
                                cq = float(item)
                            elif imptnn < 2:
                                q0 = float(item)
                                imptset = True
                            imptnn += 1
                    if imptset is True:
                        if self._dihedralparams is None:
                            self._dihedralparams = IndexedList()
                        self._dihedralparams.append(GromosDihedralParam(
                            improper=True, force_const=cq, 
                            eq_phase_angle=q0))
                        if self._improper_param_ids is None:
                            self._improper_param_ids = []
                        self._improper_param_ids.append(
                                len(self._dihedralparams)-1)
                        imptnn=0
                        cq = None
                        q1 = None
                        imptset = False
                    imptok=True
                if imptok:
                    imptn += 1
            if section == 'IMPDIHEDRALH':
                imphok=False
                if imphn == 0:
                    imphok=True
                    try:
                        val = cmdLine.split()[0]
                        imphnum = int(val)
                        imphnumset = True
                    except(TypeError, ValueError):
                        imphnumset = False
                else:
                    if imphset is False:
                        for item in cmdLine.split():
                            if imphnn < 1:
                                iqh = int(item)
                            elif imphnn < 2:
                                jqh = int(item)
                            elif imphnn < 3:
                                kqh = int(item)
                            elif imphnn < 4:
                                lqh = int(item)
                            elif imphnn < 5:
                                icqh = int(item)
                                imphset = True
                            imphnn += 1
                    if imphset is True:
                        if self.impropertypes is None:
                            self.impropertypes = IndexedList()
                        self._dihedralparams[
                                self._improper_param_ids[
                                    icqh-1]].used=True
                        self.impropertypes.append(GromosDihedralType(
                            atom1=self.atomtypes[iqh-1], 
                            atom2=self.atomtypes[jqh-1], 
                            atom3=self.atomtypes[kqh-1], 
                            atom4=self.atomtypes[lqh-1], 
                            dihedral_param=self._dihedralparams[
                                self._improper_param_ids[icqh-1]],
                            involve_h=True
                            ))
                        imphnn=0
                        iqh = None
                        jqh = None
                        kqh = None
                        lqh = None
                        icqh = None
                        imphset = False
                    imphok=True
                if imphok:
                    imphn += 1
            if section == 'IMPDIHEDRAL':
                impok=False
                if impn == 0:
                    impok=True
                    try:
                        val = cmdLine.split()[0]
                        impnum = int(val)
                        impnumset = True
                    except(TypeError, ValueError):
                        impnumset = False
                else:
                    if impset is False:
                        for item in cmdLine.split():
                            if impnn < 1:
                                iq = int(item)
                            elif impnn < 2:
                                jq = int(item)
                            elif impnn < 3:
                                kq = int(item)
                            elif impnn < 4:
                                lq = int(item)
                            elif impnn < 5:
                                icq = int(item)
                                impset = True
                            impnn += 1
                    if impset is True:
                        if self.impropertypes is None:
                            self.impropertypes = IndexedList()
                        self._dihedralparams[
                                self._improper_param_ids[
                                    icq-1]].used=True
                        self.impropertypes.append(GromosDihedralType(
                            atom1=self.atomtypes[iq-1], 
                            atom2=self.atomtypes[jq-1], 
                            atom3=self.atomtypes[kq-1], 
                            atom4=self.atomtypes[lq-1], 
                            dihedral_param=self._dihedralparams[
                                self._improper_param_ids[icq-1]]
                            ))
                        impnn=0
                        iq = None
                        jq = None
                        kq = None
                        lq = None
                        icq = None
                        impset = False
                    impok=True
                if impok:
                    impn += 1
            if section == 'TORSDIHEDRALTYPE':
                dihtok=False
                if dihtn == 0:
                    dihtok=True
                    try:
                        val = cmdLine.split()[0]
                        dihtnum = int(val)
                        dihtnumset = True
                    except(TypeError, ValueError):
                        dihtnumset = False
                else:
                    if dihtset is False:
                        for item in cmdLine.split():
                            if dihtnn < 1:
                                cq = float(item)
                            elif dihtnn < 2:
                                q0 = float(item)
                            elif dihtnn < 3:
                                np = int(item)
                                dihtset = True
                            dihtnn += 1
                    if dihtset is True:
                        if self._dihedralparams is None:
                            self._dihedralparams = IndexedList()
                        self._dihedralparams.append(GromosDihedralParam(
                            force_const=cq, eq_phase_angle=q0,
                            multiplicity=np))
                        if self._dihedral_param_ids is None:
                            self._dihedral_param_ids = []
                        self._dihedral_param_ids.append(
                                len(self._dihedralparams)-1)
                        dihtnn=0
                        cq = None
                        q1 = None
                        np = None
                        dihtset = False
                    dihtok=True
                if dihtok:
                    dihtn += 1
            if section == 'DIHEDRALH':
                dihhok=False
                if dihhn == 0:
                    dihhok=True
                    try:
                        val = cmdLine.split()[0]
                        dihhnum = int(val)
                        dihhnumset = True
                    except(TypeError, ValueError):
                        dihhnumset = False
                else:
                    if dihhset is False:
                        for item in cmdLine.split():
                            if dihhnn < 1:
                                iqh = int(item)
                            elif dihhnn < 2:
                                jqh = int(item)
                            elif dihhnn < 3:
                                kqh = int(item)
                            elif dihhnn < 4:
                                lqh = int(item)
                            elif dihhnn < 5:
                                icqh = int(item)
                                dihhset = True
                            dihhnn += 1
                    if dihhset is True:
                        if self.dihedraltypes is None:
                            self.dihedraltypes = IndexedList()
                        self._dihedralparams[
                                self._dihedral_param_ids[
                                    icqh-1]].used=True
                        self.dihedraltypes.append(GromosDihedralType(
                            atom1=self.atomtypes[iqh-1], 
                            atom2=self.atomtypes[jqh-1], 
                            atom3=self.atomtypes[kqh-1], 
                            atom4=self.atomtypes[lqh-1], 
                            dihedral_param=self._dihedralparams[
                                self._dihedral_param_ids[icqh-1]],
                            involve_h=True
                            ))
                        dihhnn=0
                        iqh = None
                        jqh = None
                        kqh = None
                        lqh = None
                        icqh = None
                        dihhset = False
                    dihhok=True
                if dihhok:
                    dihhn += 1
            if section == 'DIHEDRAL':
                dihok=False
                if dihn == 0:
                    dihok=True
                    try:
                        val = cmdLine.split()[0]
                        dihnum = int(val)
                        dihnumset = True
                    except(TypeError, ValueError):
                        dihnumset = False
                else:
                    if dihset is False:
                        for item in cmdLine.split():
                            if dihnn < 1:
                                iq = int(item)
                            elif dihnn < 2:
                                jq = int(item)
                            elif dihnn < 3:
                                kq = int(item)
                            elif dihnn < 4:
                                lq = int(item)
                            elif dihnn < 5:
                                icq = int(item)
                                dihset = True
                            dihnn += 1
                    if dihset is True:
                        if self.dihedraltypes is None:
                            self.dihedraltypes = IndexedList()
                        self._dihedralparams[
                                self._dihedral_param_ids[
                                    icq-1]].used=True
                        self.dihedraltypes.append(GromosDihedralType(
                            atom1=self.atomtypes[iq-1], 
                            atom2=self.atomtypes[jq-1], 
                            atom3=self.atomtypes[kq-1], 
                            atom4=self.atomtypes[lq-1], 
                            dihedral_param=self._dihedralparams[
                                self._dihedral_param_ids[icq-1]]
                            ))
                        dihnn=0
                        iq = None
                        jq = None
                        kq = None
                        lq = None
                        icq = None
                        dihset = False
                    dihok=True
                if dihok:
                    dihn += 1
            if section == 'LJPARAMETERS':
                ljok=False
                if ljn == 0:
                    ljok=True
                    try:
                        val = cmdLine.split()[0]
                        ljnum = int(val)
                        ljnumset = True
                    except(TypeError, ValueError):
                        ljnumset = False
                else:
                    if ljset is False:
                        for item in cmdLine.split():
                            if ljnn < 1:
                                ilj = int(item)
                            elif ljnn < 2:
                                jlj = int(item)
                            elif ljnn < 3:
                                c12 = float(item)
                            elif ljnn < 4:
                                c6 = float(item)
                            elif ljnn < 5:
                                cs12 = float(item)
                            elif ljnn < 6:
                                cs6 = float(item)
                                ljset = True
                            ljnn += 1
                    if ljset is True:
                        if self._ljparameters is None:
                            self._ljparameters = IndexedList()
                        self._ljparameters.append(
                                GromosNonbondedParam(
                                    type1=ilj-1, type2=jlj-1,
                                    type1_name=self._atomtypes[ilj-1], 
                                    type2_name=self._atomtypes[jlj-1],
                                    lj_r12=c12, lj_r6=c6,
                                    lj14_r12=cs12, lj14_r6=cs6)
                                )
                        ljnn=0
                        ilj = None
                        jlj = None
                        c12 = None
                        c6 = None
                        cs12 = None
                        cs6 = None
                        ljset = False
                    ljok=True
                if ljok:
                    ljn += 1
            if section == 'SOLUTEMOLECULES':
                smok=False
                if smn == 0:
                    smok=True
                    try:
                        val = cmdLine.split()[0]
                        smnum = int(val)
                        smnumset = True
                    except(TypeError, ValueError):
                        smnumset = False
                else:
                    if self.solmolecules is None:
                        self.solmolecules = []
                    for item in cmdLine.split():
                        self.solmolecules.append(int(item))
                    smok=True
                if smok:
                    smn += 1
            if section == 'TEMPERATUREGROUPS':
                tgok=False
                if tgn == 0:
                    tgok=True
                    try:
                        val = cmdLine.split()[0]
                        tgnum = int(val)
                        tgnumset = True
                    except(TypeError, ValueError):
                        tgnumset = False
                else:
                    if self.temperature_groups is None:
                        self.temperature_groups = []
                    for item in cmdLine.split():
                        self.temperature_groups.append(int(item))
                    tgok=True
                if tgok:
                    tgn += 1
            if section == 'PRESSUREGROUPS':
                pgok=False
                if pgn == 0:
                    pgok=True
                    try:
                        val = cmdLine.split()[0]
                        pgnum = int(val)
                        pgnumset = True
                    except(TypeError, ValueError):
                        pgnumset = False
                else:
                    if self.pressure_groups is None:
                        self.pressure_groups = []
                    for item in cmdLine.split():
                        self.pressure_groups.append(int(item))
                    pgok=True
                if pgok:
                    pgn += 1
            if section == 'LJEXCEPTIONS':
                ljeok=False
                if ljen == 0:
                    ljeok=True
                    try:
                        val = cmdLine.split()[0]
                        ljenum = int(val)
                        ljenumset = True
                    except(TypeError, ValueError):
                        ljenumset = False
                else:
                    if ljeset is False:
                        for item in cmdLine.split():
                            if ljenn < 1:
                                ilj = int(item)
                            elif ljenn < 2:
                                jlj = int(item)
                            elif ljenn < 3:
                                c12 = float(item)
                            elif ljnn < 4:
                                c6 = float(item)
                                ljeset = True
                            ljenn += 1
                    if ljeset is True:
                        if self._ljparameters is None:
                            self._ljparameters = IndexedList()
                            self._ljparameters.append(
                                    GromosNonbondedParam(
                                        type1=ilj-1, type2=jlj-1,
                                        type1_name=self._atomtypes[ilj-1], 
                                        type2_name=self._atomtypes[jlj-1],
                                        lj_r12=c12, lj_r6=c6,
                                        lj14_r12=cs12, lj14_r6=cs6)
                                    )
                        else:
                            nbfound=False
                            for ni, nb in enumerate(self._ljparameters):
                                if nb.type1 == ilj and nb.type2 == jlj:
                                    nbfound=True
                                    self._ljparameters[ni].lj_r12 = c12
                                    self._ljparameters[ni].lj_r6 = c6
                            if nbfound is False:
                                self._ljparameters.append(
                                        GromosNonbondedParam(
                                            type1=ilj-1, type2=jlj-1,
                                            type1_name=self._atomtypes[ilj-1], 
                                            type2_name=self._atomtypes[jlj-1],
                                            lj_r12=c12, lj_r6=c6)
                                        )
                        ljenn=0
                        ilj = None
                        jlj = None
                        c12 = None
                        c6 = None
                        ljeset = False
                    ljeok=True
                if ljeok:
                    ljen += 1
            if section == 'SOLVENTATOM':
                solvok=False
                if solvn == 0:
                    solvok=True
                    try:
                        val = cmdLine.split()[0]
                        solvnum = int(val)
                        solvnumset = True
                    except(TypeError, ValueError):
                        solvnumset = False
                else:
                    if solvset is False:
                        for item in cmdLine.split():
                            if solvnn == 0:
                                atnm = int(item)
                            elif solvnn == 1:
                                panm = item
                            elif solvnn == 2:
                                iac = int(item)
                            elif solvnn == 3:
                                mass = float(item)
                            elif solvnn == 4:
                                cg = float(item)
                                solvset = True
                            solvnn += 1
                    if solvset is True:
                        self.solvent = True
                        # check if SOLV in RESIDUES
                        # if not add it.
                        resfound = False
                        for resi, res in enumerate(self._residues):
                            if 'SOLV' in res:
                                mres = resi
                                resfound = True
                        if resfound is False:
                            self._residues.append('SOLV')
                            mres = len(self._residues)-1
                        if self.atomtypes is None:
                            self.atomtypes = IndexedList()
                        self.atomtypes.append(GromosAtomType(name=panm, 
                            _type=self._atomtypes[iac-1], type_id=iac, 
                            atom_number=atnm, mass=mass, charge=cg, 
                            solvent_atom=True,
                            res_id=mres, residue=self._residues[mres]))
                        solvnn=0
                        atnm = None
                        panm = None
                        iac = None
                        mass = None
                        cg = None
                        solvset = False
                    solvok=True
                if solvok:
                    solvn += 1
            if section == 'SOLVENTCONSTR':
                scok=False
                if scn == 0:
                    scok=True
                    try:
                        val = cmdLine.split()[0]
                        scnum = int(val)
                        scnumset = True
                    except(TypeError, ValueError):
                        scnumset = False
                else:
                    if scset is False:
                        for item in cmdLine.split():
                            if scnn < 1:
                                ic = int(item)
                            elif scnn < 2:
                                jc = int(item)
                            elif scnn < 3:
                                icb = float(item)
                                scset = True
                            scnn += 1
                    if scset is True:
                        if self.solvent_constraints is None:
                            self.solvent_constraints = IndexedList()
                        if self._bondparams is None:
                            self._bondparams = IndexedList()
                        # Adding this as a special bond type
                        # with solvent_constraint = True so
                        # we can track it as a pair constraint.
                        self._bondparams.append(GromosBondParam(
                            eq_length=icb, solvent_constraint=True))
                        # We are sure that it is used as a type.
                        self._bondparams[-1].used=True
                        # While we add the type as a bond, we do not 
                        # add the constraint as a bond so it will not 
                        # be miss classified as a bond. 
                        # First, we need to find the correct atom-type 
                        # indices
                        ib = None
                        jb = None
                        for atmi, atmt in enumerate(self.atomtypes):
                            if atmt.solvent_atom is True:
                                if ic == atmt.atom_number:
                                    ib = atmi+1
                                if jc == atmt.atom_number:
                                    jb = atmi+1
                        if ib is not None and jb is not None:
                            self.solvent_constraints.append(GromosBondType(
                                atom1=self.atomtypes[ib-1], 
                                atom2=self.atomtypes[jb-1], 
                                bond_param=self._bondparams[-1]
                                ))
                        scnn=0
                        ib = None
                        jb = None
                        icb = None
                        scset = False
                    scok=True
                if scok:
                    scn += 1

    def read_cnf(self, fin, decode=False):
        emptyLine = re.compile(r"^\s*$")
        gbline=0
        for line in fin:
            if decode:
                line = line.decode('utf-8')
            if emptyLine.findall(line):
                continue
            cmdLine = ' '.join([x for x in line.strip().split() if x])
            cmdLineUp = cmdLine.upper()
            if cmdLineUp.startswith('#'):
                continue
            if cmdLineUp.startswith('END'):
                section = None
            elif cmdLineUp.startswith('TITLE'):
                section = 'TITLE'
                continue
            elif cmdLineUp.startswith('TIMESTEP'):
                section = 'TIMESTEP'
                continue
            elif cmdLineUp.startswith('POSITION'):
                section = 'POSITION'
                continue
            elif cmdLineUp.startswith('LATTICESHIFTS'):
                section = 'LATTICESHIFTS'
                continue
            elif cmdLineUp.startswith('VELOCITY'):
                section = 'VELOCITY'
                continue
            elif cmdLineUp.startswith('GENBOX'):
                section = 'GENBOX'
                continue
            if section == 'TITLE':
                if self.conf_title is None:
                    self.conf_title = ''
                self.conf_title = self.conf_title + ' ' + cmdLine
            if section == 'TIMESTEP':
                if self.step is None:
                    self.step = [int(cmdLine.split()[0])]
                else:
                    self.step.append(int(cmdLine.split()[0]))
                if self.time is None:
                    self.time = [float(cmdLine.split()[1])]
                else:
                    self.time.append(float(cmdLine.split()[1]))
            if section == 'POSITION':
                if self.atoms is None:
                    self.atoms = IndexedList()
                if self.positions is None:
                    self.positions = []
                resid = int(cmdLine.split()[0])
                resname = cmdLine.split()[1].strip()
                rresname = resname.upper().replace('-', '').replace('+', '')
                atyp = cmdLine.split()[2].strip()
                aatyp = atyp.upper().replace('-', '').replace('+', '')
                anum = int(cmdLine.split()[3])
                self.positions.append([float(x) for x in cmdLine.split()[4:]])
                if 'SOLV' == resname.upper():
                    if self.atomtypes:
                        at = [atmi for atmi, atmt in enumerate(
                            self.atomtypes) if(aatyp == atmt.name.upper(
                                ).replace('-', '').replace('+', '') and 
                                resname == atmt.residue)]
                        self.atoms.append(GromosAtom(name=atyp, 
                            atom_type=self.atomtypes[at[0]], 
                            atom_number=anum, 
                            solvent_mol=True,
                            residue=resname, mol_id=resid))
                    else:
                        self.atoms.append(GromosAtom(name=atyp, 
                            atom_type=None, 
                            atom_number=anum, 
                            solvent_mol=True,
                            residue=resname, mol_id=resid))
                else:
                    if self.atomtypes:
                        at = [atmi for atmi, atmt in enumerate(
                            self.atomtypes) if(aatyp == atmt.name.upper(
                                ).replace('-', '').replace('+', '') and 
                                rresname == atmt.residue.upper(
                                    ).replace('-', '').replace('+', '') and 
                                resid == atmt.res_id)]
                        self.atoms.append(GromosAtom(name=atyp, 
                            atom_type=self.atomtypes[at[0]], 
                            atom_number=anum, 
                            residue=resname, mol_id=resid))
                    else:
                        self.atoms.append(GromosAtom(name=atyp, 
                            atom_type=None, 
                            atom_number=anum, 
                            residue=resname, mol_id=resid))
            if section == 'VELOCITY':
                if self.velocities is None:
                    self.velocities = []
                self.velocities.append([float(x) for x in cmdLine.split()[4:]])
                if self.positions is None:
                    if self.atoms is None:
                        self.atoms = IndexedList()
                    resid = int(cmdLine.split()[0])
                    resname = cmdLine.split()[1]
                    rresname = resname.upper().replace('-', '').replace('+', '')
                    atyp = cmdLine.split()[2]
                    aatyp = atyp.upper().replace('-', '').replace('+', '')
                    anum = int(cmdLine.split()[3])
                    if 'SOLV' in resname:
                        if self.atomtypes:
                            at = [atmi for atmi, atmt in enumerate(
                                self.atomtypes) if(aatyp == atmt.name.upper(
                                    ).replace('-', '').replace('+', '') and 
                                    resname == atmt.residue)]
                            self.atoms.append(GromosAtom(name=atyp, 
                                atom_type=self.atomtypes[at[0]], 
                                atom_number=anum, 
                                solvent_mol=True,
                                residue=resname, mol_id=resid))
                        else:
                            self.atoms.append(GromosAtom(name=atyp, 
                                atom_type=None, 
                                atom_number=anum, 
                                solvent_mol=True,
                                residue=resname, mol_id=resid))
                    else:
                        if self.atomtypes:
                            at = [atmi for atmi, atmt in enumerate(
                                self.atomtypes) if(aatyp == atmt.name.upper(
                                    ).replace('-', '').replace('+', '') and 
                                    rresname == atmt.residue.upper(
                                        ).replace('-', '').replace('+', '') and 
                                    resid == atmt.res_id)]
                            self.atoms.append(GromosAtom(name=atyp, 
                                atom_type=self.atomtypes[at[0]], 
                                atom_number=anum, 
                                residue=resname, mol_id=resid))
                        else:
                            self.atoms.append(GromosAtom(name=atyp, 
                                atom_type=None, 
                                atom_number=anum, 
                                residue=resname, mol_id=resid))
            if section == 'GENBOX':
                if gbline == 0:
                    self.nbounds = int(cmdLine.split()[0])
                elif gbline == 1:
                    self.a = float(cmdLine.split()[0])
                    self.b = float(cmdLine.split()[1])
                    self.c = float(cmdLine.split()[2])
                elif gbline == 2:
                    self.alpha = float(cmdLine.split()[0])
                    self.beta  = float(cmdLine.split()[1])
                    self.gamma = float(cmdLine.split()[2])
                elif gbline == 3:
                    self.phi = float(cmdLine.split()[0])
                    self.theta  = float(cmdLine.split()[1])
                    self.psi = float(cmdLine.split()[2])
                elif gbline == 4:
                    self.X = float(cmdLine.split()[0])
                    self.Y  = float(cmdLine.split()[1])
                    self.Z = float(cmdLine.split()[2])
                    gbline=0
                gbline += 1

    def __repr__(self):
        repname = '<GromosTopology' 
        if self.atoms:
            repname = repname + ' %d Atoms' % (len(self.atoms))
        if self.atomtypes:
            repname = repname + ', %d Atom Types' % (len(self.atomtypes))
        if self._residues:
            repname = repname + ', %d Residues' % (len(self._residues))
        if self.solvent:
            repname = repname + ' with Solvent'
        return repname + '>'

class GromosTrajOpen(object):

    def __init__(self, filename, handle):
        self.handle = handle
        self.filename = filename
        self.ftype = None
        self.ztype = None
        self.zfile = None
        self.ftype, self.zfile, self.ztype = get_zipType(self.filename)

    def __enter__(self):
        if self.handle.finTraj is not None:
            return self
        else:
            if self.ztype is not None:
                if 'zip' in self.ztype:
                    self.handle.zinTraj = zipfile.ZipFile(self.filename)
                    self.handle.finTraj = zin.open(self.zfile, 'r')
                    self.handle.decode = True
                    return self
                elif 'gz' in self.ztype:
                    self.handle.finTraj = gzip.open(self.filename, 'rt')
                    return self
            else:
                self.handle.finTraj = open(self.filename, 'r')
                return self

    def __exit__(self, type, value, traceback):
        if self.handle.finEOF:
            if self.handle.finTraj is not None:
                self.handle.finTraj.close()
            if self.handle.zinTraj is not None:
                self.handle.zinTraj.close()

    def iread_traj(self):
        title = None
        nbounds = None
        time = None
        step = None
        trajDict = None
        positions = None
        velocities = None
        a = None # box dimensions
        b = None
        c = None
        alpha = None # box angles
        beta = None
        gamma = None
        phi = None # rotation angles
        theta= None
        psi = None
        X = None
        Y = None
        Z = None
        section = None
        emptyLine = re.compile(r"^\s*$")
        gbline=0
        rtnOK=[False, False, False]
        while True:
            line = self.handle.finTraj.readline()
            if line:
                pass
            else:
                self.handle.finEOF = True
                break
            if self.handle.decode:
                line = line.decode('utf-8')
            if emptyLine.findall(line):
                continue
            cmdLine = ' '.join([x for x in line.strip().split() if x])
            cmdLineUp = cmdLine.upper()
            if cmdLineUp.startswith('#'):
                continue
            if cmdLineUp.startswith('END'):
                if trajDict is None:
                    trajDict = {}
                if section == "TITLE":
                    trajDict.update({'TITLE' : title})
                if section == "TIMESTEP":
                    trajDict.update({'STEP' : step})
                    trajDict.update({'TIME' : time})
                    rtnOK[0] = True
                if section == "POSITION":
                    trajDict.update({'POSITIONS' : np.asarray(positions)})
                    rtnOK[1] = True
                if section == "VELOCITY":
                    trajDict.update({'VELOCITIES' : np.asarray(velocities)})
                    rtnOK[1] = True
                if section == "GENBOX":
                    trajDict.update({'PHI' : phi})
                    trajDict.update({'THETA' : theta})
                    trajDict.update({'PSI' : psi})
                    trajDict.update({'X' : X})
                    trajDict.update({'Y' : Y})
                    trajDict.update({'Z' : Z})
                    trajDict.update({
                        'UNITCELL' : [
                            a, b, c, 
                            alpha, beta, gamma
                            ]
                        })
                    gbline=0
                    rtnOK[2] = True
                section = None
            elif cmdLineUp.startswith('TITLE'):
                section = 'TITLE'
                continue
            elif cmdLineUp.startswith('TIMESTEP'):
                section = 'TIMESTEP'
                continue
            elif cmdLineUp.startswith('POSITION'):
                section = 'POSITION'
                continue
            elif cmdLineUp.startswith('VELOCITY'):
                section = 'VELOCITY'
                continue
            elif cmdLineUp.startswith('GENBOX'):
                section = 'GENBOX'
                continue
            if section == 'TITLE':
                if title is None:
                    title = ''
                title = title + ' ' + cmdLine
            elif section == 'TIMESTEP':
                if step is None:
                    step = [int(cmdLine.split()[0])]
                else:
                    step.append(int(cmdLine.split()[0]))
                if time is None:
                    time = [float(cmdLine.split()[1])]
                else:
                    time.append(float(cmdLine.split()[1]))
            elif section == 'POSITION':
                if positions is None:
                    positions = []
                positions.append([float(x) for x in cmdLine.split()])
            elif section == 'VELOCITY':
                if velocities is None:
                    velocities = []
                velocities.append([float(x) for x in cmdLine.split()])
            elif section == 'GENBOX':
                if gbline == 0:
                    nbounds = int(cmdLine.split()[0])
                elif gbline == 1:
                    a = float(cmdLine.split()[0])
                    b = float(cmdLine.split()[1])
                    c = float(cmdLine.split()[2])
                elif gbline == 2:
                    alpha = float(cmdLine.split()[0])
                    beta  = float(cmdLine.split()[1])
                    gamma = float(cmdLine.split()[2])
                elif gbline == 3:
                    phi = float(cmdLine.split()[0])
                    theta  = float(cmdLine.split()[1])
                    psi = float(cmdLine.split()[2])
                elif gbline == 4:
                    X = float(cmdLine.split()[0])
                    Y  = float(cmdLine.split()[1])
                    Z = float(cmdLine.split()[2])
                gbline += 1
            if section is None:
                if all(rtnOK) == True:
                    return trajDict

class GromosTrajectory(GromosTopoObject):

    def __init__(self, trc=None, trv=None,
                 trc_format=None, trv_format=None): 
        self.trcfile = trc
        self.trvfile = trv
        self.trctype = trc_format
        self.trvtype = trv_format
        self.step = None
        self.time = None
        self.positions = None
        self.velocities = None
        self.unitcell = None
        self.trajiter = None
        self.trajDict = None
        self.decode = False
        self.finEOF = None
        self.finTraj = None
        self.zinTraj = None
        self.load()

    def load(self, clear=True):
        self.readTRCOK=False
        self.readTRVOK=False
        self.zctype=None
        self.zcfile=None
        self.zvtype=None
        self.zvfile=None
        if self.trctype is None:
            self.trctype, self.zcfile, self.zctype = get_zipType(self.trcfile)
        if('gromostrc' in self.trctype or
           'trc' in self.trctype):
            typtrc, self.zcfile, self.zctype = get_zipType(self.trcfile)
            if self.trcfile is not None:
                self.readTRCOK=True
        #if self.trvtype is None:
        #    self.trvtype, self.zvfile, self.zvtype = get_zipType(self.trcfile)
        #if('gromostrc' in self.trctype or
        #   'trc' in self.trctype):
        #    self.typtrv, self.zvfile, self.zvtype = get_zipType(self.trcfile)
        #    if self.trvfile is not None:
        #        self.readTRVOK=True

    def iter_read(self):
        theDict = None
        if self.readTRCOK:
            with GromosTrajOpen(self.trcfile, self) as f:
                while True:
                    theDict = f.iread_traj()
                    if theDict is not None:
                        yield theDict 
                    else:
                        break
        return None

    def iread(self):
        """Returns an iterator that goes through the given trajectory file one
           configuration at a time.
        """
        iterator_object = iter(self.iter_read())
        try:
            self.trajiter = next(iterator_object)
            if self.trajiter is not None:
                if 'STEP' in self.trajiter:
                    self.step = self.trajiter['STEP'][0]
                if 'TIME' in self.trajiter:
                    self.time = self.trajiter['TIME'][0]
                if 'POSITIONS' in self.trajiter:
                    self.positions = self.trajiter['POSITIONS']
                if 'VELOCITIES' in self.trajiter:
                    self.velocities = self.trajiter['VELOCITIES']
                if 'UNITCELL' in self.trajiter:
                    self.unitcell = self.trajiter['UNITCELL']
            return self.trajiter
        except StopIteration:
            pass
        finally:
            del iterator_object

if __name__ == "__main__":
    conffile=None
    conftype=None
    topofile=sys.argv[1]
    topotype='gromostop'
    if len(sys.argv)>2:
        conffile=sys.argv[2]
        conftype='gromoscnf'
    if len(sys.argv)>3:
        trcfile=sys.argv[3]
        trctype='gromostrc'
    #gtop = GromosTopology(top=topofile, cnf=conffile, 
    #        top_format=topotype, cnf_format=conftype)
    gtop = GromosTopology()
    gtop.load_top(topofile, topotype)
    gtop.load_cnf(conffile, conftype)
    print(gtop)
    #print(gtop.atomtypes)
    #print(gtop._residues)
    #if gtop.atomtypes:
    #    for atype in gtop.atomtypes:
    #        print([atype.unique_name]) 
    #        #print([atype.atom_number, atype.res_id, 
    #        #       atype.name, atype.type_id, 
    #        #       atype.mass, atype.charge,
    #        #       atype.charge_group, atype.exclude_atoms,
    #        #       atype.lj14_atoms])
#
#    print(gtop.atomtypes[10].list)
#    #print(gtop._bondtypes)
#    print("BONDTYPES:",gtop.bondtypes)
#    #if gtop._bondtypes:
#    #    for btype in gtop._bondtypes:
#    #        if btype.atoms_list:
#    #            print("BONDTYPE:",btype._id,btype.atoms_list) 
#    print("--------------------------------------------")
#    print("BONDPARAMS:",gtop.bondparams)
#    print("--------------------------------------------")
#    print("ANGLEPARAMS:",gtop.angleparams)
#    print("--------------------------------------------")
#    print("ANGLETYPES:",gtop.angletypes)
#    print("--------------------------------------------")
#    print("DIHPARAMS:",gtop.dihedralparams)
#    print("--------------------------------------------")
#    print("IMPPARAMS:",gtop.improperparams)
#    print("--------------------------------------------")
#    print("DIHEDRALTYPES:",gtop.dihedraltypes)
#    print("--------------------------------------------")
#    print("IMPROPERTYPES:",gtop.impropertypes)
#    print("--------------------------------------------")
#    print(gtop.atomtypes)
#    print("LJPARAMS:",gtop._ljparameters)
#    print("--------------------------------------------")
#    print("SOLVCONST:",gtop.solvent_constraints)
#    print("--------------------------------------------")
#    print("LJPARAMS:",gtop.ljparameters)
#    print("--------------------------------------------")
#    for res, rlist in gtop.residues.items():
#        print(res, rlist)
#    print("--------------------------------------------")
#    print(gtop.atoms)
    print("--------------------------------------------")
    print(np.asarray(gtop.positions))
    print("--------------------------------------------")
    gtraj = GromosTrajectory(trc=trcfile, trc_format=trctype) 
    pos = None
    while True:
        pos = gtraj.iread()
        if pos is not None:
            print("-- STEP: %s --" % (gtraj.step))
            print(gtraj.positions)
            print("UNITCELL:", gtraj.unitcell)
            print("---------------")
        else:
            break
             
