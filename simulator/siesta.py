#from __future__ import print_function
from ..atoms import *
from .. import io
from .. io import cleansymb, get_unique_symbs, convert_xyz2abc, ang2bohr
from .. units import ang2bohr
from glob import glob
import re
import sys
import subprocess

#
# SIESTA Simulation Object
#
run_keys = [
    'SystemName',
    'SystemLabel',

    # grid, scf
    'XC.functional',
    'XC.authors',
    'MeshCutoff',
    'MaxSCFIterations',
    'DM.MixingWeight',
    'DM.PulayOnFile',
    'SCF.DM.Tolerance',
    'SCF.H.Tolerance',
    'DM.UseSaveDM',
    'SCFMustConverge',
    'NeglNonOverlapInt',

    #Eigenvalue Problem
    'SolutionMethod',
    'ElectronicTemperature',

    #molecular dynamics and relaxation
    'MD.TypeOfRun',
    'MD.VariableCell',
    'MD.NumCGsteps',
    'MD.MaxCGDispl',
    'MD.MaxForceTol',
    'MD.MaxStressTol',
    'GeometryConstraints', #block
    'Diag.ParallelOverK',
    'BandLinesScale',
    'SlabDipoleCorrection',

    #output option
    'WriteCoorInitial',
    'WriteKpoints',
    'WriteEigenvalues',
    'WriteKbands',
    'WriteBands',
    'WriteDM.NetCDF',
    'WriteDMHS.NetCDF',
    'AllocReportLevel',
    'WriteMullikenPop',
    'SaveHS',
    'SaveRho',
    'SaveDeltaRho',
    'SaveElectrostaticPotential',
    'SaveTotalPotential',
    'WriteMDXmol',
    'WriteCoorXmol',
    'WriteCoorStep',
    'TS.DE.Save',
    'TS.HS.Save'
]

basis_keys = [
    'PAO.BasisType',
    'PAO.BasisSize',
    'PAO.EnergyShift',
    'PAO.SplitNorm',
    'PAO.SplitNormH',
]

kpt_keys = [
    'kgrid_Monkhorst_Pack', #block
    'PDOS.kgrid_Monkhorst_Pack',    #block
    'ProjectedDensityOfStates', #block
    'LocalDensityOfStates', #block
]

TS_keys = [
    'TS.Voltage',
    'TS.Forces',
    'TS.ChemPots',
    'TS.ChemPot.Left',
    'TS.ChemPot.Right',
    'TS.Elecs.Bulk',
    'TS.Elecs.DM.Update',
    'TS.Elecs.Eta',
    'TS.Elecs.Neglect.Principal',
    'TS.Elecs',
    'TS.Elec.Left',
    'TS.Elec.Right',
    'TS.Contours.Eq.Pole',
    'TS.Contour.c-Left',
    'TS.Contour.t-Left',
    'TS.Contour.c-Right',
    'TS.Contour.t-Right',
    'TS.Contours.nEq',
    'TS.Contour.nEq.neq',
    'TBT.Elecs.Eta',
    'TBT.Contours',
    'TBT.Contour.neq',
    'TBT.DOS.A',
    'TBT.T.Out',
    'TBT.Atoms.Device',
]

block_keys =[
    'GeometryConstraints',
    'LocalDensityOfStates',
    'kgrid_Monkhorst_Pack',
    'PDOS.kgrid_Monkhorst_Pack',
    'ProjectedDensityOfStates',
    'TS.ChemPots',
    'TS.ChemPot.Left',
    'TS.ChemPot.Right',
    'TS.Elecs',
    'TS.Elec.Left',
    'TS.Elec.Right',
    'TS.Contour.c-Left',
    'TS.Contour.t-Left',
    'TS.Contour.c-Right',
    'TS.Contour.t-Right',
    'TS.Contours.nEq',
    'TS.Contour.nEq.neq',
    'TBT.Contours',
    'TBT.Contour.neq',
    'TBT.Atoms.Device',
]

class Siesta(object):

    """
    Siesta(atoms)
    
    Class for management of SIESTA simulation.

    Parameters
    ----------
    symbol : AtomsSystem
        Class instance of AtomsSystem

    Optional parameters
    -------------------

    Example
    --------
    >>> O1 = Atom('O', Vector(0,0,0))
    >>> H1 = Atom('H', Vector(-0.6, 0.6, 0))
    >>> H2 = Atom('H', Vector( 0.6, 0.6, 0))
    >>> basis = [O1, H1, H2]
    >>> atoms = AtomsSystem(basis)
    >>> sim = s2.Siesta(atoms)
    """

    #__slots__ = ['_params', '_atoms', '_inputs']

    def __init__(self, atoms):

        if isinstance(atoms, AtomsSystem):
            self._atoms = atoms
        else:
            raise ValueError("Invaild AtomsSystem")

        self._inputs = {}
        self._run_params = {}
        self._basis_params = {}
        self._kpt_params = {}
        self._TS_params = {}
        self._block_params = {}
        for key in run_keys:
            self._run_params[key] = None
        for key in basis_keys:
            self._basis_params[key] = None
        for key in kpt_keys:
            self._kpt_params[key] = None
        for key in TS_keys:
            self._TS_params[key] = None
        for key in block_keys:
            self._block_params[key] = None


    def get_options(self, key):

        """
        print the list of available options and their default values
 
        Parameters
        ----------
 
        Optional parameters
        -------------------
 
        Example
        --------
        >>> sim.get_options()
        """
        try:
            if key in self._run_params:
                if key in self._block_params:
                    return (self._run_params[key], self._block_params[key])
                else:
                    return self._run_params[key]
            if key in self._kpt_params:
                if key in self._block_params:
                    return (self._kpt_params[key], self._block_params[key])
                else:
                    return self._kpt_params[key]
            if key in self._basis_params:
                if key in self._block_params:
                    return (self._basis_params[key], self._block_params[key])
                else:
                    return self._basis_params[key]
            if key in self._TS_params:
                if key in self._TS_params:
                    return (self._TS_params[key], self._block_params[key])
                else:
                    return self._TS_params[key]
        except KeyError:
            raise IOError('Keyword "%s" in RUN.fdf is'
                            'not known.' % key)
        except IndexError:
            raise IOError('Value missing for keyword "%s".' % key)

    def set_option(self, key, value, blockValue = None):

        """
        change the options

        available key and default values
        --------------------------------

        #1. Name and basic options
        _params = {'Name'       :'siesta',  # text
                  'Label'      :'siesta',  # text
                  'Optimization' :0,       # integer
                  'MD'           :0,       # integer
                  'Run'          :'CG',    # CG or MD
                  'cell_relax'   :0,       # integer
                  'CGsteps'      :100,     # integer
                  'ForceTol'     :0.04,    # float
                  'MDsteps'      :100,     # integer
                  'MDTimeStep'   :1.0,     # float
                  'MDInitTemp'   :0.0,     # float
                  'MDTargTemp'   :300,     # float
                  'WriteCoorStep':'.false.', # bool
         
        #2. SCF/kgrid/functional parameters
                  'kgrid'      :[1,1,1],       # 3-vector
                  'kshift'     :[0,0,0],       # 3-vector
                  'BasisType'  :'split',       # split, splitgauss, nodes, nonodes
                  'BasisSize'  :'SZ',          # SZ or MINIMAL, DZ, SZP, DZP or STANDARD
                  'EnergyShift':100,           # default: 0.02 Ry
                  'Splitnorm'  :0.15,          # default: 0.15
                  'XCfunc'     :'GGA',         # GGA or LDA
                  'XCauthor'   :'PBE',         # PBE or CA
                  'MeshCutoff' :100.0,         # float
                  'Solution'   :'Diagon',      # Diagon or OrderN
                  'MaxIt'      :300,           # integer
                  'MixingWt'   :0.2,           # float
                  'Npulay'     :0,             # integer
                  'Temp'       :300.0,         # float
         
        #3. Option for post=process
                  'LDOS'    :0,
                  'LDOSE'   :(-0.10, 0.1),
                  'Denchar' :0,
                  'PDOS'    :0,
                  'PDOSE'   :(-5,5,0.1,1001),
                  'DOS'     :0,
                  'DOSE'    :(-5,5),
                  'RHO'     :0,
                 }

        Parameters
        ----------
        key : str
            option name
        value : (various)
            option value
 
        Optional parameters
        -------------------
 
        Example
        --------
        >>> sim.set_options('kgrid', [10,10,1])
        """

        if value is True and blockValue is None and key in self._block_params:
            raise IOError('trying to change block value but block value is not given')
        
        if key in self._run_params:
            self._run_params[key] = value
            if key in self._block_params:
                self._block_params[key] = blockValue
        if key in self._basis_params:
            self._basis_params[key] = value
            if key in self._block_params:
                self._block_params[key] = blockValue
        if key in self._kpt_params:
            self._kpt_params[key] = value
            if key in self._block_params:
                self._block_params[key] = blockValue
        if key in self._TS_params:
            self._TS_params[key] = value
            if key in self._block_params:
                self._block_params[key] = blockValue
            

    def read_fdf(self, filename):
        with open(filename, 'r') as fd:
            lines = fd.readlines()
        fname = filename.split('/')[-1]
        self._inputs[fname] = lines
        self.parse_fdf(filename)

    def write_fdf(self, filename):
        if filename == 'STRUCT.fdf':
            self.write_struct()
            return
        self.generate_fdf(filename)
        fd = open(filename, 'w')
        fname = filename.split('/')[-1]
        fd.writelines(self._inputs[fname])
        fd.close()


    def parse_fdf(self, fname):
        dict_input = {'KPT.fdf' : (kpt_keys,self._kpt_params) , 'BASIS.fdf' : (basis_keys,self._basis_params), 'RUN.fdf' : (run_keys,self._run_params), 'TS.fdf':(TS_keys,self._TS_params)}
        if self._inputs[fname] is None:
            raise IOError("cannot find readable fdf file")
        
        for i, line in enumerate(self._inputs[fname]):
            try:
                line = line.replace("#", "# ")
                data = line.split()
                if len(data) == 0:
                    continue
                elif data[0][0] in ['#', '!']:
                    continue

                if "#" in data:
                    data = data[:data.index("#")]
                key = data[0]

                if key in dict_input[fname][0] and key not in self._block_params:
                    value = str(' '.join(data[1:]))
                    if value in ['T', '.true.']:
                        dict_input[fname][1][key] = True
                    elif value in ['F', '.false.']:
                        dict_input[fname][1][key] = False
                    else:
                        dict_input[fname][1][key] = value
                elif key == '%block':
                    if len(data) >= 2:
                        key = data[1]
                    if key in dict_input[fname][0]:
                        dict_input[fname][1][key] = True
                        if key in ['kgrid_Monkhorst_Pack', 'PDOS.kgrid_Monkhorst_Pack']:
                            kpt = [0, 0, 0]
                            for j in range(len(kpt)):
                                if len(self._inputs[fname][i+1+j].split()) > j:
                                    num = self._inputs[fname][i+1+j].split()[j]
                                    if num.isdecimal():
                                        kpt[j] = int(num)
                            self._block_params[key] = tuple(kpt)
                        else:
                            tmp = []
                            while not re.search(rf"%endblock", self._inputs[fname][i]):
                                i += 1
                                tmpline = self._inputs[fname][i].strip()
                                tmpline += '\n'
                                tmp.append("\t"+tmpline)
                            del tmp[-1]
                            self._block_params[key] = tmp
            except KeyError:
                raise IOError('Keyword "%s" in RUN.fdf is'
                              'not known.' % key)
            except IndexError:
                raise IOError('Value missing for keyword "%s".' % key)

    def generate_fdf(self, fname):
        dict_input = {'KPT.fdf' : self._kpt_params , 'BASIS.fdf' : self._basis_params, 'RUN.fdf' : self._run_params, 'TS.fdf' : self._TS_params}

        lines = []
        lines.append("#%s generated by NanoCore\n" %fname)
        lines.append("#(1) General system descriptors\n\n")
        
        for key, val in dict_input[fname].items():
            if val is not None:
                if val == True:
                    if key in self._block_params and self._block_params[key] is not None:
                        if key in ['kgrid_Monkhorst_Pack', 'PDOS.kgrid_Monkhorst_Pack']:
                            lines.append("\n%block {}\n".format(key))
                            lines.append("   %i   0   0   0.0\n" %(self._block_params[key][0]))
                            lines.append("   0   %i   0   0.0\n" %(self._block_params[key][1]))
                            lines.append("   0   0   %i   0.0\n" %(self._block_params[key][2]))
                            lines.append("%endblock {}\n".format(key))
                        else:
                            lines.append("\n%block {}\n".format(key))
                            for line in self._block_params[key]:
                                lines.append(line)
                            lines.append("%endblock {}\n".format(key))
                    else:
                        lines.append("%s    T\n" %(key))
                elif val == False:
                    lines.append("%s    F\n" %(key))
                else:
                    lines.append("%s    %s\n" %(key, val))
            if key == 'SystemLabel':
                if self._atoms:
                    lines.append("%include STRUCT.fdf\n")
                if any(item is not None for item in self._basis_params.values()):
                    lines.append("%include BASIS.fdf\n")
                if any(item is not None for item in self._kpt_params.values()):
                    lines.append("%include KPT.fdf\n")
                if any(item is not None for item in self._TS_params.values()):
                    lines.append("%include TS.fdf\n")
        
        self._inputs[fname] = lines
        return lines

    def run(self, mode, nproc, **option):

        """
        Run a simulation based on the information saved in this simulation object
 
        Parameters
        ----------

        Optional parameters
        -------------------
        mode : 'SCF', 'MD', 'Optimization', or 'POST'
            simulation type 
        cellparameter : float
            cell expansion or compression
        log : 0 or 1
            if true, standard outputs are saved in 'stdout.txt'.
        psf : 0 or 1
            if true, pseudopotentials are copied from pre-defined path.
 
        Example
        --------
        >>> sim.get_options()
        """

        # get the location of executable
        from nanocore.env import siesta_calculator as siesta_exec
        from nanocore.env import siesta_util_tbtrans as tbtrans

        def set_transiesta_option(**opt):
            _, left_elec = self.get_options('TS.Elec.Left')
            _, right_elec = self.get_options('TS.Elec.Right')
            for i, line in enumerate(left_elec):
                if 'HS' in line:
                    left_elec[i] = f'\tHS {opt["label_L"]}\n'
                if 'used-atoms' in line:
                    left_elec[i] = f'\tused-atoms {opt["n_left"]}\n'
            
            for i, line in enumerate(right_elec):
                if 'HS' in line:
                    right_elec[i] = f'\tHS {opt["label_R"]}\n'
                if 'used-atoms' in line:
                    right_elec[i] = f'\tused-atoms {opt["n_right"]}\n'
            
            self.set_option('TS.Elec.Left', True, left_elec)
            self.set_option('TS.Elec.Right', True, right_elec)

        def set_tbtrans_option(**opt):
            _, tbt = self.get_options('TBT.Contour.neq')
            for i, line in enumerate(tbt):
                if 'from' in line:
                    tbt[i] = f'\tfrom {opt["Emin"]} to {opt["Emax"]}\n'
                if 'delta' in line:
                    tbt[i] = f'\tdelta {opt["dE"]}\n'
            
            self.set_option('TBT.Contour.neq', True, tbt)
            

        if mode == 'elec': 
            exec = siesta_exec
            self.set_option('SolutionMethod', 'Diagon')
            self.set_option('TS.DE.SAVE', True)
            self.set_option('TS.HS.SAVE', True)

        elif mode == 'scatter':
            exec = siesta_exec
            self.set_option('SolutionMethod', 'Transiesta')
            self.set_option('TS.Voltage', option["Voltage"])
            self.set_option('TS.Elecs.Eta', option["ts_eta"])
            self.set_option('TS.DE.SAVE', None)
            self.set_option('TS.HS.SAVE', None)
            set_transiesta_option(**option)

        elif mode == 'tbtrans':
            exec = tbtrans
            self.set_option('SolutionMethod', 'Transiesta')
            self.set_option('TBT.Elecs.Eta', option["tbt_eta"])
            self.set_option('TBT.Atoms.Device', True, option["device"])
            set_tbtrans_option(**option)

        # write fdf files
        if not mode == 'POST':
            files = ['RUN.fdf', 'BASIS.fdf', 'TS.fdf', 'STRUCT.fdf', 'KPT.fdf']
            for f in files:
                self.write_fdf(f)

        cmd = f'mpirun -np {nproc}  {exec} < RUN.fdf > stdout.txt'
        result = subprocess.run(cmd, shell=True, check=True)
        return result
    
    def model_opt(self, electrode_structure, channel_structure, supercell=(1,1), d_x=(0, 0, 0), d_y=(0, 0, 0), d_z=(0, 0, 0)):
        # read structures from external files
        mos2 = io.read_xsf(channel_structure)    # 2d pbc on xz plane
        slab = io.read_xsf(electrode_structure) # 3d pbc 

        # get cell matching
        n1, n2 = slab.get_cell_match(mos2, 'y', 'y', max_unit=5)

        # set repeat unit (USER INPUT)
        m1, m2 = supercell[0], supercell[1]

        # generate supercells
        slab_2 = slab * (m1, n1, 1)
        mos2_2 = mos2 * (1, n2, m2)

        # adjust cell size
        ratio = mos2_2.get_cell()[1][1] / slab_2.get_cell()[1][1]
        slab_3 = slab_2.adjust_cell_size(ratio, 7)
        slab_4 = slab_2.adjust_cell_size(ratio, 7)

        for d in np.arange(d_z):
            os.system('mkdir %.2f'%d)
            zmax = slab_3.get_zmax()       # location of the surface
            mos2_2.select_all()
            mos2_2.translate(0, 0, zmax+d) # move MoS2 on the surface
            zmax = mos2_2.get_zmax()
            slab_4.select_all()
            slab_4.translate(0, 0, zmax+d)

            # add two structures
            new = slab_3 + mos2_2 + slab_4  # merge two systems with PBC of "slab"

            # set vacuum along z axis
            zmax = slab_4.get_zmax()
            new.set_vacuum(zmax-new._cell[2][2]+15, 'z')

            # sort atoms by z coordinates
            new.select_all()
            new.sort('z')


            # save the interface model
            io.write_xsf('model.xsf', new)
            os.system('mv model.xsf %.2f'%d)    

    def write_struct(self, cellparameter=1.0):

        if self._atoms.get_cell() is not 'None':
            cell1 = self._atoms.get_cell()[0]
            cell2 = self._atoms.get_cell()[1]
            cell3 = self._atoms.get_cell()[2]

        #---------------STRUCT.fdf----------------
        fileS = open('STRUCT.fdf', 'w')
        natm = len(self._atoms)
        fileS.write("NumberOfAtoms    %d           # Number of atoms\n" % natm)
        unique_symbs = get_unique_symbs(self._atoms)
        fileS.write("NumberOfSpecies  %d           # Number of species\n\n" % len(unique_symbs))
        fileS.write("%block ChemicalSpeciesLabel\n")
    
        for symb in unique_symbs:
            sym = ''.join(re.findall('[a-zA-Z]', symb))
            fileS.write(" %d %d %s\n" % (unique_symbs.index(symb)+1,atomic_number(sym),symb) )
        fileS.write("%endblock ChemicalSpeciesLabel\n")
    
        #Lattice
        fileS.write("\n#(3) Lattice, coordinates, k-sampling\n\n")
        fileS.write("LatticeConstant   %15.9f Ang\n" % cellparameter)
        fileS.write("%block LatticeVectors\n")
        if self._atoms.get_cell() is not 'None':
            va, vb, vc = cell1, cell2, cell3
            fileS.write("%15.9f %15.9f %15.9f\n" % tuple(va))
            fileS.write("%15.9f %15.9f %15.9f\n" % tuple(vb))
            fileS.write("%15.9f %15.9f %15.9f\n" % tuple(vc))
        fileS.write("%endblock LatticeVectors\n\n")
    
        #Coordinates
        fileS.write("AtomicCoordinatesFormat Ang\n")
        fileS.write("%block AtomicCoordinatesAndAtomicSpecies\n")
    
        for atom in self._atoms:
            x,y,z = atom.get_position(); symb = atom.get_symbol()
            if atom.get_groupid() != 0:
                symb = symb + str(atom.get_groupid())
            fileS.write(" %15.9f %15.9f %15.9f %4d %4d\n" %\
                       (x,y,z,unique_symbs.index(symb)+1, atom.get_serial()))
            
        fileS.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")
        fileS.close()


    def save_simulation(self):

        import pickle
        name = 'sim_%s.dat' % self._params['Label']
        pickle.dump([self._params, 
                     self._atoms.get_symbols(),
                     np.array(self._atoms.get_positions()),
                     np.array(self._atoms.get_cell())], 
                    open(name,'w'))
        print ("simulation information is saved as %s." % name)

#
# REload saved simulation
#

def load_simulation(filename):
    import pickle
    params, symbols, positions, cell = pickle.load(open(filename))
    atoms = []
    i = 0
    for symbol in symbols:
        atoms.append(Atom(symbol, [positions[i][0], positions[i][1], positions[i][2]]))
        i += 1
    atoms = AtomsSystem(atoms, cell=cell)
    sim = Siesta(atoms)
    for key, value in params.items():
        sim.set_option(key, value)    
    return sim


#
# SIESTA UTIL interface
#

def calc_pldos(nmesh, emin, emax, npoints, orbital_index, label = 'siesta', mpi = 0, nproc = 1):

    """
    Interface to PyProjection of siesta utils

    Parameters
    ----------
    nmesh  : list
        number of grid points along  each lattice vector
    emin   : float
        minimum energy value of PLDOS plot
    emax   : float
        maximum energy value of PLDOS plot

    Optional parameters
    -------------------
    orbital_index : int
        orbital index for PDOS plot (e.g. C_2_1_1, 1_5_2, ...)
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints

    Example
    --------
    >>> calc_pldos(nmesh = [10, 10 10], emin = -10, emax =5, npoints = 1001,
                            orbital_index= 'C_2_1_0, label = 'siesta')
    """

    # write input file
    file_INP = open('PyProjection.fdf', 'w')
    file_INP.write('SystemLable  %s\n' % label)
    file_INP.write('PyProjection.TypeOfRun    PLDOS\n')
    file_INP.write('PyProjection.MinE         %4.2f\n'%emin)
    file_INP.write('PyProjection.MaxE         %4.2f\n'%emax)
    file_INP.write('PyProjection.NumE         %4.2f\n'%npoints)
    file_INP.write('PyProjection.NumA         %d\n'%nmesh[0])
    file_INP.write('PyProjection.NumB         %d\n'%nmesh[1])
    file_INP.write('PyProjection.NumC         %d\n'%nmesh[2])
    file_INP.write('PyProjection.TargetOrbital    %s\n'%orbital_index)

    from NanoCore.env import siesta_pyprojection as pdos
    cmd = 'python %s' % pdos
    if mpi:
        cmd = 'mpirun -np %i ' % nproc + cmd
    os.system(cmd)


def calc_pdos(nmesh, emin, emax, npoints, orbital_index, label = 'siesta', mpi = 0, nporc = 1):

    """
    Interface to PyProjection of siesta utils

    Parameters
    ----------
    nmesh  : list
        number of grid points along  each lattice vector
    emin   : float
        minimum energy value of PLDOS plot
    emax   : float
        maximum energy value of PLDOS plot

    Optional parameters
    -------------------
    orbital_index : int
        orbital index for PDOS plot (e.g. C_2_1_1, 1_5_2, ...)
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints

    Example
    --------
    >>> calc_pdos(nmesh = [10, 10 10], emin = -10, emax =5, npoints = 1001,
                            orbital_index= 'C_2_1_0, label = 'siesta')
    """

    # write input file
    file_INP = open('PyProjection.fdf', 'w')
    file_INP.write('SystemLable  %s\n' % label)
    file_INP.write('PyProjection.TypeOfRun    PDOS\n')
    file_INP.write('PyProjection.MinE         %4.2f\n'%emin)
    file_INP.write('PyProjection.MaxE         %4.2f\n'%emax)
    file_INP.write('PyProjection.NumE         %4.2f\n'%npoints)
    file_INP.write('PyProjection.NumA         %d\n'%nmesh[0])
    file_INP.write('PyProjection.NumB         %d\n'%nmesh[1])
    file_INP.write('PyProjection.NumC         %d\n'%nmesh[2])
    file_INP.write('PyProjection.TargetOrbital    %s\n'%orbital_index)

    from NanoCore.env import siesta_pyprojection as pdos
    cmd = 'python %s' % pdos
    if mpi:
        cmd = 'mpirun -np %i ' % nproc + cmd
    os.system(cmd)


def calc_fatband(nmesh, emin, emax, npoints, orbital_index, label = 'siesta', mpi = 0, nporc = 1):

    """
    Interface to PyProjection of siesta utils

    Parameters
    ----------
    nmesh  : list
        number of grid points along  each lattice vector
    emin   : float
        minimum energy value of PLDOS plot
    emax   : float
        maximum energy value of PLDOS plot

    Optional parameters
    -------------------
    orbital_index : int
        orbital index for PDOS plot (e.g. C_2_1_1, 1_5_2, ...)
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints

    Example
    --------
    >>> calc_fatband(nmesh = [10, 10 10], emin = -10, emax =5, npoints = 1001,
                            orbital_index= 'C_2_1_0, label = 'siesta')
    """

    # write input file
    file_INP = open('PyProjection.fdf', 'w')
    file_INP.write('SystemLable  %s\n' % label)
    file_INP.write('PyProjection.TypeOfRun    FAT\n')
    file_INP.write('PyProjection.MinE         %4.2f\n'%emin)
    file_INP.write('PyProjection.MaxE         %4.2f\n'%emax)
    file_INP.write('PyProjection.NumE         %4.2f\n'%npoints)
    file_INP.write('PyProjection.NumA         %d\n'%nmesh[0])
    file_INP.write('PyProjection.NumB         %d\n'%nmesh[1])
    file_INP.write('PyProjection.NumC         %d\n'%nmesh[2])
    file_INP.write('PyProjection.TargetOrbital    %s\n'%orbital_index)

    from NanoCore.env import siesta_pyprojection as fat
    cmd = 'python %s' % fat
    if mpi:
        cmd = 'mpirun -np %i ' % nproc + cmd
    os.system(cmd)



def get_dos(emin, emax, npoints=1001, broad=0.05, label='siesta'):

    """
    Interface to Eig2dos of siesta utils

    Parameters
    ----------

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)
    emin : float
        minimum value of DOS plot
    emax : float
        maximum value of DOS plot
    npoints : int
        the number of datapoints
    broad : float
        broadening factor for DOS plot

    Example
    --------
    >>> E, dos, dos1, dos2 = s2.get_dos(-10, 10, npoints=1001, broad=0.1)
    """

    # Eig2DOS script
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_dos as sud
    os.system('%s/%s -f -s %f -n %i -m %f -M %f %s.EIG > DOS' % (sul, sud, 
                                                                 broad, npoints, 
                                                                 emin, emax, label))

    # reload DOS
    f_dos = open('DOS').readlines()

    energy = []; dos_1 = []; dos_2 = []; dos = []
    for line in f_dos:
        if not line.startswith('#'):
            e, du, dn, dt = line.split()
            e = float(e); du = float(du); dn = float(dn); dt = float(dt)
            energy.append(e); dos_1.append(du); dos_2.append(dn); dos.append(dt)

    return energy, dos, dos_1, dos_2


def get_band(simobj, pathfile, label='siesta', rerun=0):

    """
    Interface to new.gnubands of siesta utils

    Parameters
    ----------
    simobj : simulation object
        SIESTA simulation objects, need for re-run
    pathfile : str
        the location of bandpath file (SIESTA format)

        example)

            BandLinesScale    pi/a
            %block BandLines
             1  0.000  0.000  0.000  \Gamma   # Begin at gamma
            50  0.816  0.000  0.000  M        # 50 points from gamma to M
            25  0.816  0.471  0.000  K        # 25 points from M to K
            60  0.000  0.000  0.000  \Gamma   # 60 points from K to gamma
            %endblock BandLines

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> # for spin-unpolarized cases
    >>> path, eigs = get_band(sim, './bandline')
    >>>
    >>> # for spin-polarized cases
    >>> path1, eig1, path2, eigs2 = get_band(sim, './bandline')
    """

    if rerun:
        # attach path file
        f = open('RUN.fdf', 'a')
        path = open(pathfile).readlines()
        for line in path: f.write(line)
        f.write('WriteBands            T')
        f.close()
        # re-run
        simobj.run(mode='POST')

    # gnuband script
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_band as sub
    os.system('%s/%s < %s.bands > BAND' % (sul, sub, label))

    # read band data
    from . import DataTnBand as dtb
    from glob import glob
    dtb.DataTnBand('%s.bands' % label)
    fs  = glob('band???.oneD'); fs.sort()
    fs1 = glob('band_*spin1.oneD'); fs1.sort()
    fs2 = glob('band_*spin2.oneD'); fs2.sort()

    # case 1: spin-unpolarized
    if not fs1:
        kptss = []; eigss = []
        for f in fs:
            temp = open(f).readlines()
            kpts = []; eigs = []
            for line in temp:
                if not line.startswith('#'):
                    kpt, eig = line.split()
                    kpt = float(kpt); eig = float(eig)
                    kpts.append(kpt); eigs.append(eig)
            kptss.append(kpts); eigss.append(eigs)

        os.system('rm *.oneD')
        return kptss, eigss

    # case 2: spin-polarized
    else:
        # spin 1
        kptss1 = []; eigss1 = []
        for f in fs1:
            temp = open(f).readlines()
            kpts = []; eigs = []
            for line in temp:
                if not line.startswith('#'):
                    kpt, eig = line.split()
                    kpt = float(kpt); eig = float(eig)
                    kpts.append(kpt); eigs.append(eig)
            kptss1.append(kpts); eigss1.append(eigs)

        # spin 2
        kptss2 = []; eigss2 = []
        for f in fs2:
            temp = open(f).readlines()
            kpts = []; eigs = []
            for line in temp:
                if not line.startswith('#'):
                    kpt, eig = line.split()
                    kpt = float(kpt); eig = float(eig)
                    kpts.append(kpt); eigs.append(eig)
            kptss2.append(kpts); eigss2.append(eigs)

        os.system('rm *.oneD')
        return kptss1, eigss1, kptss2, eigss2


def siesta_xsf2cube(f_in, grid_type):

    from NanoCore.io import ang2bohr

    # read file
    lines = open(f_in).readlines()

    # data
    atoms_block = []
    data_grid_blocks = []
    mesh_size = []
    orgin_point = []
    cell = []
    grid_data = []

    i_data = 0
    i = 0

    for line in lines:
        line_sp = line.split()

        # symbols, positions
        if len(line_sp) == 4:
            if i < 1000: atoms_block.append(line)

        else:
            if 'BEGIN_DATAGRID' in line:

                # initialize grid data
                grid_data = []
                i_data += 1

                # origin points: not used
                orgx, orgy, orgz = lines[i+2].split()
                orgx = float(orgx); orgy = float(orgy); orgz = float(orgz)
                origin_point = [orgx, orgy, orgz]

                # mesh: not used
                ngridx, ngridy, ngridz = lines[i+1].split()
                ngridx = int(ngridx); ngridy = int(ngridy); ngridz = int(ngridz)
                mesh_size = [ngridx, ngridy, ngridz]

                # cell vertors
                v11, v12, v13 = lines[i+3].split(); v11 = float(v11); v12 = float(v12); v13 = float(v13)
                v21, v22, v23 = lines[i+4].split(); v21 = float(v21); v22 = float(v22); v23 = float(v23)
                v31, v32, v33 = lines[i+5].split(); v31 = float(v31); v32 = float(v32); v33 = float(v33)

                cell = [[v11, v12, v13],
                        [v21, v22, v23],
                        [v31, v32, v33]]

                # write atoms
                atoms = []
                for atom_line in atoms_block:
                    symb, x, y, z = atom_line.split()
                    symb = int(symb)
                    x = float(x); y = float(y); z = float(z)
                    atoms.append( Atom(atomic_symbol[symb], [x,y,z]) )

                atoms = AtomsSystem(atoms, cell=cell)
                if grid_type =='LDOS':  filename_out = 'LDOS_%i.xsf' % i_data
                elif grid_type =='RHO': filename_out = 'RHO_%i.xsf' % i_data 
                io.write_xsf(filename_out, atoms)

                # grid data
                npoints = mesh_size[0] * mesh_size[1] * mesh_size[2]
                remain  = npoints % 6
                nlines = 0
                if remain:
                    nlines = npoints/6 + 1
                else:
                    nlines = npoints/6

                # write head
                f_out = open(filename_out, 'a')
                f_out.write('BEGIN_BLOCK_DATAGRID_3D\n')
                f_out.write('DATA_from:siesta.%s\n' % grid_type)
                f_out.write('BEGIN_DATAGRID_3D_RHO:spin_%i\n' % i_data)
                f_out.write('%8i %8i %8i\n' % tuple(mesh_size))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(origin_point))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[0]))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[1]))
                f_out.write('%12.8f%12.8f%12.8f\n' % tuple(cell[2]))

                # select data block
                for line_temp in lines[i+6:i+6+nlines]:
                    f_out.write(line_temp)
                f_out.close()

                #end for line_temp in lines[i+6:i+6+nlines+1]:
            #end if 'BEGIN_DATAGRID' in line:
        #end else:
        i += 1
    #end for line in lines:
#end def


def get_ldos(v1, v2, v3, origin, nmesh, label='siesta'):

    """
    Interface to rho2xsf of siesta utils: LDOS

    Parameters
    ----------
    v1, v2, v3 : Vector object or (3,) float array
        define the space
    origin : Vector object or (3,) float array
        define the origin of the space
    nmesh : (3,) int array
        define the number of mesh points along v1, v2, and v3

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> s2.get_ldos(v1, v2, v3, origin, nmesh)
    """

    # add block
    #f = open('RUN.fdf', 'a')
    #f.write('%block LocalDensityOfStates\n')
    #f.write('%8.4f %8.4f  eV\n' % (emin, emax))
    #f.write('%endblock LocalDensityOfStates\n')
    #f.close()

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s\n' % label)                           # 1.label
    file_INP.write('A\n')                                    # 2.unit: Ang
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(origin)) # 3.origin point
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v1))     # 4.spaning vector1
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v2))     # 5.spaning vector2
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v3))     # 6.spaning vector3
    file_INP.write('%-5i %-5i %-5i\n' % tuple(nmesh))        # 7.grid points
    file_INP.write('LDOS\n')                                 # 8-1.LDOS
    file_INP.write('BYE\n')
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_rho as sur
    os.system('%s/%s < INP' % (sul, sur))

    # convert 
    os.system('rm INP')
    os.system('mv %s.XSF LDOS.XSF' % label)
    #siesta_xsf2cube('siesta.XSF', grid_type)


def get_rho(v1, v2, v3, origin, nmesh, label='siesta'):

    """
    Interface to rho2xsf of siesta utils: RHO

    Parameters
    ----------
    v1, v2, v3 : Vector object or (3,) float array
        define the space
    origin : Vector object or (3,) float array
        define the origin of the space
    nmesh : (3,) int array
        define the number of mesh points along v1, v2, and v3

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)

    Example
    --------
    >>> s2.get_rho(v1, v2, v3, origin, nmesh)
    """

    # add keyword
    #f = open('RUN.fdf', 'a')
    #f.write('SaveRho   .true.\n')
    #f.close()

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s\n' % label)                           # 1.label
    file_INP.write('A\n')                                    # 2.unit: Ang
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(origin)) # 3.origin point
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v1))     # 4.spaning vector1
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v2))     # 5.spaning vector2
    file_INP.write('%-8.4f %-8.4f %-8.4f\n' % tuple(v3))     # 6.spaning vector3
    file_INP.write('%-5i %-5i %-5i\n' % tuple(nmesh))        # 7.grid points
    file_INP.write('RHO\n')                                  # 8-2.RHO
    file_INP.write('BYE\n')
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_rho as sur
    os.system('%s/%s < INP > OUT' % (sul, sur))

    # convert 
    os.system('rm INP OUT')
    os.system('mv %s.XSF RHO.XSF' % label)


def get_pdos(simobj, emin, emax, by_atom=1, atom_index=[], species=[], broad=0.1, npoints=1001, label='siesta'):

    """
    Interface to fmpdos of siesta utils

    Parameters
    ----------
    simobj : simulation object
        SIESTA simulation objects, need for re-run
    emin : float
        minimum value of DOS plot
    emax : float
        maximum value of DOS plot

    Optional parameters
    -------------------
    by_atom : bool
        if true,  set atom_index=[]
        if false, set species=[]
    atom_index : list
        serial numbers for PDOS plot
    species : list
        atomic symbols for PDOS plot
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints
    broad : float
        broadening factor for DOS plot

    Example
    --------
    >>> E, dos1, dos2 = get_pdos(simobj, -10, 10, by_atom=1,
                                 atom_index=[1,2,3,4], 
                                 broad=0.05, npoints=1001)

    >>> E, dos1, dos2 = get_pdos(simobj, -10, 10, by_atom=0,
                                 species=['C','H'], 
                                 broad=0.05, npoints=1001)
    """

    # re-run
    #simobj.run(mode='POST')

    # temp. input file for rho2xsf
    file_INP = open('INP', 'w')
    file_INP.write('%s.PDOS\n' % label)                      # 1.input PDOS
    file_INP.write('PDOS\n')                                 # 2.output file

    if atom_index:
        tmp_str = ''
        for ind in atom_index: tmp_str += '%i ' % ind
        file_INP.write(tmp_str + '\n')

    if species:
        tmp_str = ''
        for spe in species: tmp_str += '%s ' % spe
        file_INP.write(tmp_str + '\n')

    file_INP.write('0 \n')                                   # all quantum numbers
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_pdos as sup
    os.system('%s/%s < INP > OUT' % (sul, sup))
    os.system('rm INP OUT')
   
    # read PDOS
    lines = open('PDOS').readlines()
    energy = []; dos_1 = []; dos_2 = []

    for line in lines:
        if not line.startswith('#'):
            tmp = line.split()
            if len(tmp) == 2:
                e = float(tmp[0]); d = float(tmp[1])
                energy.append(e); dos_1.append(d)
            if len(tmp) == 3:
                e = float(tmp[0]); d1 = float(tmp[1]); d2 = float(tmp[2])
                energy.append(e); dos_1.append(d1); dos_2.append(d2)
    os.system('rm PDOS')

    return energy, dos_1, dos_2


def get_pldos(simobj, emin, emax, broad=0.1, npoints=1001, label='siesta'):

    """
    Interface to fmpdos of siesta utils

    Parameters
    ----------
    simobj : simulation object
        SIESTA simulation objects, need for re-run
    emin : float
        minimum value of DOS plot
    emax : float
        maximum value of DOS plot

    Optional parameters
    -------------------
    label : string
        label name (*.DM, *.XV, ...)
    npoints : int
        the number of datapoints
    broad : float
        broadening factor for DOS plot

    Example
    --------
    >>> z_coords, Z, E = s2.get_pldos(sim, -5, 5, 
                                      broad=0.05, 
                                      npoints=1001)
    """

    # AtomsSystem from simulation object
    atoms = simobj._atoms.copy()
    
    # slice atoms by z coordinates
    z_coords = []; indice = []
    for atom in atoms:
        if not atom[2] in z_coords: z_coords.append(atom[2])

    for z in z_coords:
        temp = []
        for atom in atoms:
            if abs(z-atom[2]) < 0.01: temp.append(atom.get_serial())
        indice.append(temp)

    # get pdos
    Z = []; E = []
    for ind in indice:
        E1, dos11, dos12 = get_pdos(simobj, emin, emax, by_atom=1, 
                                    atom_index=ind, broad=broad, npoints=npoints, label=label)
        E = np.array(E1)
        Z.append(np.array(dos11))

    return z_coords, np.abs(Z).T, E


def get_hartree_pot_z(label='siesta'):

    # temp. input file for macrove
    file_INP = open('macroave.in', 'w')
    file_INP.write('Siesta\n')     # Which code have you used to get the input data?
    file_INP.write('Potential\n')  # Which is the input data used to compute the band offset?
    file_INP.write('%s\n' % label) # Name of the file where the input data is stored
    file_INP.write('2 \n')         # Number of convolutions required to calculate the macro. ave.
    file_INP.write('0 \n')         # First length for the filter function in macroscopic average
    file_INP.write('0 \n')         # Second length for the filter function in macroscopic average
    file_INP.write('0 \n')         # Total charge
    file_INP.write('spline \n')    # Type of interpolation
    file_INP.close()

    # run rho2xsf
    from NanoCore.env import siesta_util_location as sul
    from NanoCore.env import siesta_util_vh as sv
    os.system('%s/%s < macroave.in' % (sul, sv))
    os.system('rm macroave.in')

    # read PDOS
    lines = open('%s.PAV' % label).readlines()
    energy = []; pot = []
                                                                      
    for line in lines:
        if not line.startswith('#'):
            tmp = line.split()
            e = float(tmp[0]); p = float(tmp[1])
            energy.append(e); pot.append(p)

    return energy, pot


def get_total_energy(output_file='stdout.txt'):

    # from standard output file
    os.system("grep 'siesta:         Total =' %s > OUT" % output_file)
    lines = open('OUT').readlines()
    e = float(lines[0].split()[-1])
    os.system('rm OUT')
    return e


#
# OLD SIESTA UTILS
#
bohr2ang = 1./ang2bohr

def read_struct(file_name):
    vec_block = []; atoms_block = []; abc_cell_block = []
    atoms_length = 0; species = []
    n_of_species = 0; name = ''; atoms = []; cell = []; cell_scale = ''
    lattice_constant = 0.
    _is_ang_scale = 0; _is_bohr_scale = 0; _is_scaled_ang_scale = 0
    _is_fraction_scale = 0

    f = open(file_name)
    lines = f.readlines()

    i = 0
    for line in lines:
        #print i
        
        line_s = line.split(); keyword = ''

        if line_s:
            keyword = line_s[0].lower()
            #print keyword

        if keyword == "systemlabel":
            name = line_s[1]

        elif keyword == "latticeconstant":
            lattice_constant = float(line_s[1])
            #print lattice_constant
            try:
                cell_scale = line_s[2]
            except:
                cell_scale = 'Ang'

        elif keyword == "atomiccoordinatesformat":
            if line_s[1].lower() == 'ang':
                _is_ang_scale = 1
            elif line_s[1].lower() == 'bohr':
                _is_bohr_scale = 1
            elif line_s[1].lower() == 'scaledcartesian':
                #print "ON"
                _is_scaled_ang_scale = 1
            elif line_s[1].lower() == 'fractional':
                _is_fraction_scale = 1
            else:
                #print 'Warning : Default atomic scale, "Ang".\n'
                pass

        elif keyword == "numberofatoms":
            atoms_length = int(line_s[1])
            #print "natms", atoms_length

        elif keyword == "numberofspecies":
            n_of_species = int(line_s[1])
            #print "nspec", n_of_species

        elif keyword =="%block":
            keyword_ = line_s[1].lower()
            #print keyword_
            
            if keyword_ == "latticeparameters":
                abc_cell_block = lines[i+1].split()

            elif keyword_ == "latticevectors":
                vec_block = lines[i+1:i+4]

            elif keyword_ == "atomiccoordinatesandatomicspecies":
                atoms_block = lines[i+1:i+1+atoms_length]
                #print "atoms_block", atoms_block

            elif keyword_ == "chemicalspecieslabel":
                temp = lines[i+1:i+1+n_of_species]
                for spec in temp:
                    species.append(spec.split()[2])
                #print species
        i +=1

    # cell converting
    va = 0; vb = 0; vc = 0
    if (not abc_cell_block) and vec_block:
        a1, a2, a3 = vec_block[0].split()
        a1 = float(a1); a2 = float(a2); a3 = float(a3)
        b1, b2, b3 = vec_block[1].split()
        b1 = float(b1); b2 = float(b2); b3 = float(b3)
        c1, c2, c3 = vec_block[2].split()
        c1 = float(c1); c2 = float(c2); c3 = float(c3)
        va = np.array([a1, a2, a3])
        vb = np.array([b1, b2, b3])
        vc = np.array([c1, c2, c3])
        if cell_scale == 'Ang':
            va = lattice_constant * va
            vb = lattice_constant * vb
            vc = lattice_constant * vc
        elif cell_scale == 'Bohr':
            va = lattice_constant * bohr2ang * va
            vb = lattice_constant * bohr2ang * vb
            vc = lattice_constant * bohr2ang * vc
        else:
            #print "Can`t find cell scale"
            pass

        #a, b, c, alpha, beta, gamma = convert_xyz2abc(va, vb, vc)
        cell = np.array([va,vb,vc])

    elif abc_cell_block and (not vec_block):
        a, b, c, alpha, beta, gamma = abc_cell_block.split()
        a = float(a); b = float(b); c = float(c)
        alpha = float(alpha); beta = float(beta); gamma = float(gamma)
        cell = [a, b, c, alpha, beta, gamma]

    # atoms
    for atm in atoms_block:

        if len(atm.split()) == 4:
            x, y, z, spec = atm.split()
        elif len(atm.split()) == 5:
            x, y, z, spec, serial = atm.split()

        x = float(x); y = float(y); z = float(z); spec = int(spec)

        if _is_ang_scale:
            pass
        elif _is_bohr_scale:
            x = bohr2ang * x; y = bohr2ang * y; z = bohr2ang * z

        elif _is_scaled_ang_scale:
            #if vec_cell:
            x = lattice_constant*x
            y = lattice_constant*y
            z = lattice_constant*z
            #elif not vec_cell:
            #    print "Can`t guess cell scale and type\n"
#        elif _is_fraction_scale:
        
        group = re.findall('\d', species[spec-1])
        sym = ''.join(re.findall('[a-zA-Z]', species[spec-1]))
        if group:
            group = int(''.join(group))
        else:
            group = None
        
        atom = Atom(sym, (x, y, z), groupid=group)
        atoms.append(atom)

    if cell.shape == (3,3):
        #XYZ.write_xyz(file_name.replace('fdf','xyz'), atoms, cell)
        return AtomsSystem(atoms, cell=cell)
    else:
        #XYZ.write_xyz(file_name.replace('fdf','xyz'), atoms)
        return AtomsSystem(atoms, cell=None)


def read_struct_out(file_name):
    f = open(file_name)
    lines = f.readlines()
    v1 = Vector(float(lines[0].split()[0]),float(lines[0].split()[1]),float(lines[0].split()[2]))
    v2 = Vector(float(lines[1].split()[0]),float(lines[1].split()[1]),float(lines[1].split()[2]))
    v3 = Vector(float(lines[2].split()[0]),float(lines[2].split()[1]),float(lines[2].split()[2]))
    num_at = int(lines[3].split()[0])
    atoms = []
    for line in lines[4:num_at+4]:
        spec, atn, sx, sy, sz = line.split()
        sx, sy, sz = float(sx), float(sy), float(sz)
        symb = atomic_symbol[int(atn)]
        position = sx*v1 + sy*v2 + sz*v3
        atoms.append(Atom(symb, position))
    return AtomsSystem(atoms, cell = [v1,v2,v3])


def get_eos(pattern='*', struct_file='STRUCT.fdf'):
    # should be replaced by something using xml parser...
    dirs = glob(pattern)
    dirs.sort()

    volume = []
    for f in dirs:
        os.chdir(f)
        os.chdir('input') #
        atoms = read_fdf("%s" % struct_file)
        cell = atoms.get_cell()
        v = Vector(cell[0]).dot( Vector(cell[1]).cross(Vector(cell[2])) )
        volume.append(v)
        os.chdir('..')
        os.chdir('..') #
    os.system("grep 'siesta:         Total =' */stdout.txt > OUT")
    lines = open('OUT').readlines()
    volume = np.array(volume)

    energy = []
    for line in lines:
        e = float(line.split()[-1])
        energy.append(e)
    energy = np.array(energy)

    import pylab as plb # this includes numpy as np!
    from scipy.optimize import leastsq

    # make a vector to evaluate fits on with a lot of points so it looks smooth                    
    vfit = np.linspace(min(volume),max(volume),100)
 
    ### fit a parabola to the data
    # y = ax^2 + bx + c
    a,b,c = plb.polyfit(volume, energy, 2) #this is from pylab
 
    # now here are our initial guesses.
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4

    # now we have to create the equation of state function
    def Murnaghan(parameters,vol):
        '''
        given a vector of parameters and volumes, return a vector of energies.
        equation From PRB 28,5480 (1983)
        '''
        E0 = parameters[0]
        B0 = parameters[1]
        BP = parameters[2]
        V0 = parameters[3]
        E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)
        return E
 
    # and we define an objective function that will be minimized
    def objective(pars,y,x):
        # we will minimize this function
        err =  y - Murnaghan(pars,x)
        return err

    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function

    murnpars, ier = leastsq(objective, x0, args=(energy, volume)) #this is from scipy

    # now we make a figure summarizing the results
    plb.plot(volume,energy,'ro')
    plb.plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
    plb.plot(vfit, Murnaghan(murnpars,vfit), label='Murnaghan fit')
    plb.xlabel('Volume ($\AA^3$)')
    plb.ylabel('Energy (eV)')
    plb.legend(loc='best')

    # add some text to the figure in figure coordinates
    ax = plb.gca()
    plb.text(0.4,0.5,'Min volume = %1.2f $\AA^3$' % murnpars[3],
             transform = ax.transAxes)
    plb.text(0.4,0.4,'Bulk modulus = %1.2f eV/$\AA^3$ = %1.2f GPa' % (murnpars[1],
                                                                      murnpars[1]*160.21773)
             , transform = ax.transAxes)
    plb.savefig('a-eos.png')
    plb.show()

    print ('initial guesses  : ',x0)
    print ('fitted parameters: ', murnpars)

from ase.calculators.vasp import vasp


