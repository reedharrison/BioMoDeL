###################################################################################################################################
# Class:        Universe
# Purpose:      Class to pass for analysis. Contains information for analysis
# Arguments:   
###################################################################################################################################
class structData:
    # import prody as prody

    def __init__(self, pdb, pqr=None, selstr=None):
        self.pdb = pdb
        self.pqr = pqr
        self.selstr = selstr

    def setPDB(self, newpdb):
        self.pdb = newpdb
    def setPQR(self, newpqr):
        self.pqr = newpqr
    def setSel(self, newsel):
        self.selstr = newsel

    def getPDB(self):
        return self.pdb
    def getPQR(self):
        return self.pqr
    def getSel(self):
        return self.selstr

###################################################################################################################################
# Function:     calc_features
# Purpose:      Calculate features from available analyses
# Arguments:   
###################################################################################################################################
def calc_feature(analysis_files, feature, cutoff=4, path_dssp='C:\Python27\Scripts\dssp-2.0.4-win32.exe'):
    import prody as prody
    import numpy as np
    import os as os
    # import scipy.spatial.distance as sci

    prevdir = os.getcwd()

    pdb = analysis_files.getPDB()
    pqr = analysis_files.getPQR()
    selstr = analysis_files.getSel()

    print str(pdb)+': '+feature
    container = np.zeros(0)

    if feature is 'del_sasa': # SASA lost upon binding, requires 2 chain selections. Otherwise use 'sasa' to calculate each selection
        val_1 = calc_sasa(pdb, sel=selstr[0], per_res=False, path_dssp=path_dssp)
        val_2 = calc_sasa(pdb, sel=selstr[-1], per_res=False, path_dssp=path_dssp)
        val_12 = calc_sasa(pdb, sel=selstr[0]+' or '+selstr[-1], per_res=False, path_dssp=path_dssp)
        container = np.append(container, val_12 - val_1 - val_2)
    elif feature is 'sasa': # SASA for each element of selstr
        for phrase in selstr:
            container = np.append(container, calc_sasa(pdb, sel=phrase, per_res=False, path_dssp=path_dssp))

    elif feature is 'num_contact': # Number of contacts between first and last element of selstr
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2)
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1)
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_acidic': # Number of contacts between first and last element of selstr containing acidic AA
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).acidic
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).acidic
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_basic': # Number of contacts between first and last element of selstr containing basic AA
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).basic
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).basic
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_ionize': # Number of contacts between first and last element of selstr containing ionizable AA
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).select('acidic or basic')
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).select('acidic or basic')
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_hbond': # Number of hbonds between first and last element of selstr
        container = np.append(container, find_int(pdb, selstr[0], selstr[-1], mode='hbond'))
    elif feature is 'num_contact_sb': # Number of saltbridges between first and last element of selstr
        container = np.append(container, find_int(pdb, selstr[0], selstr[-1], mode='sb'))
    elif feature is 'num_contact_aliphatic': # Number of aliphatic contacts between first and last element of selstr
        container = np.append(container, find_int(pdb, selstr[0], selstr[-1], mode='aliphatic'))
    elif feature is 'num_contact_aromatic': # Number of aromatic contacts between first and last element of selstr
        container = np.append(container, find_int(pdb, selstr[0], selstr[-1], mode='aromatic'))
    elif feature is 'num_contact_arg':
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).select('resname ARG')
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).select('resname ARG')
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_lys':
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).select('resname LYS')
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).select('resname LYS')
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_his':
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).select('resname HIS')
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).select('resname HIS')
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_trp':
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).select('resname TRP')
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).select('resname TRP')
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_polar':
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).polar
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).polar
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)
    elif feature is 'num_contact_nonpolar':
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).hydrophobic
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).hydrophobic
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1 + val2)

    elif feature is 'num_acidic': # Number of acidic contacts, per selection, in between first and last selections
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).acidic
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).acidic
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1)
        container = np.append(container, val2)
    elif feature is 'num_basic': # Number of basic contacts, per selection, in between first and last selections
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).basic
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).basic
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1)
        container = np.append(container, val2)
    elif feature is 'num_ionize': # Number of ionizable contacts, per selection, in between first and last selections
        pdb1 = pdb.select(selstr[0])
        pdb2 = pdb.select(selstr[-1])
        try:
            con1 = prody.Contacts(pdb1).select(cutoff, pdb2).select('acidic or basic')
        except:
            con1 = None
        try:
            con2 = prody.Contacts(pdb2).select(cutoff, pdb1).select('acidic or basic')
        except:
            con2 = None
        val1 = np.zeros(0)
        val2 = np.zeros(0)
        if con1 is None: val1 = 0
        else: val1 = len(con1)
        if con2 is None: val2 = 0
        else: val2 = len(con2)
        container = np.append(container, val1)
        container = np.append(container, val2)
    elif feature is 'num_hbond': # Number of intra-selection hbonds
        for phrase in selstr:
            container = np.append(container, find_int(pdb, phrase, mode='hbond'))
    elif feature is 'num_sb': # Number of intra-selection salt bridges
        for phrase in selstr:
            container = np.append(container, find_int(pdb, phrase, mode='sb'))
    elif feature is 'num_aliphatic': # Number of aliphatic contacts in each selection
        for phrase in selstr:
            container = np.append(container, find_int(pdb, phrase, mode='aliphatic'))
    elif feature is 'num_aromatic': # Number of aromatic contacts in each selection
        for phrase in selstr:
            container = np.append(container, find_int(pdb, phrase, mode='aromatic'))

    elif feature is 'num_res_aliphatic': # Number of aliphatic residues in each selection
        for phrase in selstr:
            container = np.append(container, len(np.unique(pdb.select(phrase).aliphatic.getResnums())))
    elif feature is 'num_res_aromatic': # Number of aromatic residues in each selection
        for phrase in selstr:
            container = np.append(container, len(np.unique(pdb.select(phrase).aromatic.getResnums())))
    elif feature is 'num_res_acidic': # Number of acidic residues in each selection
        for phrase in selstr:
            container = np.append(container, len(np.unique(pdb.select(phrase).acidic.getResnums())))
    elif feature is 'num_res_basic': # Number of basic residues in each selection
        for phrase in selstr:
            container = np.append(container, len(np.unique(pdb.select(phrase).basic.getResnums())))
    elif feature is 'num_res_ionize': # Number of ionizable amino acids in each selection
        for phrase in selstr:
            container = np.append(container, len(np.unique(pdb.select(phrase).select('acidic or basic').getResnums())))
    
    elif feature is 'mass': # Mass of each selection
        mass_dict = {'C':12,'N':14,'S':32,'O':16,'H':1}
        for phrase in selstr:
            val = [mass_dict.get(atom.getElement()) for atom in pdb.select(phrase).iterAtoms()]
            val = [x for x in val if x is not None]
            container = np.append(container, sum(val))
    elif feature is 'dMass': # Sum of mass differences between first strsel and all other strsel
        mass_dict = {'C':12,'N':14,'S':32,'O':16,'H':1}
        val1 = [mass_dict.get(atom.getElement()) for atom in pdb.select(selstr[0]).iterAtoms()]
        val1 = [x for x in val1 if x is not None]
        counter = 0
        for phrase in selstr:
            val2 = [mass_dict.get(atom.getElement()) for atom in pdb.select(phrase).iterAtoms()]
            val2 = [x for x in val2 if x is not None]
            counter = counter + (sum(val1) - sum(val2))
        container = np.append(container, counter)
    elif feature is 'tMass': # Total mass of all selstr
        mass_dict = {'C':12,'N':14,'S':32,'O':16,'H':1}
        counter = 0
        for phrase in selstr:
            val = [mass_dict.get(atom.getElement()) for atom in pdb.select(phrase).iterAtoms()]
            val = [x for x in val if x is not None]
            counter = counter + sum(val)
        container = np.append(container, counter)

    elif feature is 'charge': # Per-atom charge of each selection
        for phrase in selstr:
            container = np.vstack(container, pqr.select(phrase).getCharges())
    elif feature is 'dCharge': # Sum of charge difference between first strsel and all other strsel
        val1 = sum(pqr.select(selstr[0]).getCharges())
        counter = 0
        for phrase in selstr:
            val2 = sum(pqr.select(phrase).getCharges())
            counter = counter + (val1 - val2)
        container = np.append(container, counter)
    elif feature is 'tCharge': # Sum of charges in all strsel
        for phrase in selstr:
            container = np.append(container, pqr.select(phrase).getCharges())

    elif feature is 'num_contact_elec_attract': # Total number of attractive electrostatic contacts between first and last strsel
        container = np.append(container, find_int(pqr, selstr[0], selstr[-1], mode='coulomb', cutoff=6))
    elif feature is 'num_contact_elec_repulse': # Total number of repulsive electrostatic interactions between first and all other strsel
        container = np.append(container, find_int(pqr, selstr[0], selstr[-1], mode='coulomb_repulse', cutoff=6))
    elif feature is 'num_elec_attract':
        for phrase in selstr:
            container = np.append(container, find_int(pqr, phrase, mode='coulomb', cutoff=6))
    elif feature is 'num_elec_repulse':
        for phrase in selstr:
            container = np.append(container, find_int(pqr, phrase, mode='coulomb_repulse', cutoff=6))

    os.chdir(prevdir)
    # print container
    return container

###################################################################################################################################
# Function:     calc_sasa
# Purpose:      Calculate solvent accessible surface 
# Arguments:   
###################################################################################################################################
def calc_sasa(pdb, sel='protein', per_res=False, path_dssp='C:\Python27\Scripts\dssp-2.0.4-win32.exe'):
    import os as os
    import tempfile as tempfile
    import shutil as shutil
    import prody as prody
    import numpy as np

    dssp_in = 'dssp.in'

    prevdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(os.path.expanduser(tmpdir))

        prody.writePDB(dssp_in+'.pdb', pdb.select(sel))
        pdb_dssp = prody.parsePDB(dssp_in+'.pdb')
        execDSSP(dssp_in+'.pdb', dssp=path_dssp)
        prody.parseDSSP(dssp_in+'.dssp', pdb_dssp)

        if per_res:
            return(pdb_dssp._data['dssp_acc'])
        else:
            return(np.sum(pdb_dssp._data['dssp_acc']))

    os.chdir(prevdir)
    # return pdb

###################################################################################################################################
# Function:     find_int
# Purpose:      To identify a particular kind of biological interaction
# Arguments:    pdb = Prody PDB object with atomic information, coordinate info only used if dcd not supplied
#               dcd = Ensemble or trajectory with coordinate information
#               mode = Type of interaction to find, options are 'sb', 'hbond', 'aliphatic', 'aromatic'
###################################################################################################################################
def find_int(pdb, sel1, sel2=None, dcd=None, mode='hbond', cutoff=4):
    import scipy.spatial.distance as sci
    if sel2 is None:
        sel2 = sel1 
    if dcd is None:
        pdb1 = pdb.select(sel1)
        pdb2 = pdb.select(sel2)
        if mode is 'sb':
            sb1 = sci.cdist(pdb1.basic.getCoords(), pdb2.acidic.getCoords())
            sb2 = sci.cdist(pdb1.acidic.getCoords(), pdb2.basic.getCoords())

            try:
                num_int1 = len(sb1[sb1<cutoff])
            except:
                num_int1 = 0
            try:
                num_int2 = len(sb2[sb2<cutoff])
            except:
                num_int2 = 0
            num_int = num_int1 + num_int2

            return num_int

        elif mode is 'hbond':
            sel_hba = 'name O OH2 OW OD1 OD2 SG OE1 OE1 OE2 ND1 NE2 SD OG OG1 OH'
            sel_hbd = 'name N OH2 OW NE NH1 NH2 ND2 SG NE2 ND1 NZ OG OG1 NE1 OH'
            pdb1_hba = pdb1.select(sel_hba)
            pdb1_hbd = pdb1.select(sel_hbd)
            pdb2_hba = pdb2.select(sel_hba)
            pdb2_hbd = pdb2.select(sel_hbd)

            hb1 = sci.cdist(pdb1_hba.getCoords(), pdb2_hbd.getCoords())
            hb2 = sci.cdist(pdb1_hbd.getCoords(), pdb2_hba.getCoords())

            try:
                num_int1 = len(hb1[hb1<cutoff])
            except:
                num_int1 = 0
            try:
                num_int2 = len(hb2[hb2<cutoff])
            except:
                num_int2 = 0
            num_int = num_int1 + num_int2
            
            return num_int

        elif mode is 'aliphatic':
            try:
                int1 = sci.cdist(pdb1.aliphatic.getCoords(), pdb2.aliphatic.getCoords())
            except:
                print 'No aliphatic interactions'

            try:
                num_int = len(int1<cutoff)
            except:
                num_int = 0

            return num_int

        elif mode is 'aromatic':
            try:
                int1 = sci.cdist(pdb1.aromatic.getCoords(), pdb2.aromatic.getCoords())
            except:
                print 'No aromatic interactions'

            try:
                num_int = len(int1<cutoff)
            except:
                num_int = 0

            return num_int

        elif mode is 'coulomb': # make sure pdb is the result of prody.parsePQR() that contains charge info
            pdb1_coul_pos = pdb1.select('charge>0')
            pdb1_coul_neg = pdb1.select('charge<0')
            pdb2_coul_pos = pdb2.select('charge>0')
            pdb2_coul_neg = pdb2.select('charge<0')

            coul1 = sci.cdist(pdb1_coul_pos.getCoords(), pdb2_coul_neg.getCoords())
            coul2 = sci.cdist(pdb1_coul_neg.getCoords(), pdb2_coul_pos.getCoords())

            try:
                num_int1 = len(coul1[coul1<cutoff])
            except:
                num_int1 = 0
            try:
                num_int2 = len(coul2[coul2<cutoff])
            except:
                num_int2 = 0
            num_int = num_int1 + num_int2
            
            return num_int

        elif mode is 'coulomb_repulse': # make sure pdb is the result of prody.parsePQR() that contains charge info
            pdb1_coul_pos = pdb1.select('charge>0')
            pdb1_coul_neg = pdb1.select('charge<0')
            pdb2_coul_pos = pdb2.select('charge>0')
            pdb2_coul_neg = pdb2.select('charge<0')

            coul1 = sci.cdist(pdb1_coul_pos.getCoords(), pdb2_coul_pos.getCoords())
            coul2 = sci.cdist(pdb1_coul_neg.getCoords(), pdb2_coul_neg.getCoords())

            try:
                num_int1 = len(coul1[coul1<cutoff])
            except:
                num_int1 = 0
            try:
                num_int2 = len(coul2[coul2<cutoff])
            except:
                num_int2 = 0
            num_int = num_int1 + num_int2
            
            return num_int
    

###################################################################################################################################
# Function:     find_int
# Purpose:      To identify a particular kind of biological interaction
# Arguments:    pdb = Prody PDB object with atomic information, coordinate info only used if dcd not supplied
#               dcd = Ensemble or trajectory with coordinate information
#               type = Type of interaction to find, options are 'sb', 'hbond', 'aliphatic', 'aromatic'
###################################################################################################################################
def import_sel_macro():
    import prody

    defSelectionMacro('acceptor', 'name O OH2 OW OD1 OD2 SG OE1 OE1 OE2 ND1 NE2 SD OG OG1 OH')
    defSelectionMacro('donor', 'name N OH2 OW NE NH1 NH2 ND2 SG NE2 ND1 NZ OG OG1 NE1 OH')

###################################################################################################################################
# Function:     execDSSP
# Purpose:      Runs DSSP on given pdb_file
# Arguments:    This function is essentially the same as the one in Prody, except that I allow the path to dssp to be specified
###################################################################################################################################
def execDSSP(pdb, outputname=None, outputdir=None, stderr=True, dssp=None):
    """Execute DSSP for given *pdb*.  *pdb* can be a PDB identifier or a PDB
    file path.  If *pdb* is a compressed file, it will be decompressed using
    Python :mod:`gzip` library.  When no *outputname* is given, output name
    will be :file:`pdb.dssp`.  :file:`.dssp` extension will be appended
    automatically to *outputname*.  If :file:`outputdir` is given, DSSP
    output and uncompressed PDB file will be written into this folder.
    Upon successful execution of :command:`dssp pdb > out` command, output
    filename is returned.  On Linux platforms, when *stderr* is false,
    standard error messages are suppressed, i.e.
    ``dssp pdb > outputname 2> /dev/null``.

    For more information on DSSP see http://swift.cmbi.ru.nl/gv/dssp/.
    If you benefited from DSSP, please consider citing [WK83]_.

    .. [WK83] Kabsch W, Sander C. Dictionary of protein secondary structure:
       pattern recognition of hydrogen-bonded and geometrical features.
       *Biopolymers* **1983** 22:2577-2637."""

    import prody
    import os as os

    if dssp is None:
        raise EnvironmentError('command not found: dssp executable is not '
                               'found in one of system paths')
    assert outputname is None or isinstance(outputname, str),\
        'outputname must be a string'
    assert outputdir is None or isinstance(outputdir, str),\
        'outputdir must be a string'
    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir is None:
            pdb = gunzip(pdb, os.path.splitext(pdb)[0])
        else:
            pdb = gunzip(pdb, os.path.join(outputdir,
                         os.path.split(os.path.splitext(pdb)[0])[1]))
    if outputdir is None:
        outputdir = '.'
    if outputname is None:
        out = os.path.join(outputdir,
                           os.path.splitext(os.path.split(pdb)[1])[0] +
                           '.dssp')
    else:
        out = os.path.join(outputdir, outputname + '.dssp')

    if not stderr and PLATFORM != 'Windows':
        status = os.system('{0} {1} > {2} 2> /dev/null'.format(
                           dssp, pdb, out))
    else:
        status = os.system('{0} {1} > {2}'.format(dssp, pdb, out))

    if status == 0:
        return out

###################################################################################################################################
# Function:     execPDB2PQR
# Purpose:      Run PDB2PQR on supplied pdbfile
# Arguments:    
###################################################################################################################################
def execPDB2PQR(path_pdb2pqr, pdbfile, optargs='--ff charmm --chain', outfile=None):
    import os as os
    if outfile is None:
        outfile = os.path.splitext(pdbfile)[0]+'.pqr'
    os.system('{0} {1} {2} {3}'.format(path_pdb2pqr, optargs, pdbfile, outfile))
    return outfile

###################################################################################################################################
# Function:     get_pkas
# Purpose:      Calculate pkas with PROPKA
# Arguments:    
###################################################################################################################################
def get_pkas(pdbfile, optargs=None):
    """Run a single PROPKA calculation using *pdbfile* as input.
    Commandline options can be passed as a **list** in *optargs*.
    .. rubric:: Example
    ::
       single("protein.pdb", optargs=["--mutation=N25R/N181D", "-v", "--pH=7.2"])
    """
    import propka.lib, propka.molecular_container
    import numpy as np

    optargs = optargs if optargs is not None else []
    options, ignored_pdbfiles = propka.lib.loadOptions(*optargs)

    # Run PROPKA
    my_molecule = propka.molecular_container.Molecular_container(pdbfile, options)
    my_molecule.calculate_pka()

    # Extract pKa calculations
    summary = propka.output.getSummarySection(my_molecule, conformation = 'AVR', parameters=my_molecule.version.parameters)
    summary = summary.split('\n')[3:-1]
    resname = []
    resnum = np.zeros(0)
    chain = []
    pka = np.zeros(0)
    model_pka = np.zeros(0)
    for line in summary:
        resname.append(line.split()[0])
        resnum = np.append(resnum, int(line.split()[1]))
        chain.append(line.split()[2])
        pka = np.append(pka, float(line.split()[3]))
        model_pka = np.append(model_pka, float(line.split()[4]))

    # Return results
    return(resname, resnum, pka, model_pka)
    # my_molecule.write_pka()
    #return my_molecule

###################################################################################################################################
# Function:     modelLoops
# Purpose:      Model missing atoms in PDB file
# Arguments:    
###################################################################################################################################
def modelLoops(pdbid, chids, alnfile='temp.ali'):
    from modeller import environ, model, alignment
    from modeller.automodel import loopmodel
    import os as os
    import tempfile as tempfile
    import shutil as shutil
    import prody as prody
    import numpy as np

    prevdir = os.getcwd()
    pdb = []
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(os.path.expanduser(tmpdir))

        prody.fetchPDB(pdbid)
        e = environ()
        for chid in chids:
            knowns = pdbid+'_'+chid
            sequence = pdbid+'_'+chid+'_full'

            try:    # Try to model structure
                aln = alignment(e)
                m = model(e, file=pdbid, model_segment=('FIRST:'+chid, 'LAST:'+chid))
                aln.append_model(m, atom_files=pdbid, align_codes=knowns)
                aln.append_sequence(getSeqres(pdbid, chid)[0])
                aln[-1].code = sequence

                aln.malign()
                aln.write(file=alnfile, alignment_format='PIR')

                a = loopmodel(e, alnfile=alnfile, knowns=knowns, sequence=sequence)
                a.make()
                pdbfile = a.outputs[0]['name']

                h = model(e, file=pdbfile)
                aln = alignment(e)
                aln.append_model(m, atom_files=pdbid, align_codes=knowns)
                aln.append_model(h, atom_files=pdbfile, align_codes=sequence)
                h.res_num_from(m, aln)  # Restore old residue numbering and chain indexing
                h.write(file=pdbfile)

                if not pdb:
                    pdb = prody.parsePDB(pdbfile)
                else:
                    pdb = pdb + prody.parsePDB(pdbfile)

            except: # If it fails, return original PDB file. Likely the original file has no gaps (i.e. NOT EFFICIENT).
                print 'PDB %s chain %s could not be modelled' %(pdbid, chid)
                ref = prody.parsePDB(pdbid)
                sel = 'chain %s' %chid
                atom = ref.select(sel)
                reffile = knowns+'.pdb'
                prody.writePDB(reffile, atom)

                if not pdb:
                    pdb = prody.parsePDB(reffile)
                else:
                    pdb = pdb + prody.parsePDB(reffile)

    os.chdir(prevdir)
    return pdb

###################################################################################################################################
# Function:     findGaps
# Purpose:      Find residues missing from a PDB file by parsing the PDB header
# Arguments:    
###################################################################################################################################

###################################################################################################################################
# Function:     getSeqres
# Purpose:      Get sequence for biomolecule from original PDB file
# Arguments:    
###################################################################################################################################
def getSeqres(pdbid, chids):
    import prody
    sequences = []
    polymers = prody.parsePDBHeader(pdbid, 'polymers')
    for polymer in polymers:
        if polymer.chid in chids:
            sequences.append(polymer.sequence)
    return sequences