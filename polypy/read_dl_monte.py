import numpy as np
from polypy import read_utils as rd_ut


def check_tmp_string(tmp_str, lst_require):
    '''Check whether tmp_string is in list of tmp_strings'''

    ErrStr = 'Read tmp_string '+tmp_str+' is not key word '
    for require in lst_require:
        ErrStr=ErrStr+require+','

    if( tmp_str.lower() in lst_require ):
        return
    else:
        raise ValueError(ErrStr)


def read_config_sysmols(fp):
    '''Returns lattice vectors as np.array((3,3)) from open CONFIG
      If not able to convert to np.array raise exception and return 'Not valid DLMONTE CONFIG lattice vectors' '''

    require = 'nummol'

    data = check_line_begins_with(fp, require)

    np_nummols = np.array( int( data.pop(0) ) )
    np_maxmols = np.array( data, dtype = np.int64 )
    np_moltypes = np_maxmols.size

    
    return [ np_nummols, np_moltypes, np_maxmols ]


def read_atom(fp):
    '''Read atom name, check type and read coordinates from CONFIG'''
    dict_atom = {}

    lst_tmp_str = fetch_line_as_lst_tmp_strings(fp)

    require = ['core', 'c']

    check_tmp_string(lst_tmp_str[1],require)

    dict_atom['label'] = lst_tmp_str[0]
    dict_atom['type'] = lst_tmp_str[1].lower()

    dict_atom['coor'] = fetch_line_as_floats(fp)

    return dict_atom


def read_molecule(fp):
    '''Returns molecule as dictionary from open CONFIG'''
    
    dict_molecule = {}

    require = 'molecule'

    data = check_line_begins_with(fp, require)

    dict_molecule['name'] = data.pop(0)
    dict_molecule['numatoms'] = np.array( int( data[0] ) )
    dict_molecule['maxatoms'] = np.array( int( data[1] ) )

    list_atoms = []

    for iatom in range( dict_molecule['numatoms'] ):
        list_atoms.append( read_atom(fp) )

    dict_molecule['atoms'] = list_atoms

    return dict_molecule


def read_dlmonte_config(fp):
    '''Returns entire config as dictionary from open CONFIG'''

    dict_config = {}

    try:
        dict_config['title'] = rd_ut.read_config_title(fp)

        dict_config['style'] = rd_ut.read_config_style(fp)

        dict_config['lvs'] = rd_ut.read_config_lvs(fp)

        [ dict_config['nummols'], dict_config['moltypes'], dict_config['maxmols'] ] = read_config_sysmols(fp)

        list_molecules = []

        for imol in range( dict_config['nummols'] ):
            list_molecules.append( read_molecule(fp) )

        dict_config['mols'] = list_molecules

    except:
        return "Error"

    return dict_config

def open_config(config_file):
    '''Open and read in an entire CONFIG file
       return as dictionary'''

    [fp,tmp_str] = open_file(config_file)
    
    dict_config = read_config(fp)

    fp.close()

    if( dict_config != "Error"):
        return dict_config
    else:
        print("Unable to read config")

def read_trajectory(archive_file):
    '''Open and read in an entire ARCHIVE file
    return trajectory as dictionary'''

    [fp,tmp_str] = open_file(archive_file)
    
    iconfig = 0

    dict_traj = {}

    while True:
        config = read_config(fp)
        if(config != "Error"):
            dict_traj[iconfig] = config
            iconfig += 1
        else:
            break

    fp.close()

    if( iconfig != 0 ):
        dict_traj['numconfigs'] = iconfig
        print("Read trajectory with "+str(iconfig)+" configurations")
        return dict_traj
    else:
        print("Unable to read trajectory")

