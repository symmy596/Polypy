import numpy as np
from polypy import read_utils as rd_ut

def read_atom(fp, output_style):
    '''Read atom name, check type and read coordinates from CONFIG'''
    dict_atom = {}

    try:

        lst_tmp_str = rd_ut.fetch_line_as_lst_tmp_strings(fp)

        dict_atom['label'] = lst_tmp_str[0]
        dict_atom['id'] = int( lst_tmp_str[1] )
        dict_atom['coor'] = rd_ut.fetch_line_as_floats(fp)

        if( output_style > 0):
            dict_atom['vel'] = rd_ut.fetch_line_as_floats(fp)

        if( output_style > 1):
            dict_atom['force'] = rd_ut.fetch_line_as_floats(fp)

    except:
         return("Error")

    return dict_atom

def read_config_atoms(fp,style):
    '''Read numatoms atoms from open CONFIG
    or read from file until there are no more atoms
    Return atoms as dictionary of labels and coordinates
    Currently presume that this is EOF
    i.e. to read trajectory, numatoms must be passed'''

    iatom = 0
    atoms = {}

    if( len(style) == 2 ):
        while True:
            atom = read_atom( fp, style[0] )
            if(atom != "Error"):
                atoms[iatom] = atom
                iatom += 1
            else:
                 break
    else:
        for iatom in range( style[2] ):
            atom = read_atom(fp,style[0])
            if(atom != "Error"):
                atoms[iatom] = atom
                iatom += 1
            else:
                print("Unable to read "+str(style[2])+" atoms from CONFIG")
                return "Error" 

    atoms['numatoms'] = iatom

    return atoms


def read_dlpoly_config(fp):
    '''Returns entire config as dictionary from open CONFIG'''

    dict_config = {}

    try:
        dict_config['title'] = rd_ut.read_config_title(fp)
        dict_config['style'] = rd_ut.read_config_style(fp)
        dict_config['lvs'] = rd_ut.read_config_lvs(fp)
        dict_config['atoms'] = read_config_atoms( fp, dict_config['style'] ) 

        if ( dict_config['atoms'] == "Error" ):
            return "Error"
            
    except:
        return "Error"

    return dict_config

def open_config(config_file):
    '''Open and read in an entire CONFIG file
       return as dictionary'''

    [fp,tmp_str] = rd_ut.open_file(config_file)
    dict_config = read_dlpoly_config(fp)

    fp.close()

    if( dict_config != "Error"):
        return dict_config
    else:
        print("Unable to read config")

def read_traj_timestep(fp):
    '''Reads timestep keyword and returns line from open HISTORY'''

    require = "timestep"
    data = rd_ut.check_line_begins_with(fp, require)

    return data

def read_traj_config(fp, style):
    '''Returns entire config as dictionary from open HISTORY'''

    dict_config = {}

    try:
        dict_config['timestep'] = read_traj_timestep(fp)
        dict_config['lvs'] = rd_ut.read_config_lvs(fp)
        dict_config['atoms'] = read_config_atoms(fp,style) 

        if ( dict_config['atoms'] == "Error" ):
            return "Error"
            
    except:
        return "Error"

    return dict_config

def read_trajectory(archive_file):
    '''Open and read in an entire ARCHIVE file
    return trajectory as dictionary'''

    [fp,tmp_str] = rd_ut.open_file(archive_file)

    dict_traj = {}
    iconfig = 0

    dict_traj['title'] = rd_ut.read_config_title(fp)
    dict_traj['style'] = rd_ut.read_config_style(fp)

    while True:
        config = read_traj_config( fp, dict_traj['style'] )
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

