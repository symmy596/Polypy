import numpy as np


def open_file(file_name):
    '''Open specificed 'file_name'.  If existant return 
       list=[fp, file_name]'''
    fp = open(file_name,'r')
    this_list = []
    this_list.append(fp)
    this_list.append(file_name)
    return this_list

def fetch_line_as_tmp_string(fp):
    '''Returns next line in fp'''
    tmp_str = fp.readline()
    this_str = tmp_str.rstrip('\n')
    return this_str

def fetch_line_as_lst_tmp_strings(fp):
    '''Returns next line in fp as list of tmp_strings'''
    tmp_str = fp.readline()
    this_data = tmp_str.split()
    return this_data

def fetch_line_as_ints(fp):
    '''Returns next line as np array of ints'''
    
    tmp_str = fp.readline()
    data = tmp_str.split()
    this_ints = []
    for value in data:
        try:
            this_ints.append( int(value) )
        except:
            "Unable to convert: "+value+" from: "+tmp_str+"into int"
    np_ints=np.array(this_ints)
    return np_ints


def fetch_line_as_floats(fp):
    '''Returns next line as np array of floats'''
    
    tmp_str = fp.readline()
    this_data = tmp_str.split()
    this_floats = []
    for value in this_data:
        try:
            this_floats.append( float(value) )
        except:
            "Unable to convert: "+value+" from: "+tmp_str+"into float"
    np_floats=np.array(this_floats)
    return np_floats


def check_line_begins_with(fp, require):
    '''Check that next line in fp begins with 'require'
    Return rest of line as tmp_string'''

    tmp_str = fp.readline()
    this_data = tmp_str.split()

    key_word = this_data.pop(0)

    ErrStr = 'In file: read tmp_string '+key_word+' is not key word '+require

    if( key_word.lower() == require ):
        return this_data
    else:
        raise ValueError(ErrStr)


def check_tmp_string(tmp_str, lst_require):
    '''Check whether tmp_string is in list of tmp_strings'''

    ErrStr = 'Read tmp_string '+tmp_str+' is not key word '
    for require in lst_require:
        ErrStr=ErrStr+require+','

    if( tmp_str.lower() in lst_require ):
        return
    else:
        raise ValueError(ErrStr)

def read_config_title(fp):
    '''Returns title from open CONFIG'''
    title = fetch_line_as_tmp_string(fp)
    return title


def read_config_style(fp):
    '''Returns np_style from open CONFIG
       np_style[0] is style_coor: direct(=0) or fractional(=1)
       np_style[1] is style_cell: cubic(=0) orthogonal(=1) other(>=2)'''
    np_style = fetch_line_as_ints(fp)
    return np_style


def read_config_lvs(fp):
    '''Returns lattice vectors as np.array((3,3)) from open CONFIG
      If not able to convert to np.array raise exception and return 'Not valid DLMONTE CONFIG lattice vectors' '''

    this_vectors = []
    for i in range (0, 3):
        this_vectors.append(fetch_line_as_floats(fp) )
    np_lvs = np.array(this_vectors)
    return np_lvs
