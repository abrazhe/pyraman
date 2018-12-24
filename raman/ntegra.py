## Reading images and spectra exported as txt from Ntegra system

import numpy as np

def parse_txt_line(line):
    splitted = line.strip().split(',')
    floats = [np.float(el) for el in splitted]
    x = floats[0]
    y = floats[1]
    nu = np.array(floats[2::2])
    v = np.array(floats[3::2])
    return x,y,nu,v
    
def load_image(fname):
    "Loads image exported as text from Ntegra software (fixme)"
    with open(fname, 'r') as file_handle:
        lines = file_handle.readlines()
    content = [parse_txt_line(l) for l in lines]
    all_x = [c[0] for c in content]
    all_y = [c[1] for c in content]
    
    uni_x = np.unique(all_x)
    uni_y = np.unique(all_y)
    img_shape = (len(uni_x),len(uni_y))
        
    nu = content[0][2]
    knu = np.argsort(nu)
    all_spectra = np.array([c[3][knu] for c in content])
    
    # todo: re-formulate as reshape
    out = np.zeros((img_shape)+(len(nu),))
    k = 0
    for r in range(img_shape[0]):
        for c in range(img_shape[1]):
            out[r,c] = all_spectra[k]
            k +=1
    
    return out, nu[knu], uni_x, uni_y
