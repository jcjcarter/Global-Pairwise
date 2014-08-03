### You may paste any/all of the code in this file that you find
### useful into your solutions.

import random

def permute_string(s):
    """
    Return a new string with the characters in s randomly permuted.

    Arguments:
    s -- string

    Returns:
    Random permutation of s
    """
    charlist = list(s)
    random.shuffle(charlist)
    newstr = "".join(charlist)
    return newstr

def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    with open(filename) as f:
        p = f.read()
    p = p.rstrip()
    return p

def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    M = {}
    with open(filename) as f:
        ykeys = f.readline()
        ykeychars = ykeys.split()
        for line in f.readlines():
            vals = line.split()
            xkey = vals.pop(0)
            M[xkey] = {}
            for ykey, val in zip(ykeychars, vals):
                M[xkey][ykey] = int(val)
    return M
