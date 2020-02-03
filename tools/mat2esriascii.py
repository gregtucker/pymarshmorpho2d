#!/bin/python
"""mat2esriascii.py: read data from a matlab .mat file, write any 2D arrays to
an esri ascii format file.

Examples
--------
>>> from scipy.io import savemat
>>> import numpy as np
>>> a = np.array([[1, 2, 3], [4, 5, 6]])
>>> mdict = {'__version__': 1.0, '__globals__':[], 'mydata':a}
>>> savemat('test.mat', mdict)
>>> main(['me', 'test.mat'])
Wrote array to file "mydata.asc" with 2 rows and 3 columns.
>>> print(open('mydata.asc').read())
NCOLS 3
NROWS 2
XLLCORNER 0
YLLCORNER 0
CELLSIZE 1
NODATA_VALUE -9999
1.000000 2.000000 3.000000
4.000000 5.000000 6.000000
"""

import sys
import numpy
from scipy.io import loadmat


def parse_header_info(args):
    header = {}
    if len(args) > 2:
        header['cellsize'] = args[2]
    else:
        header['cellsize'] = '1'
    if len(args) > 3:
        header['nodata_value'] = args[3]
    else:
        header['nodata_value'] = '-9999'
    if len(args) > 5:
        header['xllc'] = args[4]
        header['yllc'] = args[5]
    else:
        header['xllc'] = '0'
        header['yllc'] = '0'
    return header

def make_header_string(header, myarray):
    header_string = ""
    header_string += 'NCOLS ' + str(myarray.shape[1]) + '\n'
    header_string += 'NROWS ' + str(myarray.shape[0]) + '\n'
    header_string += 'XLLCORNER ' + header['xllc'] + '\n'
    header_string += 'YLLCORNER ' + header['yllc'] + '\n'
    header_string += 'CELLSIZE ' + header['cellsize'] + '\n'
    header_string += 'NODATA_VALUE ' + header['nodata_value']
    return header_string

def write_to_file(basename, myarray, header):
    header_string = make_header_string(header, myarray)
    numpy.savetxt(basename + '.asc', myarray, header=header_string,
                  comments='', fmt='%f')
    print('Wrote array to file "' + basename + '.asc" with '
          + str(myarray.shape[0]) + ' rows and '
          + str(myarray.shape[1]) + ' columns.')

def find_matrices_and_write_them_to_file(matdict, header):
    for item in matdict:
        a = matdict[item]
        if type(a) is numpy.ndarray and a.ndim == 2:
            write_to_file(item, a, header)

def main(args):
    filename = args[1]
    header = parse_header_info(args)
    d = loadmat(filename)  # read dict from .mat file
    find_matrices_and_write_them_to_file(d, header)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python mat2esriascii.py <filename> {cellsize} {nodata} {xllc} {yllc}')
    else:
        main(sys.argv)
