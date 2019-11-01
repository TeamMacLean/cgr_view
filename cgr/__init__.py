#!/usr/bin/env python

import os
import subprocess
import sys
import math
import numpy
from scipy.sparse import dok_matrix
import scipy.misc
import re
import matplotlib.pyplot as plt
from PIL import Image


def estimate_genome_size(fasta):
    """guesses genome size from fasta file size, assumes 1 byte ~= 1 base"""
    return(os.path.getsize(fasta))
    
def run_jellyfish(fasta, k, out):
    """runs jellyfish on fasta using k kmer size, produces db file out"""
    genome_size = estimate_genome_size(fasta)
    cmd = ["jellyfish", "count", "-m", str(k), "-s", str(genome_size), fasta, "-o", out] 
    result = subprocess.run(cmd)
    return(result)
    
def get_kmer_list(jellyfish):
    """runs jellyfish dump on jellyfish. captures output as a generator stream. each item = [kmer, count]"""
    cmd = ["jellyfish", "dump", "-c", jellyfish]
    proc = subprocess.Popen(cmd, stdout = subprocess.PIPE )
    for line in proc.stdout:
        yield line.decode("utf-8").rstrip().split(" ")
    proc.wait()
    proc.stdout.close()

def get_grid_size(k):
    return(int(math.sqrt(4**k)))

def get_coord(kmer):
    """given a kmer gets its box position in the grid, returns [x,y]"""
    grid_size = get_grid_size(len(kmer))
    maxx = grid_size
    maxy = grid_size
    posx = 1
    posy = 1
    for char in kmer:
        if char == "C":
            posx += (maxx / 2)
        elif char == "T":
            posy += (maxy / 2)
        elif char == "G":
            posx += (maxx / 2)
            posy += (maxy / 2)
        maxx = (maxx / 2)
        maxy /= 2
    return([int(posx)-1, int(posy)-1])
    
def get_k(jellyfish):
    """asks the jellyfish file what value was used for k"""
    cmd = ["jellyfish", "info", jellyfish]
    result = subprocess.run(cmd, capture_output=True)
    r = re.match(r".*count\s-m\s(\d+)", result.stdout.decode("utf-8"))
    return(int(r.group(1)))
    
def cgr_matrix(jellyfish):
    """runs the cgr process on a jellyfish file and returns a scipy.sparse.dok_matrix object of the CGR"""
    k = get_k(jellyfish)
    grid_size = get_grid_size(k)
    cgr_mat = dok_matrix( (grid_size, grid_size), dtype=numpy.float32)
    for kmer, count in get_kmer_list(jellyfish):
        x, y = get_coord(kmer)
        cgr_mat[x,y] = count
    return(cgr_mat)

def fcgr_matrix(cgr_matrix):
    """ converts a cgr_matrix to a frequency cgr_matrix"""
    return(cgr_matrix / cgr_matrix.sum() )

def log_matrix(cgr_matrix):
    if not is_numpy_array(cgr_matrix):    
        cgr = numpy.log(cgr_matrix.toarray())
        return(dok_matrix(cgr, dtype=numpy.float32))
    else:
        return(numpy.log(cgr_matrix))

def is_numpy_array(cgr_matrix):
    return(type(cgr_matrix) == numpy.ndarray)

def join_cgr(a,b):
    """takes two cgr of shape (n,n) and returns one array of shape (n,n,2)"""
    stacked = numpy.stack((a.toarray(), b.toarray() ), axis = -1)
    return(stacked)

def show_image(cgr_matrix, cmap = "Greys", log = True):
    """displays a cgr matrix in the viewer - is likely scaled by the image generator. works only on 1channel cgr"""
    if log:
        cgr_matrix = log_matrix(cgr_matrix)
    if not is_numpy_array(cgr_matrix):
        cgr_matrix = cgr_matrix.toarray()
    plt.imshow(cgr_matrix, interpolation='nearest', cmap = cmap)
    plt.show()

        
def pad_two_channel_to_three(cgr2chan):
    ''' takes a two channel cgr and pads a third channel containing zero values'''
    dim = numpy.prod(cgr2chan.shape)
    return numpy.insert(cgr2chan, tuple(range(2, dim + 2, 2)), 0).reshape(cgr2chan.shape[0],cgr2chan.shape[1],3)
    
def save_as_csv(cgr_matrix, file = "cgr_matrix.csv", delimiter = ",", fmt='%.3f'):
    """writes 1 channel cgr matrix to file"""
    numpy.savetxt(file, cgr_matrix.toarray(), delimiter = delimiter, fmt = fmt)
    


# TODO
# Improve rendering. Scale array to 0..1. Use HSV, H for single array (others are constant and definable. S = 0.5 default, V = 0.5 default)
# For two array use H and S or H and V, set other as above
# Allow output scaling
# save functions in cgr.view. Merging made internal
# docs