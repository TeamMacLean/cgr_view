#!/usr/bin/env python

import os
import subprocess
import math
import numpy
from scipy.sparse import dok_matrix
import re
import matplotlib.pyplot as plt
import skimage.color
import skimage.io
import skimage.transform
import tempfile




def estimate_genome_size(fasta):
    """guesses genome size from fasta file size, assumes 1 byte ~= 1 base"""
    return (os.path.getsize(fasta))


def run_jellyfish(fasta, k, out):
    """runs jellyfish on fasta using k kmer size, produces db file out"""
    genome_size = estimate_genome_size(fasta)
    cmd = ["jellyfish", "count", "-m", str(k), "-s", str(genome_size), fasta, "-o", out]
    result = subprocess.run(cmd)
    return (result)


def get_kmer_list(jellyfish):
    """runs jellyfish dump on jellyfish. captures output as a generator stream. each item = [kmer, count]"""
    cmd = ["jellyfish", "dump", "-c", jellyfish]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in proc.stdout:
        yield line.decode("utf-8").rstrip().split(" ")
    proc.wait()
    proc.stdout.close()


def get_grid_size(k):
    return (int(math.sqrt(4 ** k)))


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
    return ([int(posx) - 1, int(posy) - 1])


def get_k(jellyfish):
    """asks the jellyfish file what value was used for k"""
    cmd = ["jellyfish", "info", jellyfish]
    result = subprocess.run(cmd, capture_output=True)
    r = re.match(r".*count\s-m\s(\d+)", result.stdout.decode("utf-8"))
    return (int(r.group(1)))


def cgr_matrix(jellyfish):
    """runs the cgr process on a jellyfish file and returns a scipy.sparse.dok_matrix object of the CGR"""
    k = get_k(jellyfish)
    grid_size = get_grid_size(k)
    cgr_mat = dok_matrix((grid_size, grid_size), dtype=numpy.float32)
    for kmer, count in get_kmer_list(jellyfish):
        x, y = get_coord(kmer)
        cgr_mat[x, y] = count
    return (cgr_mat)



def join_cgr(cgrs):
    """takes tuple of cgr of shape (n,n) and returns one stacked array of size (n,n, len(cgrs) )"""
    cgrs = (cgr.toarray for cgr in cgrs)
    return numpy.dstack(cgrs)

def save_as_csv(cgr_matrix, file="cgr_matrix.csv", delimiter=",", fmt='%.3f'):
    """writes 1 channel cgr matrix to file"""
    numpy.savetxt(file, cgr_matrix.toarray(), delimiter=delimiter, fmt=fmt)


def make_blanks_like(a, h=1.0,s=1.0,v=1.0):
    return numpy.full_like(a, h), numpy.full_like(a, s), numpy.full_like(a, v)


def scale_cgr(cgr_matrix):
    return (cgr_matrix / max(cgr_matrix.values())).toarray()


def blocky_scale(im, nR, nC):
    nR0 = len(im)     # source number of rows
    nC0 = len(im[0])  # source number of columns
    return [[ im[int(nR0 * r / nR)][int(nC0 * c / nC)]
             for c in range(nC)] for r in range(nR)]


def resize_rgb_out(rgb,resize):
    r = blocky_scale(rgb[:, :, 0], resize, resize)
    g = blocky_scale(rgb[:, :, 1], resize, resize)
    b = blocky_scale(rgb[:, :, 2], resize, resize)
    return numpy.dstack((r, g, b))


def is_cgr_matrix(obj):
    return type(obj) == scipy.sparse.dok.dok_matrix


def draw_cgr(cgr_matrices, h=0.8, s=0.5, v=1.0, main="s", show = True, write = True, out = "cgr.png", resize = False):
    if is_cgr_matrix(cgr_matrices): #one channel
        draw_single_cgr(cgr_matrices, h=h, s=s, v=v, main=main, show = show, write = write, out = out, resize = resize)
    elif all( [is_cgr_matrix(o) for o in cgr_matrices] ) and len(cgr_matrices) ==2: #all cgr matrices
        draw_two_cgrs(cgr_matrices, v = v, show = show, write = write, out = out, resize = resize)
    elif all( [is_cgr_matrix(o) for o in cgr_matrices] ) and len(cgr_matrices) == 3 :
        draw_three_cgrs(cgr_matrices)
    else:
        raise Exception("don't know what to do, cgr_matrices must be one cgr_matrix or a tuple of 2 or 3 cgr_matrices.")


def draw_single_cgr(cgr_matrix, h=0.8, s=0.5, v=1.0, main="s", show = True, write = True, out = "cgr.png", resize = False):

    scaled = scale_cgr( cgr_matrix )

    h_blank, s_blank, v_blank = make_blanks_like(scaled, h,s,v)

    hsv = None
    if main == "h":
        hsv = numpy.dstack((scaled, s_blank, v_blank))
    elif main == "s":
        hsv = numpy.dstack((h_blank, scaled, v_blank))
    elif main == "v":
        hsv = numpy.dstack((h_blank, s_blank, scaled))

    rgb = skimage.color.hsv2rgb(hsv)

    if show:
        plt.imshow(rgb)
        plt.show()

    if write:
        if resize:
            rgb = resize_rgb_out(rgb, resize )
        skimage.io.imsave(out, rgb)


def draw_two_cgrs(cgr_matrices, v = 1.0, show = True, write = True, out = "cgr.png", resize = False ):
    scaled_l = [scale_cgr(cgrm) for cgrm in cgr_matrices]
    v_blank = make_blanks_like(scaled_l[0], v=v)[2]
    hsv_stack = numpy.dstack((scaled_l[0], scaled_l[1], v_blank))
    rgb = skimage.color.hsv2rgb(hsv_stack)

    if show:
        draw(rgb)

    if write:
        write_out(rgb, out, resize)


def draw_three_cgrs(cgr_matrices,show = True, write = True, out = "cgr.png", resize = False):
    scaled_t = (scale_cgr(cgrm) for cgrm in cgr_matrices)
    hsv_stack = numpy.dstack(scaled_t)
    rgb = skimage.color.hsv2rgb(hsv_stack)

    if show:
        draw(rgb)

    if write:
        write_out(rgb, out, resize)

def draw(rgb):
    plt.imshow(rgb)
    plt.show()

def write_out(rgb, out, resize):
    if resize:
        rgb = resize_rgb_out(rgb, resize)
    skimage.io.imsave(out, rgb)


def stack_cgrs(cgr_matrices):
    return numpy.stack(cgr_matrices, axis=-1)


def save_cgr(cgr_obj, outfile = "cgr"):
    '''can be a one or more dimensional array'''
    numpy.save(outfile, cgr_obj)


def many_seq_record_to_one_cgr(fa_file, k):
    temp_jf = tempfile.NamedTemporaryFile()
    run_jellyfish(fa_file, k, temp_jf.name)
    cgr1 = cgr_matrix(temp_jf)
    temp_jf.close()
    return cgr1



def many_seq_record_to_many_cgr(seq_record, k):
    temp_fa = tempfile.NamedTemporaryFile()
    temp_jf = tempfile.NamedTemporaryFile()
    SeqIO.write(seq_record, temp_fa.name, "fasta")
    run_jellyfish(temp_fa, k, temp_jf.name)
    cgr1 = cgr_matrix(temp_jf)
    temp_fa.close()
    temp_jf.close()
    return cgr1


def from_fasta(fasta_file, outfile="my_cgrs", as_single = False, k = 7):
    if as_single:
        cgr_t = many_seq_record_to_one_cgr(fasta_file, k)
        save_cgr(cgr_t, outfile=outfile)
    else:
        cgr_t = ( many_seq_record_to_many_cgr(seq_record, k) for seq_record in SeqIO.parse(fasta_file, "fasta") )
        save_cgr(cgr_t, outfile = outfile)

# TODO
# test fasta saving code.
