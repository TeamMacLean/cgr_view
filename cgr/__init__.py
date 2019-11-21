"""
cgr

A module for creating, saving and drawing k-mer matrices and Chaos Game Representations (CGRs) of nucleotide sequences

Prerequisites
-------------

- Jellyfish

An external program for counting k-mers. Must be accessible on the path. You can install from conda as follows:

    .. highlight:: bash
    .. code-block:: bash

        conda install -c bioconda jellyfish


Quickstart
----------

+ Input fasta file, get cgr

    * one cgr for each entry in the fasta file

    .. highlight:: python
    .. code-block:: python

        cgr.from_fasta("my_seqs.fa", outfile = "my_cgrs", k = 7)

    * just one cgr with all entries in the fasta file (eg for genomes and contigs)

    .. highlight:: python
    .. code-block:: python

        cgr.from_fasta("my_genome.fa", outfile = "genome_cgr", k = 7, as_single = True)



Workflow:
---------

1. make kmer count db in Jellyfish from fasta -> generate cgr from db.
2. optionally merge cgrs into single cgr as separate channels
3. stack all composed cgrs into an array of cgrs
4. save as numpy binary (.npy) files


Usage:
------

1. Import module

    .. highlight:: python
    .. code-block:: python

        import cgr


2. Make kmer count db

    .. highlight:: python
    .. code-block:: python

        cgr.run_jellyfish("test_data/NC_012920.fasta", 11, "11mer.jf")
        cgr.run_jellyfish("test_data/NC_012920.fasta", 10, "10_mer.jf")


2. Load CGRs from kmer count db

    .. highlight:: python
    .. code-block:: python

        cgr1 = cgr.cgr_matrix("/Users/macleand/Desktop/athal-5-mers.jf")
        cgr2 = cgr.cgr_matrix("test_data/five_mer.jf")

3. Draw a cgr and save to file

    * just one cgr, can choose colour (value of 'h') and which channel to put cgr in

    .. highlight:: python
    .. code-block:: python

        cgr.draw_cgr(cgr1, h = 0.64, v = 1.0, out = "my_cgr.png", resize = 1000, main = "s" )

    * two cgrs, first in tuple goes in 'h', second goes in 's'. Can set 'v'

    .. highlight:: python
    .. code-block:: python

        cgr.draw_cgr( (cgr1, cgr1), v = 1.0, out = "two_cgrs.png")

    * three cgrs 'h','s' and 'v' are assigned as order in tuple

    .. highlight:: python
    .. code-block:: python

        cgr.draw_cgr( (cgr1, cgr1, cgr1) )

4. Save a single cgr into a text file

    .. highlight:: python
    .. code-block:: python

        cgr.save_as_csv(cgr1, file = "out.csv")

5. Join n cgrs into one, extending the number of channels ...

    .. highlight:: python
    .. code-block:: python

        merged_cgr = cgr.join_cgr( (cgr1, cgr2, ... ) )

6. Write to file (numpy binary)

    .. highlight:: python
    .. code-block:: python

        cgr.save_cgr("my_cgr, merged_cgr )

7. Input fasta file, get cgr
    * one cgr for each entry in the fasta file

    .. highlight:: python
    .. code-block:: python

        cgr.from_fasta("my_seqs.fa", outfile = "my_cgrs", k = 7)

    * just one cgr with all entries in the fasta file (eg for genomes and contigs)

    .. highlight:: python
    .. code-block:: python

        cgr.from_fasta("my_genome.fa", outfile = "genome_cgr", k = 7, as_single = True)

"""

import os
import subprocess
import math
import numpy
import scipy
import re
import matplotlib.pyplot as plt
import skimage.color
import skimage.io
import skimage.transform
import tempfile
from Bio import SeqIO
from typing import Generator, List, Tuple


def estimate_genome_size(fasta: str) -> int:
    """
    Guesses genome size from fasta file size, assumes 1 byte ~= 1 base

    :param: fasta str -- a fasta file
    :return: int -- approximate genome size in nucleotides

    """
    return (os.path.getsize(fasta))


def run_jellyfish(fasta: str, k: int, out: str) -> int:
    """
    runs Jellyfish on fasta file using k kmer size, produces Jellyfish db file as side effect.

    :param: fasta str -- a fasta file
    :param: k int -- size of kmers to use
    :param: out str -- file in which to save kmer db
    :return: int -- return code of Jellyfish subprocess

    """

    genome_size = estimate_genome_size(fasta)
    cmd = ["jellyfish", "count", "-m", str(k), "-s", str(genome_size), fasta, "-o", out]
    result = subprocess.run(cmd)
    return result.returncode


def get_kmer_list(jellyfish: str) -> Generator[List, str, None]:
    """
    runs jellyfish dump on a Jellyfish DB. Captures output as a generator stream.
    Each item returned is a list [kmer: str, count: str]

    :param: jellyfish str -- a Jellyfish DB file
    :return: Generator -- a list of [kmer string, times_kmer_seen]

    """

    cmd = ["jellyfish", "dump", "-c", jellyfish]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    for line in proc.stdout:
        yield line.decode("utf-8").rstrip().split(" ")
    proc.wait()
    proc.stdout.close()


def get_grid_size(k: int) -> int:
    """
    returns the grid size (total number of elements for a
    cgr of k length kmers

    :param: k int -- the value of k to be used
    :return: int -- the total number of elements in the grid

    """
    return int(math.sqrt(4 ** k))


def get_coord(kmer: str) -> List[int]:
    """
    given a kmer gets the coordinates of the box position in the cgr grid,
    returns as list [x,y] of coordinates

    :param: kmer str -- a string of nucleotides
    :return: coords [x,y] -- the x,y positions of the nucleotides in the cgr

    """
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
    return [int(posx) - 1, int(posy) - 1]


def get_k(jellyfish: str) -> int:
    """
    asks the jellyfish file what value was used for k

    :param: jellyfish str -- jellyfish DB file
    :return: int -- length of k used

    """
    cmd = ["jellyfish", "info", jellyfish]
    result = subprocess.run(cmd, capture_output=True)
    r = re.match(r".*count\s-m\s(\d+)", result.stdout.decode("utf-8"))
    return int(r.group(1))


def get_max_count(jellyfish) -> int:
    """
    estimates the count of the most represented kmer in the jellyfish file by using the last bucket of the
    :param jellyfish:
    :return: int estimated count of the most represented kmer
    """
    cmd = ["jellyfish", "histo", jellyfish ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    counts = []
    for line in proc.stdout:
        counts.append( line.decode("utf-8").rstrip().split(" ") )
    proc.wait()
    proc.stdout.close()
    return int(counts[-1][0])

def cgr_matrix(jellyfish: str) -> scipy.sparse.dok_matrix:
    """
    Main function, creates the cgr matrix, a sparse matrix of type scipy.sparse.dok_matrix

    Runs the cgr process on a jellyfish file and returns a scipy.sparse.dok_matrix object of the CGR with dtype int32
    Only observed kmers are represented, absent coordinates mean 0 counts for the kmer at that coordinate.

    :param: jellyfish str -- jellyfish DB file
    :return: scipy.sparse.dok_matrix -- sparse matrix of kmer counts
    """

    k = get_k(jellyfish)
    max_c = get_max_count(jellyfish)
    dtype_to_use = numpy.uint8
    if max_c > 255:
        dtype_to_use = numpy.uint16
    grid_size = get_grid_size(k)
    cgr_mat = scipy.sparse.dok_matrix((grid_size, grid_size), dtype=dtype_to_use)
    for kmer, count in get_kmer_list(jellyfish):
        x, y = get_coord(kmer)
        cgr_mat[x, y] = count
    return cgr_mat


def join_cgr(cgrs: tuple) -> numpy.ndarray:
    """
    Takes tuple of cgrs of shape (n,n) and returns one stacked array of size (n,n, len(cgrs) )

    :param: cgrs tuple -- tuple of cgrs to be joined
    :return: numpy.ndarray

    """
    return numpy.dstack(cgrs)


def save_as_csv(cgr_matrix: scipy.sparse.dok_matrix, file: str = "cgr_matrix.csv", delimiter: str = ",", fmt: str = '%d'):
    """
    Writes simple 1 channel cgr matrix to CSV file.

    See also numpy.savetxt

    :param: cgr_matrix scipy.sparse.dok_matrix -- cgr_matrix to save
    :param: file str -- filename to write to
    :param: delimiter str -- column separator character
    :param: fmt str -- text format string
    :return: None
    """
    numpy.savetxt(file, cgr_matrix.toarray(), delimiter=delimiter, fmt=fmt)


def make_blanks_like(a: scipy.sparse.dok_matrix , h: float=1.0, s: float=1.0, v: float=1.0) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
    """
    returns tuple of numpy.ndarrays with default values of h,s and v of shape of a

    :param: a scipy.sparse.dok_matrix -- a cgr matrix to make blanks like
    :param: h float -- the values with which to fill the first numpy.ndarray
    :param: s float -- the values with which to fill the second numpy.ndarray
    :param: v float -- the values with which to fill the third numpy.ndarray
    :return: Tuple of numpy.ndarray
    """
    return numpy.full_like(a, h), numpy.full_like(a, s), numpy.full_like(a, v)


def scale_cgr(cgr_matrix: scipy.sparse.dok_matrix) -> scipy.sparse.dok_matrix:
    """
    returns scaled version of cgr_matrix in range 0..1

    :param: cgr_matrix scipy.sparse.dok_matrix -- matrix to scale
    :return: scaled scipy.sparse.dok_matrix

    """
    return (cgr_matrix / max(cgr_matrix.values())).toarray()


def blocky_scale(im: numpy.ndarray, nR: int, nC: int) -> numpy.ndarray:
    """
    Upscales an array in preparation for drawing. By default the array is a square with sqrt(k ** 4) wide and high.
    For many values of k this will be too small to view well on a monitor. This function does a scale operartion that
    increases the size of the image by simply increasing the pixels in each square.

    :param: im numpy.ndarray -- the image to be scaled
    :param: nR int -- the number of height pixels to be in the final image
    :param: nC int -- the number of width pixels to be in the final image
    :return: numpy.ndarray -- upscaled image

    """
    nR0 = len(im)     # source number of rows
    nC0 = len(im[0])  # source number of columns
    return numpy.asarray([[im[int(nR0 * r / nR)][int(nC0 * c / nC)]
                            for c in range(nC)] for r in range(nR)])


def resize_rgb_out(rgb: numpy.ndarray,resize: int) -> numpy.ndarray:
    """

    given an rgb image in one pixel per kmer size, increases size so that the resulting image is resize * resize pixels

    :param: rgb numpy.ndarray -- an RGB image array
    :param: resize -- pixel width (and therefore height) of resulting image
    :return: numpy.ndarray -- resized image with shape (resize, resize)

    """
    r = blocky_scale(rgb[:, :, 0], resize, resize)
    g = blocky_scale(rgb[:, :, 1], resize, resize)
    b = blocky_scale(rgb[:, :, 2], resize, resize)
    return numpy.dstack((r, g, b))


def is_cgr_matrix(obj) -> bool:
    """returns true if obj is a scipy.sparse.dok.dok_matrix object """
    return type(obj) == scipy.sparse.dok.dok_matrix


def draw_cgr(cgr_matrices: scipy.sparse.dok_matrix,
             h: float = 0.8,
             s: float = 0.5,
             v: float = 1.0,
             main: str = "s",
             show: bool = True,
             write: bool = True,
             out: str = "cgr.png",
             resize: bool = False) -> None:

    """Draws cgrs to a file. Allows user to set which of up to 3 provided cgr matrices goes in at which of the H, S or V image channels.
    Typically for one channel, select h to specify the image colour and set cgr as s to change that colour according to counts in cgr.
    Set v to 1.0 for maximum brightness.

    :param: cgr_matrices scipy.sparse.dok_matrix or tuple of scipy.sparse.dok_matrix elements, cgrs to be drawn. Tuple provides order for HSV channels of image.
    :param: h float -- (0..1) value for h channel if not used for cgr data
    :param: s float -- (0..1) value for s channel if not used for cgr data
    :param: v float -- (0..1) value for v channel if not used for cgr data
    :param: main str -- the channel to place the cgr matrix in if a single cgr matrix is passed
    :param: show bool -- render CGR picture to screen
    :param: write -- write CGR picture to file
    :param: out str -- filename to write to
    :param: resize bool or int -- if False no image resizing is done, if an int image is rescaled to resize pixels width and height
    :return: None
    """

    if is_cgr_matrix(cgr_matrices):  #one channel
        draw_single_cgr(cgr_matrices, h=h, s=s, v=v, main=main, show = show, write = write, out = out, resize = resize)
    elif all( [is_cgr_matrix(o) for o in cgr_matrices] ) and len(cgr_matrices) ==2: #all cgr matrices
        draw_two_cgrs(cgr_matrices, v = v, show = show, write = write, out = out, resize = resize)
    elif all( [is_cgr_matrix(o) for o in cgr_matrices] ) and len(cgr_matrices) == 3 :
        draw_three_cgrs(cgr_matrices,  show = show, write = write, out = out, resize = resize)
    else:
        raise Exception("don't know what to do, cgr_matrices must be one cgr_matrix or a tuple of 2 or 3 cgr_matrices.")


def draw_single_cgr(cgr_matrix, h=0.8, s=0.5, v=1.0, main="s", show = True, write = True, out = "cgr.png", resize = False):
    """
    draws a single cgr image, selecting channels and resizing as appropriate

    :param: cgr_matrix scipy.sparse.dok_matrix to be drawn.
    :param: h float -- (0..1) value for h channel if not used for cgr data
    :param: s float -- (0..1) value for s channel if not used for cgr data
    :param: v float -- (0..1) value for v channel if not used for cgr data
    :param: main str -- the channel to place the cgr matrix in
    :param: show bool -- render CGR picture to screen
    :param: write -- write CGR picture to file
    :param: out str -- filename to write to
    :param: resize bool or int -- if False no image resizing is done, if an int image is rescaled to resize pixels width and height
    :return: None

    """
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
            rgb = resize_rgb_out(rgb, resize)
        skimage.io.imsave(out, rgb)


def draw_two_cgrs(cgr_matrices, v = 1.0, show = True, write = True, out = "cgr.png", resize = False ):

    """draws two cgr matrices into a single image. first matrix of tuple becomes h channel, second of tuple becomes v channel

    :param: cgr_matrices  tuple of scipy.sparse.dok_matrix elements, cgrs to be drawn.
    :param: v float -- (0..1) value for v channel
    :param: show bool -- render CGR picture to screen
    :param: write -- write CGR picture to file
    :param: out str -- filename to write to
    :param: resize bool or int -- if False no image resizing is done, if an int image is rescaled to resize pixels width and height
    :return: None

    """
    scaled_l = [scale_cgr(cgrm) for cgrm in cgr_matrices]
    v_blank = make_blanks_like(scaled_l[0], v=v)[2]
    hsv_stack = numpy.dstack((scaled_l[0], scaled_l[1], v_blank))
    rgb = skimage.color.hsv2rgb(hsv_stack)

    if show:
        draw(rgb)

    if write:
        write_out(rgb, out, resize)


def draw_three_cgrs(cgr_matrices,show = True, write = True, out = "cgr.png", resize = False):

    """Draws a tuple of 3 cgr matrices as an image

    :param: cgr_matrices tuple of scipy.sparse.dok_matrix elements, cgrs to be drawn. Tuple provides order for HSV channels of image
    :param: show bool -- render CGR picture to screen
    :param: write -- write CGR picture to file
    :param: out str -- filename to write to
    :param: resize bool or int -- if False no image resizing is done, if an int image is rescaled to resize pixels width and height
    :return: None

    """

    scaled_t = (scale_cgr(cgrm) for cgrm in cgr_matrices)
    hsv_stack = numpy.dstack(scaled_t)
    rgb = skimage.color.hsv2rgb(hsv_stack)

    if show:
        draw(rgb)

    if write:
        write_out(rgb, out, resize)

def draw(rgb: numpy.ndarray) -> None:
    """
    renders RGB array on the screen.

    :param: rgb numpy.ndarray -- RGB channel image
    """
    plt.imshow(rgb)
    plt.show()

def write_out(rgb: numpy.ndarray, out: str, resize: int) -> None:
    """
    writes RGB array as image

    :param rgb: numpy.ndarray -- RGB channel image
    :param out: str file to write to
    :param resize: bool or int. If False will not resize, if int will resize image up to that size
    :return: None
    """
    if resize:
        rgb = resize_rgb_out(rgb, resize)
    skimage.io.imsave(out, rgb)


def stack_cgrs(cgr_matrices: Tuple) -> numpy.ndarray:
    """
    stacks cgrs of tuple of N numpy.ndarrays of shape (w,h)
    returns ndarray of ndarrays of shape (w,h,N)

    :param cgr_matrices: tuple of cgr_matrices
    :return: numpy.ndarray

    """
    cgr_t = tuple(c.toarray() for c in cgr_matrices)
    return numpy.stack(cgr_t, axis=-1)


def save_cgr(cgr_obj: numpy.ndarray, outfile: str = "cgr") -> None:
    """
    Saves cgr_obj as numpy .npy file.
    cgr_obj one or more dimensional numpy.ndarray.
    saves as ndarray not dokmatrix, so can be loaded in regular numpy as collections of cgrs

    :param cgr_obj: numpy.ndarray constructed cgr_object to save
    :param outfile: str file
    :return: None

    """
    numpy.save(outfile, cgr_obj, allow_pickle=True)

def load_npy(file: str) -> numpy.ndarray:
    """
    loads numpy .npy file as ndarray.
    Useful for restoring collections of cgrs but resulting array is not compatible directly with
    drawing methods here.

    :param file str -- numpy .npy file to load
    :return: numpy.ndarray

    """
    return numpy.load(file, allow_pickle=True)


def many_seq_record_to_one_cgr(fa_file: str, k: int) -> scipy.sparse.dok_matrix:
    """
    Reads many sequence records in a FASTA file into a single CGR matrix, treating all sequence records as if they are one sequence, EG as if for a genome sequence in Chromosomes.
    :param fa_file: str FASTA FILE name
    :param k: int length of k to use
    :return: scipy.sparse.dok_matrix
    """
    temp_jf = tempfile.NamedTemporaryFile()
    run_jellyfish(fa_file, k, temp_jf.name)
    cgr1 = cgr_matrix(temp_jf.name)
    temp_jf.close()
    return cgr1


def many_seq_record_to_many_cgr(seq_record: SeqIO.FastaIO, k: int) -> scipy.sparse.dok_matrix:
    """
    :param seq_record: Bio.SeqIO FASTA record
    :param k: int size of k to use
    :return: scipy.sparse.dok_matrix
    """
    temp_fa = tempfile.NamedTemporaryFile()
    temp_jf = tempfile.NamedTemporaryFile()
    SeqIO.write(seq_record, temp_fa.name, "fasta")
    run_jellyfish(temp_fa.name, k, temp_jf.name)
    cgr1 = cgr_matrix(temp_jf.name)
    temp_fa.close()
    temp_jf.close()
    return cgr1


def from_fasta(fasta_file: str, outfile: str = "my_cgrs", as_single: bool=False, k: int = 7) -> None:
    """
    Factory function to load in a FASTA file and generate a binary .npy of CGRs

    :param fasta_file: str FASTA file to load
    :param outfile: str outfile to save
    :param as_single: bool If True treats all entries as single sequence and return one CGR. If False, treats all entries individually and returns many CGR
    :param k: int length of kmer to use
    :return: None
    """

    if as_single:
        cgr1 = many_seq_record_to_one_cgr(fasta_file, k)

        save_cgr(cgr1.toarray(), outfile=outfile )
    else:
        cgr_t = tuple( many_seq_record_to_many_cgr(seq_record, k) for seq_record in SeqIO.parse(fasta_file, "fasta") )
        cgr1 = stack_cgrs(cgr_t)
        save_cgr(cgr1, outfile = outfile )

# TODO
# test new dtype switching cgr matrix function - try using from_fasta
