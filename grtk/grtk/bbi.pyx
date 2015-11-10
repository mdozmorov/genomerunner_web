"""
A read-only Python API for BigWig and BigBED files.
"""
# TODO:
# add a "normalize" kwarg to "mean"
# handle invalid BBI files
# handle exons (in Python or C++?)
# API for getting BBI header information
# Parallelize the "same region across many files" case in BigWigSet
# Make BigBED class
# make BigWigSet read metadata
# Make both BigBED and BigWig iterable

cdef extern from "mm.hpp":
    cdef cppclass MMFile:
        pass

cdef extern from "bbi.hpp":
    cdef struct BED:
        char* chrom
        int start, end
        char* name
        float score
        char strand
        char* rest

    cdef struct Summary:
        int length, covered
        float sum, mean0, mean

    cdef cppclass BigWigFile:
        BigWigFile(char* path)
        Summary summary(BED)

cdef class BigWig:
    cdef BigWigFile *h

    def __cinit__(self, str path):
        bpath = path.encode("UTF-8")
        cdef const char* pPath = bpath
        self.h = new BigWigFile(pPath)

    def __dealloc__(self):
        del self.h

    def mean(self, str chrom, int start, int end):
        cdef BED bed
        chr = chrom.encode("UTF-8")
        bed.chrom = chr
        bed.start = start
        bed.end = end
        cdef Summary s = self.h.summary(bed)
        return s.mean0

import os
import numpy

class BigWigSet(object):
    def __init__(self, dir):
        self._handles = []
        dir = os.path.abspath(dir)
        for p in os.listdir(dir):
            if p.endswith(".bw"):
                self._handles.append(BigWig(os.path.join(dir,p)))
    
    def mean(self, str chrom, int start, int end):
        return numpy.array([h.mean(chrom,start,end) for h in self._handles])
       
def test():
    bws = BigWigSet("data/")
    print bws.mean("chr1", 0, 100000)
