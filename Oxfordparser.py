#!/usr/bin/env python

import numpy as np
import collections
import gzip

SNPINFO_FIELDS = ['chrm', 'rsid', 'bp_location', 'ref_allele', 'alt_allele', 'freq']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()

SNP_COMPLEMENT = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

class DosageParser:

    _read_dosage_once = False
    _read_sample_once = False
    _nsample = 0


    def __init__(self, dosage_filename, sample_filename, start, end, chrm):
        self._dosage_filename = dosage_filename
        self._sample_filename = sample_filename
        self._startsnp = start
        self._endsnp = end
        self._chrom = chrm

    @property
    def nsample(self):
        self._read_samples()
        return self._nsample

    @property
    def samplenames(self):
        self._read_samples()
        return self._samplenames


    @property
    def nsnps(self):
        self._read_dosage()
        return self._nsnps


    @property
    def snpinfo(self):
        self._read_dosage()
        return self._snpinfo


    @property
    def dosage(self):
        self._read_dosage()
        return self._dosage



    def _read_dosage(self):
        if self._read_dosage_once:
            return
        self._read_dosage_once = True
        _dosagelist = list()
        _snpinfo = list()
        with gzip.open(self._dosage_filename, 'r') as mfile:
            linenum = 0
            for line in mfile:
                linenum += 1
                if linenum >= self._startsnp and linenum <= self._endsnp:
                    linestrip = line.decode().strip().split()
                    if len(linestrip[3])==1 and len(linestrip[4]) ==1:
                        if linestrip[4] != SNP_COMPLEMENT[linestrip[3]]:
                            this_dosage = np.array(linestrip[5:], dtype=float)
                            genotype = np.array(([this_dosage[i+1] + 2 * this_dosage[i + 2] for i in range(0, len(this_dosage),3)])).astype(float)
                            frequency = (np.sum(genotype) / (len(genotype) * 2))
                            maf = frequency
                            if maf >= 0.20 and maf <= 0.80: 
                                this_snp = SnpInfo(chrm        = int(self._chrom),
                                                   bp_location = int(linestrip[2]),
                                                   rsid        = linestrip[1],
                                                   ref_allele  = linestrip[3],
                                                   alt_allele  = linestrip[4],
                                                   freq        = maf)
                                #this_dosage = np.array(mline[5:], dtype=float)

                                _snpinfo.append(this_snp)
                                _dosagelist.append(genotype)

        self._dosage = np.array(_dosagelist)
        self._nsnps = self._dosage.shape[0]
       # self._nsample = self._dosage.shape[1]
        self._snpinfo = _snpinfo


    def _read_samples(self):
        if self._read_sample_once:
            return
        self._read_sample_once = True
        samplelist = list()
        with open(self._sample_filename, 'r') as samfile:
            header = samfile.readline().strip().split()
            header_types = samfile.readline().strip().split()
            sample = 0
            samplenames = list()
            for line in samfile:
                sample += 1
                samplenames.append(line.strip().split()[0])

        self._nsample = sample
        self._samplenames = samplenames

