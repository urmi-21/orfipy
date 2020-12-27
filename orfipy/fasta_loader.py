#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 16:44:04 2020

@author: usingh
"""
import pyfastx
orig = "ACTG"
comp = "TGAC"
transtab = str.maketrans(orig, comp)

class FastxWrapper():
    
    def __init__(self,file,ftype='fasta'):
        """
        init object

        Parameters
        ----------
        file : infile path fasta or fastq
            DESCRIPTION.
        ftype : format, optional
            Format of in file. The default is 'fasta'.

        Returns
        -------
        None.

        """
        self.file=file
        self.ftype=ftype
        #read the file and init keys
        self.read_file()

    def read_file(self):
        """
        Read the file into pyfastx object

        Returns
        -------
        None.

        """
        if self.ftype=='fastq':
            self.seqs=pyfastx.Fastq(self.file)
        else:
            self.seqs=pyfastx.Fasta(self.file)
        
        self.keys=list(self.seqs.keys())
    
    def get_fastx_object(self):
        return self.seqs

    def get_seq(self,key,rc=False):
        """
        Return sequence by key

        Parameters
        ----------
        key : str
            sequence name.
        rc : bool, optional
            If True return reverse complement. The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        #return seq
        if not rc: return str(self.seqs[key]).upper()
        
        if self.ftype=='fasta':
            return self.seqs[key][:].antisense.upper()
        else:
            return self.revcomp(str(self.seqs[key])).upper()
        
    def revcomp(self,seq):
        return seq.translate(transtab)[::-1]