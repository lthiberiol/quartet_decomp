from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from Bio import Phylo as phy
import pandas as pd
import numpy as np
import os
from os.path import isfile, isdir
from os import listdir, mkdir, system
from sys import argv, exit
from itertools import combinations
import subprocess

class cd: 
    """
    Context manager for changing the current working directory
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

class original:
    """ 
    Methods required for the original quartet decomposition!
    """
    def __init__(self, raxml_path='raxml', mafft_path='mafft', target_folder='.', bootstrap=100, num_threads=1):
        self.raxml         = raxml_path
        self.mafft         = mafft_path
        self.boot          = bootstrap
        self.target_folder = target_folder
        self.threads       = num_threads

    def run_msa(self, fasta):
        if system( '%s --auto --reorder %s > %s.aln' %(self.mafft, fasta, fasta) ):
            exit( '**Error while running:\n\t%s --auto --reorder %s > %s.aln ' %(self.mafft, fasta, fasta) )

    def create_bootstrap_replicates(self, fasta):
        if not isdir( '%s_bs' %fasta ):
            try:
                mkdir( '%s_bs' %fasta )
            except:
                raise
                exit()

        with cd( '%s_bs' %fasta ):
            if not isfile( fasta ):
                system( 'ln -s ../%s' %fasta)

            if system( '%s -f j -b 12345 -N %i -s %s -n %s -T 2 -m PROTCATILG' %(self.raxml, self.boot, fasta, fasta) ):
                exit( '**Error while running:\n\t%s -f j -b 12345 -N %i -s %s.aln -n %s -T 2 -m PROTCATILG' %(self.raxml, self.boot, fasta, fasta) )

    def reference_ml_tree(self, fasta):
    
    def pairwise_dists(self, folder):
        if system( '%s -f x -p 12345 -m PROTGAMMAILG -s %s -n %s_pw -T 2' %(self.raxml, fasta, fasta) ):
            exit( '**Error while running:\n\traxml -f x' )

        
