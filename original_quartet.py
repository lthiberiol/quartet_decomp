from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from Bio import Phylo as phy
import pandas as pd
import numpy as np
import os
from os.path import isfile, isdir
from os import listdir, mkdir, system, unlink
from sys import argv, exit
from itertools import combinations
import subprocess
import ete2
from scipy.spatial.distance import squareform

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
    def __init__( self, raxml_path='raxml', mafft_path='mafft', target_folder='.', bootstrap=100, num_threads=1, fasttree='fasttree' ):
        self.raxml    = raxml_path
        self.mafft    = mafft_path
        self.boot     = bootstrap
        self.folder   = target_folder
        self.threads  = num_threads
        self.fasttree = fasttree
        
        self.msa_folder          = '%s/MSAs' %self.folder
        self.tree_folder         = '%s/trees' %self.folder
        self.bs_folder           = '%s/bootstraps' %self.folder
        self.quartet_ml_folder   = '%s/quartet_ml' %self.folder

        try:
            if not isdir(self.msa_folder):
                mkdir( self.msa_folder  )
            if not isdir(self.tree_folder):
                mkdir( self.tree_folder )
            if not isdir(self.bs_folder):
                mkdir( self.bs_folder   )
        except:
            raise
            exit()

    def run_msa(self, fasta):
        infile  = '%s/%s'          %( self.folder, fasta)
        outfile = '%s/%s.aln' %( self.msa_folder, fasta)

        if system( '%s --auto --reorder %s > %s' %(self.mafft, infile, outfile) ):
            exit( '**Error while running:\n\t%s --auto --reorder %s > %s.aln ' %(self.mafft, fasta, self.msa_folder, fasta) )
        else:
            return '%s.aln' %fasta

    def reference_ml_tree(self, fasta):
        infile  = '%s/%s'          %( self.msa_folder, fasta)
        outfile = '%s/%s.fastTree' %( self.tree_folder, fasta)


        if system( '%s  -wag -gamma -out %s %s' %(self.fasttree, outfile, infile) ):
            exit( '**Error while running:\n\tfastTree' )

        tree = ete2.Tree( outfile )
        tree.resolve_polytomy()
        tree.write( outfile='%s-no_polytomies' %outfile )

        return ( '%s.fastTree' %fasta, '%s.fastTree-no_polytomies' %fasta)
    
    def create_bootstrap_replicates(self, (fasta, tree) ):
        outfolder = '%s/%s' %( self.bs_folder, fasta  )
        inaln     = '%s/%s' %( self.msa_folder, fasta )
        intree    = '%s/%s' %( self.tree_folder, tree )

        if not isdir( outfolder ):
            try:
                mkdir( outfolder )
            except:
                raise
                exit()

        with cd( outfolder ):
            if not isfile( fasta ):
                system( 'ln -s ../../%s' %inaln)
            if not isfile( 'reference.tre' ):
                system( 'ln -s ../../%s ./reference.tre' %intree)

            if system( '%s -f j -b 12345 -N %i -s %s -n %s -T 2 -m PROTCATWAG -t reference.tre' %(self.raxml, self.boot, fasta, fasta) ):
                raise
                exit( '**Error while running:\n\t%s -f j -b 12345 -N %i -s %s.aln -n %s -T 2 -m PROTCATILG' %(self.raxml, self.boot, fasta, fasta) )

            if isfile( fasta ):
                unlink( fasta )
            if isfile( '%s.reduced' %fasta ):
                unlink( '%s.reduced' %fasta )

        return outfolder

    def pairwise_dists(self, folder):

        with cd( folder ):
            for replicate in listdir('.'):
                if replicate.startswith('RAxML') or replicate == 'reference.tre':
                    continue

                if system( '%s -f x -p 12345 -m PROTGAMMAWAG -s %s -n %s_pw -T 2 -t reference.tre' %(self.raxml, replicate, replicate) ):
                    exit( '**Error while running:\n\traxml -f x' )

                with open('RAxML_distances.%s_pw' % replicate) as ml_dist:
                    observed_taxa    = set()
                    taxa_order       = []
                    condensed_matrix = []
                    for line in ml_dist.readlines():
                        taxon1, taxon2, dist = line.split()

                        condensed_matrix.append( float(dist) )

                        if not taxon1 in observed_taxa:
                            observed_taxa.add( taxon1 )
                            taxa_order.append( taxon1 )
                    taxa_order.append( taxon2 )

                    dist_matrix = pd.DataFrame( index= taxa_order, columns=taxa_order, data=squareform(condensed_matrix) )
                    dist_matrix.to_csv( '%s-pw.tab' %replicate, sep='\t')

    def run_ml_mapping( self, (fasta, tree) ):
        outfolder = self.quartet_ml_folder
        inaln     = '%s/%s' %( self.msa_folder, fasta )
        intree    = '%s/%s' %( self.tree_folder, tree )

        if system( self.raxml+' -f q -p 12345 -s '+inaln+' -n '+fasta+' -T 2 -m PROTCATWAG -t '+intree+' -w '+outfolder ):
            raise
            exit( '**Error while running:\n\t'+self.raxml+' -f q -p 12345 -s '+inaln+' -n '+fasta+' -T 2 -m PROTCATWAG -t '+intree+' -w '+outfolder )

        return 'RAxML_quartets.'+fasta

    def ml_quartets(self, infile):
        infile = '%s/%s' %(self.quartet_ml_folder, infile)
        
        topologies_lnl = open(infile).xreadlines()

        #
        #throw away the first two lines
        topologies_lnl.next()
        topologies_lnl.next()

        codes = {}
        for line in topologies_lnl:
            line = line.split()
            if not line:
                break
            codes[line[1]] = line[0]

        #
        #throw away next line
        topologies_lnl.next()

        best_topology_prop = []
        while topologies_lnl:
            try:
                lines = [topologies_lnl.next(), topologies_lnl.next(), topologies_lnl.next() ]
            except:
                break

            lnl = []
            topologies = []
            for line in lines:
                line = line.split(': ')
                lnl.append( float(line[1]) )
#                topologies.append( frozenset( [ frozenset( [ codes[line[0]], codes[line[2]] ] ),
#                                                frozenset( [ codes[line[6]], codes[line[8]] ] ) ]
#                                            ) 
#                                 )

            lnl  = np.asarray(lnl)
            tmp  = np.exp( lnl - lnl.max() )
            prob = tmp - tmp.sum()

            best_topology_prop.append( prob.max() )

        return best_topology_prop


    def bs_quartets(self, dist_matrix ):

        quartets = {}

        for quartet in combinations(dist_matrix.index, 4):
            Dab = dist_matrix.loc[quartet[0], quartet[1]]
            Dac = dist_matrix.loc[quartet[0], quartet[2]]
            Dad = dist_matrix.loc[quartet[0], quartet[3]]
            Dbc = dist_matrix.loc[quartet[1], quartet[2]]
            Dbd = dist_matrix.loc[quartet[1], quartet[3]]
            Dcd = dist_matrix.loc[quartet[2], quartet[3]]
            topology = set()
            if Dab + Dcd < Dac + Dbd:
                topology.add(frozenset([quartet[0], quartet[1]]))
                topology.add(frozenset([quartet[2], quartet[3]]))
            elif Dac + Dbd < Dad + Dbc:
                topology.add(frozenset([quartet[0], quartet[2]]))
                topology.add(frozenset([quartet[1], quartet[3]]))
            else:
                topology.add(frozenset([quartet[0], quartet[3]]))
                topology.add(frozenset([quartet[1], quartet[2]]))
            topology = frozenset(topology)
            quartet = frozenset(quartet)
            if topology not in quartets[quartet]:
                quartets[quartet][topology] = 0
            quartets[quartet][topology] += 1
