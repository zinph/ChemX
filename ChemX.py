#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      AnalogueGenerator
#
# Created:     02/03/2019
# Copyright:   (c) kzphy 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
import numpy as np
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit.Chem import DataStructs

class AnalogueGenerator:
    def __init__(self, target_compound):
        self.directory = 'G:/My Drive/NCSU/DiamondHacks/'
        self.target = target_compound

        self.fragment_collector('macrocycles_IDs.csv')
        self.frag_database = [self.fcat.GetEntryDescription(i) for i in range(self.fcat.GetNumEntries())]

    def fragment_collector(self,name):
        fName= 'C:/RDKit_2017_03_2/Data/FunctionalGroups.txt'
        fparams = FragmentCatalog.FragCatParams(1,6,fName)
        self.fcat=FragmentCatalog.FragCatalog(fparams)


##        macrocycle_file = 'macrocycles_IDs.csv'
##        suppl = [i.split(',')[0] for i in open(self.directory+name,'r').read().splitlines()][1:]       # read all the macrocycle smiles from file
##        ms = [Chem.MolFromSmiles(i) for i in suppl]     # mols of macrocycles

        zinc_file = 'smiles_database.csv'
        zinc_suppl = [i.split(',')[1] for i in open(self.directory + name,'r').read().splitlines()][1:]
        zinc_ms = [Chem.MolFromSmiles(i) for i in zinc_suppl]

        self.fcgen=FragmentCatalog.FragCatGenerator()  # fragment generator

        for m in ms:
            self.fcgen.AddFragsFromMol(m,self.fcat)

##        for m in zinc_ms:
##            self.fcgen.AddFragsFromMol(m,self.fcat)

    def fragment_writer(self):
        frg_manager = open(self.directory + 'ZINC_fragments_half','w')
        to_write = '\n'.join(frag_database)
        frg_manager.write(to_write)

    def find_replaceable_analogues(self):
        fpgen = FragmentCatalog.FragFPGenerator()
        fcat = FragmentCatalog.FragCatalog(fparams)
        m = Chem.MolFromSmiles(self.target)






def main():

    sample = AnalogueGenerator()


if __name__ == '__main__':
    main()
