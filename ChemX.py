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

import re
import itertools
import numpy as np
import pandas as pd
from rdkit import Chem
from itertools import chain
from rdkit.Chem import FragmentCatalog
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import BRICS
import math


class ChemX:

    def __init__(self, target_compound,name):
        self.directory = 'G:/My Drive/NCSU/DiamondHacks/ChemX/'
        self.target = target_compound

        self.fragment_database()
##        self.frag_database = [self.fcat.GetEntryDescription(i) for i in range(self.fcat.GetNumEntries())]
        self.calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        self.chembank = open('data/'+name,'a+')
        self.templates = []


    def fragment_database(self):
        fName= 'C:/RDKit_2017_03_2/Data/FunctionalGroups.txt'
        fparams = FragmentCatalog.FragCatParams(1,6,fName)
        self.fcat=FragmentCatalog.FragCatalog(fparams)

##        macrocycle_file = 'macrocycles_IDs.csv'
##        suppl = [i.split(',')[0] for i in open(self.directory+name,'r').read().splitlines()][1:]       # read all the macrocycle smiles from file
##        ms = [Chem.MolFromSmiles(i) for i in suppl]     # mols of macrocycles

        zinc_file = 'data/smiles_database.csv'
        zinc_suppl = [i.split(',')[1] for i in open(self.directory + zinc_file,'r').read().splitlines()][1:]
        zinc_ms = [Chem.MolFromSmiles(i) for i in zinc_suppl]

        pre_synthetic_frag_database = [BRICS.BRICSDecompose(i) for i in zinc_ms]
        self.synthetic_frag_database = list(set(chain.from_iterable(pre_synthetic_frag_database)))

##        self.fcgen=FragmentCatalog.FragCatGenerator()  # fragment generator

##        for m in ms:
##            self.fcgen.AddFragsFromMol(m,self.fcat)

##        for m in zinc_ms:
##            self.fcgen.AddFragsFromMol(m,self.fcat)

    def fragment_writer(self):
        frg_manager = open(self.directory + 'data/ZINC_fragments_10k','w')
        to_write = '\n'.join(self.synthetic_frag_database)
        frg_manager.write(to_write)
        frg_manager.close()


    def eudis(self,v1, v2):
        dist = [(a - b)**2 for a, b in zip(v1, v2)]
        dist = math.sqrt(sum(dist))
        return dist

    def get_similar_fragments(self, target_cpd,n=20):
        
        score_frg= []
        for i in self.synthetic_frag_database:
            cur_frag_ds = [0 if math.isnan(x) else x for x in self.calc.CalcDescriptors(Chem.MolFromSmiles(i))]
            dst = self.eudis(cur_frag_ds,target_cpd)
            score_frg.append((i,dst))
        sorted_tuples = sorted(score_frg, key=lambda x: x[-1])
        #lowest_scores = heapq.nsmallest(n, [k[1] for k in score_frg])

        chosen_fragments = [k[0] for k in sorted_tuples[:n]]
        #chosen_fragments = [list(score_dict.keys())[list(score_dict.values()).index(j)] for j in lowest_scores]
        return chosen_fragments

    def fragment_target(self):
        self.target_fragments =  list(BRICS.BRICSDecompose(Chem.MolFromSmiles(self.target)))


    def gather_fragments_4alltargetfragments(self):
        
        self.all_replaceable_fragments = []
        counter = 0
        for i in self.target_fragments:
            print('Searching the database for replaceable chemical fragment ...' + str(counter+1))
            target_ds =  [0 if math.isnan(x) else x for x in self.calc.CalcDescriptors(Chem.MolFromSmiles(i))]
            self.all_replaceable_fragments.append(self.get_similar_fragments(target_ds,50))
            counter+=1
        print('Populating all the possible fragments for different parts of your target drug ...')

    def write_fragments_file_4Target(self,filename):
        print('Writing to file all the possible fragments for replacement ...')
        header = ['fragment '+str(u+1) for u in range(len(self.target_fragments))]
        top_row = dict(zip(header, self.target_fragments))
        
        cv_frame = pd.DataFrame(np.column_stack(self.all_replaceable_fragments),columns=header)
        data = []
        data.insert(0,top_row)
        pd.concat([pd.DataFrame(data), cv_frame], ignore_index=True)
        cv_frame.to_csv('data/'+filename+'_fragments',index=False)

    def generate_frag_templates(self):
        self.potential_cpd_templates = list(itertools.product(*self.all_replaceable_fragments))


    def collect_mini_frags_from_each_template(self,current_template):
        '''
        Tackle [##*] issues and split the pieces. The output is used to find the
        best sampler by the most shared number of mini_frags.
        '''
        mini_frags = []
        for each in current_template:
            numjoints = each.count('*')
            numbers = re.findall(r'\d+', each)
            possible_joints = ['['+str(m) +'*]' for m in numbers]
            #print(possible_joints)

            for i in range(numjoints):
                for j in possible_joints:
                    if j in each:
                        littleones = each.replace(j,'').replace('()',',').rstrip().split(',')
                        #print(littleones)
                        if littleones not in mini_frags:
                            mini_frags.append(littleones)
        mini_frags_flattened = list(itertools.chain(*mini_frags))
        return mini_frags_flattened


    def combine_frag(self):
        self.generate_frag_templates()
        print('Merging fragments together to generate compounds...')
        for current_template in self.potential_cpd_templates:
            fragms = [Chem.MolFromSmiles(x) for x in sorted(current_template)]
            ms = BRICS.BRICSBuild(fragms)
            prods = [next(ms) for x in range(1)]
#            mini_frags = self.collect_mini_frags_from_each_template(current_template)
#            percent = len(mini_frags)
#            counter = 0
            for i in range(1):
#                for j in range(len(mini_frags)):
                sampler = Chem.MolToSmiles(prods[i],True)
#                    if mini_frags[j] in sampler:
#                        counter+=1
#                        if counter == percent:
                              
                if sampler not in self.templates:
                    print(sampler)  
                    self.templates.append(sampler)
                    self.chembank.write(sampler+'\n')



def main():
    #directory = 'G:/My Drive/NCSU/DiamondHacks/'

    mefloquine  = 'OC(C1CCCCN1)C1=CC(=NC2=C(C=CC=C12)C(F)(F)F)C(F)(F)F'
    
    print('Chemical accepted into the program.')
    
    
    name = input('Name of your chemical library: ')
    sample = ChemX(mefloquine,name)
    sample.fragment_target()
    sample.gather_fragments_4alltargetfragments()
    sample.write_fragments_file_4Target(name)
    sample.combine_frag()


if __name__ == '__main__':
    main()
