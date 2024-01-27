#-------------------------------------------------------------------------------
# Name:        ChemX
# Purpose:     Virtual compound library generator for drug discovery
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
from rdkit import DataStructs
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import BRICS
from rdkit.Chem import AllChem
import math


class ChemX:
    def __init__(self, target_compound, name):
        self.target = target_compound
        self.calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        self.chembank = open('data/' + name, 'a+')
        self.templates = []

        # Initialize the fragment database
        self.fragment_database()

    def fragment_database(self):
        # Load functional groups from file
        fName = 'data/FunctionalGroups.txt'
        fparams = FragmentCatalog.FragCatParams(1, 6, fName)
        self.fcat = FragmentCatalog.FragCatalog(fparams)

        # Load smiles from the ZINC database
        zinc_file = 'data/smiles_database.csv'
        zinc_suppl = [i.split(',')[1] for i in open(zinc_file, 'r').read().splitlines()][1:]
        zinc_ms = [Chem.MolFromSmiles(i) for i in zinc_suppl]

        # Generate synthetic fragment database
        pre_synthetic_frag_database = [BRICS.BRICSDecompose(i) for i in zinc_ms]
        self.synthetic_frag_database = list(set(chain.from_iterable(pre_synthetic_frag_database)))

    def fragment_writer(self):
        # Write synthetic fragment database to file
        frg_manager = open('data/ZINC_fragments_10k', 'w')
        to_write = '\n'.join(self.synthetic_frag_database)
        frg_manager.write(to_write)
        frg_manager.close()

    def tanimoto_similarity(self, mol1, mol2):
        # Tanimoto similarity calculation
        fp1 = AllChem.GetMorganFingerprint(mol1, 2)
        fp2 = AllChem.GetMorganFingerprint(mol2, 2)
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    def get_similar_fragments(self, target_cpd, n=5, threshold=0.4):
        # Find similar fragments in the synthetic fragment database using Tanimoto similarity
        score_frg = []
        target_mol = Chem.MolFromSmiles(target_cpd)
        for i in self.synthetic_frag_database:
            cur_frag_mol = Chem.MolFromSmiles(i)
            similarity = self.tanimoto_similarity(target_mol, cur_frag_mol)
            if similarity >= threshold:
                score_frg.append((i, similarity))
        sorted_tuples = sorted(score_frg, key=lambda x: x[-1], reverse=True)

        chosen_fragments = [k[0] for k in sorted_tuples[:n]]
        return chosen_fragments

    def fragment_target(self):
        # Decompose the target compound into fragments
        self.target_fragments = list(BRICS.BRICSDecompose(Chem.MolFromSmiles(self.target)))

    def gather_fragments_4alltargetfragments(self):
        # Gather replaceable fragments for all target fragments
        self.all_replaceable_fragments = []
        counter = 0
        for i in self.target_fragments:
            print('Searching the database for replaceable chemical fragment ...' + str(counter + 1))
            target_ds = Chem.MolToSmiles(Chem.MolFromSmiles(i))
            self.all_replaceable_fragments.append(self.get_similar_fragments(target_ds, 5))
            counter += 1
        print('Populating all the possible fragments for different parts of your target drug ...')

    def write_fragments_file_4Target(self, filename):
        # Write replaceable fragments to file
        print('Writing to file all the possible fragments for replacement ...')
        header = ['fragment ' + str(u + 1) for u in range(len(self.target_fragments))]
        top_row = dict(zip(header, self.target_fragments))

        max_len = max(len(fragment_list) for fragment_list in self.all_replaceable_fragments)
        
        # Create a DataFrame with a sufficient number of columns
        df_columns = header + [f'extra_{i+1}' for i in range(max_len - len(header))]
        cv_frame = pd.DataFrame(columns=df_columns)

        # Populate the DataFrame with available data
        for idx, fragment_list in enumerate(self.all_replaceable_fragments):
            fragment_dict = dict(zip(header, fragment_list))
            fragment_dict.update({f'extra_{i+1}': '' for i in range(len(fragment_list), max_len)})
            cv_frame = cv_frame.append(fragment_dict, ignore_index=True)

        # Combine top row and DataFrame
        cv_frame = pd.concat([pd.DataFrame([top_row]), cv_frame], ignore_index=True)

        # Save DataFrame to CSV
        cv_frame.to_csv('data/' + filename + '_fragments', index=False)


    def generate_frag_templates(self):
        # Generate potential compound templates
        self.potential_cpd_templates = list(itertools.product(*self.all_replaceable_fragments))

    def collect_mini_frags_from_each_template(self, current_template):
        # Collect mini fragments from each template
        mini_frags = []
        for each in current_template:
            num_joints = each.count('*')
            numbers = re.findall(r'\d+', each)
            possible_joints = ['[' + str(m) + '*]' for m in numbers]

            for i in range(num_joints):
                for j in possible_joints:
                    if j in each:
                        little_ones = each.replace(j, '').replace('()', ',').rstrip().split(',')
                        if little_ones not in mini_frags:
                            mini_frags.append(little_ones)
        mini_frags_flattened = list(itertools.chain(*mini_frags))
        return mini_frags_flattened

    def combine_frag(self, max_compounds=10):
        # Combine fragments to generate compounds, storing only the first 50 compounds for each input
        self.generate_frag_templates()
        print('Merging fragments together to generate compounds...')

        for current_template in self.potential_cpd_templates:
            fragms = [Chem.MolFromSmiles(x) for x in sorted(current_template)]
            ms = BRICS.BRICSBuild(fragms)
            for i, prod in enumerate(ms):
                if i >= max_compounds:
                    break

                sampler = Chem.MolToSmiles(prod, True)
                print(i, sampler)
                if sampler not in self.templates:
                    self.templates.append(sampler)
                    self.chembank.write(sampler + '\n')


def main():
    # Run ChemX for a sample chemical
    mefloquine = 'OC(C1CCCCN1)C1=CC(=NC2=C(C=CC=C12)C(F)(F)F)'
    print('Chemical accepted into the program.')

    # Get the name for the chemical library
    name = input('Name of your chemical library: ')

    # Initialize and run ChemX
    sample = ChemX(mefloquine, name)
    sample.fragment_target()
    sample.gather_fragments_4alltargetfragments()
    sample.write_fragments_file_4Target(name)
    sample.combine_frag()


if __name__ == '__main__':
    main()
