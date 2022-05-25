#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 19:18:27 2022

@author: moradigd
This the packge for predicting plaasmid permissivenss from 16s rRNA data
"""

import pandas as pd
import warnings
import click
import time
import os
from Modules.Phylogenetic import Phylogenetic
from Modules.Prediction import Prediction

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

warnings.filterwarnings("ignore")
CURRENT_PATH =os.path.dirname(os.path.realpath(__file__))

def main():
    
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    requiredNamed = parser.add_argument_group('Required arguments')
    requiredNamed.add_argument('-i', '--input', help="input multifasta 16s rRNA", required=True)
    requiredNamed.add_argument('-o', '--output', help='Output file name', default='stdout',required=True)
    
    parser.add_argument("-p", "--plasmid", default="pKJK5",choices=['pB10', 'pKJK5', 'RP4'], help="Plasmid specific prediction")
    parser.add_argument("-r", "--rank", default='Order',choices=['Kingdom', 'Phylum', 'Class','Order', 'Family','Genus'], help="The closest rank to the input sequence")
    parser.add_argument("-t", "--tree", help="Produce phylogeny tree newick format (default True)")
    parser.set_defaults(tree=True)
    args = vars(parser.parse_args())
        
    plasmid=args["plasmid"]
    rank = args["rank"]
    tree=args["tree"]
    
    
    Prediction_object =Prediction(CURRENT_PATH)
    Phylogenetic_object = []
    if(tree):
        Phylogenetic_object = Phylogenetic(CURRENT_PATH)
    
    print("Reading input fasta file\n")
    input_file = open(os.path.join(CURRENT_PATH, args["input"]), "r")
    files_lines=input_file.readlines()


    unique_sequences=pd.read_csv(os.path.join(CURRENT_PATH,'Metadata',"unique_sequences_"+plasmid+".csv"))
    unique_sequences=unique_sequences.iloc[:,1].values
        
    input_sequence_list=Prediction_object.processor(files_lines,label="sequence")
    label_list=Prediction_object.processor(files_lines, label="labels")
    
    prediction_output=[]
    significance_output=[]
    taxa_output=[]
    converted_sequences_phyolgenetic_tree=[]
    
    print("Predicting permissiveness for "+ plasmid)
    for j in range(len(input_sequence_list)):
        print("       Predicting for "+ label_list[j])
        with click.progressbar(range(100), fill_char='=', empty_char=' ') as bar:
            for user in bar:
                time.sleep(0.01)

        
        unique_kmer_frequency=[]
        unique_kmer_frequency.append([input_sequence_list[j].count(i) for i in unique_sequences])
        
        prediction_output.append(Prediction_object.prediction(unique_kmer_frequency,plasmid))
        significance_output.append(Prediction_object.significance(Prediction_object.prediction(unique_kmer_frequency, plasmid), plasmid))
        
        if(tree):
            converted_sequences_phyolgenetic_tree=converted_sequences_phyolgenetic_tree+Phylogenetic_object.binary_sequence_generator(unique_kmer_frequency[0], label_list[j])
        taxa_output.append(Prediction_object.taxa_extractor(unique_kmer_frequency,plasmid,rank))
    
    prediction_output = pd.concat([ pd.DataFrame(input_sequence_list), pd.DataFrame(prediction_output), pd.DataFrame(taxa_output) , pd.DataFrame(significance_output)  ], axis=1)
    prediction_output.columns = ['Sequence Input','plasmid','Closest '+ rank,'Greater Than % Baseline Population']
    prediction_output.plasmid=plasmid
    prediction_output.index=label_list
    prediction_output=prediction_output.rename_axis('Tag')
    prediction_output.to_csv(os.path.join(CURRENT_PATH,args["output"]+".csv"))
    
    #phyogeny generator
    if(tree):
        print("\nProducing phylogenetic tree")
        Phylogenetic_object.multifasta_fille_generator(converted_sequences_phyolgenetic_tree)    
        distance_matrix=Phylogenetic_object.distance_matrix_generator()
        Phylogenetic_object.distance_tree_file_generator(distance_matrix)

    

if __name__ == "__main__":
    main()
