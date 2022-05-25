#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 12:06:11 2022

@author: moradigd
"""
from scipy.spatial import distance
import pickle
import pandas as pd
import numpy as np
import os

class Prediction:
    def __init__(self, PATH):
        self.PATH=PATH
        
    def significance(self, input_value,plasmid):
        input_df=pd.read_csv(os.path.join(self.PATH ,'Metadata',"baseline_"+plasmid+".csv"))
        return(int(len(np.where(input_value > input_df.Value)[0])/len(input_df.Value)*100))


    def processor(self,file_input, label="labels"):
        labels=[]
        sequences=[]
        for i in file_input:
            if ">" in i:
                labels.append(i.replace(">", "").replace('\n', ''))
            else:
                sequences.append(i.replace('\n', ''))
        if label=="labels":
            return labels
        elif label=="sequence":
            return sequences
    
    def prediction(self,unique_kmer_frequency, plasmid):
        pickle_file=open(os.path.join(self.PATH,'Metadata',plasmid+'_model.pkl'),'rb')
        model=pickle.load(pickle_file)
        prediction = model.predict(pd.DataFrame(unique_kmer_frequency))[0]
        return(prediction)
    
    def level_index(self, rank_input):
        ranks=['Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus']
        return(np.where(np.array(ranks)==rank_input)[0][0])
    
    def taxa_extractor(self, input_sequence, plasmid, rank_input='Order'):
        input_dataframe=pd.read_csv(os.path.join(self.PATH ,'Metadata',"Profile_sequences_"+plasmid+".csv"))
        input_taxa=pd.read_csv(os.path.join(self.PATH ,'Metadata',"Taxa_info_"+plasmid+".csv"))
        

        distance_holder=[]
        for i in range(input_dataframe.shape[0]):
            distance_holder.append(distance.euclidean(input_sequence, input_dataframe.iloc[i,:].values))
        return(input_taxa.iloc[np.argmin(distance_holder),  self.level_index(rank_input)])
    
        
        
        

        
    