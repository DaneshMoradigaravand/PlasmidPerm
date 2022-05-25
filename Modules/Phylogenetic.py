import os
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

class Phylogenetic:
    def __init__(self, PATH):
        self.PATH=PATH
    
    def binary_sequence_generator(self, input_kmer_pattern, label):
        string_inp="".join([ 'A' if x==0 else 'C' for x in input_kmer_pattern])
        return([">"+label,string_inp])

    def multifasta_fille_generator(self, converted_sequences_phyolgenetic):
        file_output = open(os.path.join(self.PATH,"binary_presence_absence_kmers.fasta"), "w")
        file_output.writelines('\n'.join(converted_sequences_phyolgenetic) + '\n' )
        file_output.close()
        
    def distance_matrix_generator(self):    
        align = AlignIO.read(os.path.join(self.PATH,"binary_presence_absence_kmers.fasta"), "fasta")
        calculator = DistanceCalculator('identity')
        distMatrix = calculator.get_distance(align)
        return(distMatrix)
 
    def distance_tree_file_generator(self,distance_matrix):   
        constructor = DistanceTreeConstructor()
        UPGMATree = constructor.upgma(distance_matrix)
        Phylo.write(UPGMATree, os.path.join(self.PATH,"binary_presence_absence_kmers.tre") , "newick")