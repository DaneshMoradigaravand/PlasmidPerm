# Machine Learning aided Prediction of Plasmid Permsiveness 

## Table of contents
1. [Citation](#citation)
2. [Introduction](#content)
3. [Installation](#installation)
4. [Manual](#manual)
5. [Output](#output)
6. [Supplemental files](#package)
7. [Contact](#contact)

----

### Citation <a name="citation"></a>

This package has been developed by the Moradigaravand and Kreft lab as part of the following paper, currently under review:

*Plasmid Permissiveness of Wastewater Microbiomes can be Predicted from 16S rRNA phylogeny by Machine Learning.*

Citation will be added upon the final publication of the manuscript.

----

### Introduction <a name="content"></a>
The package predicts a relative value for permissiveness from any arbitrary 16s rRNA input data, using a random forest model. Plasmid permissiveness is the ability of recipient bacteria to receive external DNA through the mechanism of conjugation. Prediction is made in the relative model, i.e. the permissivenss value is compared with other strains in the training dataset and reported as the % of strains having a smaller permissiveness in the training dataset. Furthermore, the package reports the closest systematic type, based on the selected rank, to the input sequence. The prediction is made for the broad-range plasmids of **pB10**, **pKJK5** and **RP4**.  

----
### Installation <a name="installation"></a>

There are three ways to run the tool:

- The package may be downloaded and run as a binary file, **plasmidperm.bin**, as ./plasmidperm.bin

- The tool is available on DockerHub and may be fetched and run using the following commmands:

```
docker pull daneshmoradigaravand/plasmidperm:tagname
docker run -v $PWD:/data --rm -it plasmidperm-docker ./app.py -i input_fasta_file -o /data/output_report
```
*input_fasta_file* and *output_report* should be names as the input fasta and output report files, respectively.

- streamlite application. The graphical interface can be run usin streamlit. You need to navigate into the **streamlit** directory and launch the application, usin the following command:

```
streamlit run app.py
```
The commmand provides a link to the following front web application:

<p align="center">
<img src="https://user-images.githubusercontent.com/35295619/170274133-a3cac858-f320-4444-9426-de689adcf249.png" width="500" />
</p>

----
### Manual <a name="installation"></a>

The tools is initiated using the binary command. The help instruction is called using -h option. Note in the multifasta file base **U** need to be replaced by **T**.   

```
usage: plasmidperm.bin [-h] -i INPUT -o OUTPUT [-p {pB10,pKJK5,RP4}] [-r {Kingdom,Phylum,Class,Order,Family,Genus}] [-t TREE]

optional arguments:
  -h, --help            show this help message and exit
  -p {pB10,pKJK5,RP4}, --plasmid {pB10,pKJK5,RP4}
                        Plasmid specific prediction (default: pKJK5)
  -r {Kingdom,Phylum,Class,Order,Family,Genus}, --rank {Kingdom,Phylum,Class,Order,Family,Genus}
                        The closest rank to the input sequence (default: Order)
  -t TREE, --tree TREE  Produce phylogeny tree newick format (default True) (default: True)

Required arguments:
  -i INPUT, --input INPUT
                        input multifasta 16s rRNA (default: None)
  -o OUTPUT, --output OUTPUT
                        Output file name (default: stdout)

 ```
 
 ----
### Output <a name="output"></a>

The tool produces the following files 

1. **output_file.csv** file with the followinig structure:
 Input	plasmid	Closest Family	Greater Than % Baseline Population

| Tag	          | Sequence      | Plasmid | Closest Systematic Rank  | Greater Than % Baseline Population|
| ------------- |:-------------:| :-----: | :-----------------------:| ---------------------------------:|
| seq1          | AGCTGTGGGTTTA |  pB10   |   Pseudomonadaceae       |               95                  |
| seq2          | AACCCGCGAGGAA |  pB10   |    Aeromonadaceae        |               65                  |


***Tag***: The tags corresponding to the tas in the multifasta input file.

***Sequence***: The sequences from the multifasta file.

***Plasmid***: The plasmid for which perissiveness of reccipients are prredicted.

***Closest Systematic Rank***: The closest rank based on the Eucledian distance between the kmer profile for the sequence and sequuences in the training dataset.

***Greater Than % Baseline Population***: The percentage of isolates in the training dataset which had a permissiveness smmaller than the predicted permissiveness.

2. **binary_presence_absence_kmers.fasta** The file contain binary semim-sequence file for the presence and absence of kmers denoted by **A** and **C** bases, respectively.

3. **binary_presence_absence_kmers.tre** The phylogenetic neighbour-joining tree in the newick format made from the sequence distance matrix of the kmer sequences. The tree can be visulaized in Figtree.

----
### Supplemental files <a name="package"></a>

These files are located in the **SupplementalFiles** folder and include:

1. **Training_Code.ipynb** The python code used for training the model. 

2. **input_file.csv** The input files used for training.

----
### Contact <a name="contact"></a>
For queries, please contact [Danesh Moradigaravand](mailto:d.moradigaravand@bham.ac.uk?subject=[GitHub]), Data-Driven Microbiology lab, Center for Computational Biology, University of Birmingham. 
 


-----
