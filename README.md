# RAGB Program, 09/04/2017
**Version 1.0.0** 
By Arnon Benshahar
---
## Introduction
We  formalize  a  new  problem  variant  in  gene-block  discovery, denoted **Reference-Anchored  Gene  Blocks(RAGB)**. Given a query sequence **Q** of length **n**, representing the gene-array of a DNA element,a  window  size  bound **d** on  the  length  of  a  substring  of  interest  in **Q**, and a set of target gene sequences **T=T1...Tc**. Our objective is to identify gene-blocks in **T** that are centered in a subset **q** of co-localized genes from **Q**, and contain genomes from **T** in which the corresponding orthologs of the genes from **q** are also co-localized. **RAGB** program is available open-source at https://github.com/benarnon/RAGB and at https://www.cs.bgu.ac.il/~negevcb/RAGB/ where you can find also executable file, supplementary materials, including omitted proofs, figures and data. 

---

## Getting Started
RAGB can be executed only on Linux OS.
In order to run RAGB you need:

- **ncbi-blast-2.6.0+** program (please go to https://www.ncbi.nlm.nih.gov/books/NBK52640/ and follow the **downloading, configuration and execution steps**)

The easiest way to run the project is to execute the executable file ***RAGB*** (https://www.cs.bgu.ac.il/~negevcb/RAGB/RAGBI/). 

```
./RAGB [-q FOLDER] [-g FOLDER] [-o FOLDER] [-d INT] [-n INT] [-iv STRING] [-min_genomes INT] [-min_genes INT] [-rank INT] [-e FLOAT]
```

### Optional Arguments:
- **-q** : folder containing the gbk (GeneBank format) or *IslandViewer* format files (see **Input Format** section) of the centroid genome query (or queries in case of multiple runs).
- **-g** : folder containing all target/reference gbk files.
- **-o** : output folder.
- **-d** : size of the sliding window, the default is 12.
- **-n** : number of processors that you want this script to run on. The default is every CPU that the system has.
- **-iv** : IslandViewer queries format, T for IslandViewer format and F for normal gbk files, the default is F.
- **-min_genomes** : minimum number of genomes in a gene-block, the default is 4.
- **-min_genes** : minimum number of genes in a gene interval, the default is 3.
- **-rank** : minimum ranking score. the default is 20.
- **-e** : eval for the BLAST search, the default is 0.01.


### Input Format
Our program requires two input folders, the references folder, and the queries folder:
1. **Reference folder:** this folder has subfolder for each species. Each subfolder has the gbk file for that species.
Example: inside ```/db``` there is ```/db/species1``` folder which has ```species1.gbk```, this file includes all the genes in *species1*. 
2. **Query folder:** for this folder we have two options:
    1. **IslandViewer Format**:  this is for centroid queries that were predicted by the [IslandViewer 3 tool](http://www.pathogenomics.sfu.ca/islandviewer/browse/) as [genomic islands](https://en.wikipedia.org/wiki/Genomic_island). In this case, the folder has a subfolder for each species that was analysed by the IslandViewer tool. Each folder has the gbk file of the species and a CSV file which includes the details of all the islands that were predicted by the *IslandViewer* tool. 
**Example**:  inside ```/query/species_1``` there are two files: ```speice_1.gbk``` and ```species_1.csv```
        - ```speice_1.gbk``` includes all the genes in *species_1*. 
        - ```species_1.csv```contains the information of all the islands that were predicted in *species_1*.
    2.  **Genebank format (gbk)** : the folder contains the gbk files of all the centroid genome queries.
    

### Data
We supply data so users can test **RAGB**. The data is available at https://www.cs.bgu.ac.il/~negevcb/RAGB/RAGBI/ragbi-data/.
For the reference folder we give the ```/genomes``` directory which contain 33 GeneBank files of different speciess.
For the query folder we give a few options:
- ```/IslandViewer_sample_data```, this folder contain 3 subfolder for 3 different speciess that were analysed by the IsalndViewer tool (this input is in IsalndViewer format).
- ```/IslandViewer_full_data```, this folder contain 29 subfolder for 29 different speciess that were analysed by the IsalndViewer tool (this input is in IsalndViewer format).
- ```/gbk_sample_data```, this folder contain 4 centroid genomes in gbk format (this input is in gbk format).
- ```/gbk_full_data```, this folder contain 29 centroid genomes in gbk format (this input is in gbk format).


## Output
Our program output a directory with the following files:
1. ```targets_file.csv```: a CSV file which contains information about the target/refernece genomes that were given as an input.
2. ```centroid_genome_file.csv```: a CSV file which contains a breif information about the centroid genomes that were given as an input.
3. ```general_results.csv```: a CSV file which contains a breif summary of the results. For each centroid genomes who outputed at least one valid **gene block** it displays the following information: it's top **gene block's** ranking score, number of cliques, number of gene blocks and an avarage number of gene blocks per clique.
4. Result's directory for each centroid genome, it contains two files:
    1.```***centroid_genome_name***_info_file.csv```: this file holds information about all the genes in the centroid genome.
    2.```***centroid_genome_name***_results_file.csv```: this file displays an extensive report for the program results of this centroid genome. The results **gene blocks** are divided into **cliques**. ***Example***:

```
Clique Number	1								
 									
Block Number	1	
Ranking Score	175.6785691449	
Number of Genes from the Query genome	4		
Number of Intervals from the Target genomes	30
 									
Group A (genes)							
Gene Number|Start|End|Gene Id|Gene Name|Attribute|Strand|			
1|001|200|g1|ID_001|att1|-1|	
2|201|500|g2|ID_002|att2|1|
3|505|800|g3|ID_003|att3|-1|
 									
Group B (Genes Intervals)								

Interval 1
Species              genome_1	
Strain	            NC_00001	
Number of Genes	    3

Gene Number|Start|End|Gene Id|Gene Name|Attribute|Strand|Target Gene Id|Target Gene Attribute|Blast E-value|
1|500|700|g1|ID_001|att1|-1|g57|att80|0.001|
2|701|1000|g2|ID_002|att2|1|g58|att73|1e-10|
3|1005|1300|g3|ID_003|att3|1|g59|att22|7e-100|

Interval 2
Species              genome_2
Strain	            NC_00002	
Number of Genes	    4

Gene Number|Start|End|Gene Id|Gene Name|Attribute|Strand|Target Gene Id|Target Gene Attribute|Blast E-value|
3|2200|2495|g3|ID_003|att3|1|g120|att80|0.001|
2|2496|2795|g2|ID_002|att2|1|g121|att73|1e-10|
1|2796|2996|g1|ID_001|att1|1|g122|att22|7e-100|
1|2997|3117|g1|ID_001|att1|1|g123|att22|7e-100|
```
---

## Contact

- Voice: 00972-524541543
- Email: arnon.benshahar@gmail.com

