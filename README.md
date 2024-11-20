# localplateau 
This tool is an adaption from [plateau](https://plateau.bcp.fu-berlin.de/).

## General purpose
The tool can be used to identify and quantify shared consensus epitopes. 

### How to use
1. In your command line go to the directory where this github repository is cloned to.
2. Specify your input as described [here](#Input-files)
3. Enter the following command:
    ```
    python3 plateau_main.py 
    ```




#### First time set up 
1. In the command line go to the directory where you want to clone this github repository to. 

2. Clone this repository with 
    ```bash 
    git clone git@github.com:AG-Walz/localplateau.git
    ```    

### Input files
#### params.yaml
In the params.yaml file the input parameters can be specified. 
| Parameter | Description |
| --- | --- |
| mhcquant_out | Defines the path to the evidence file. |
| fasta_proteome | Defines the path of the proteome. |
| max_step_size | Defines the maximal step size. If the start positions of two peptides differ by that number the peptide are only grouped together if the overlap by a minimum of min_overlap amino acids. |
| min_overlap | Defines the minimal overlap between two epitopes to still be grouped together if the start position between both epitopes differ more than max step size. |
| min_epi_length | Defines the minimum epitope length. This is the minimal length a core epitope has to have. If for a group the whole sequence is shorter than the minimum epitope length the core will be defined as the whole sequence.| 
| seq_column | Defines the column header in the input evidence file that contains the peptide sequences. |
| protacc_column | Defines the column header in the input evidence file that contains the protein accessions of proteins that contain that peptide. |
| out_dir | Defines the directory in which the results will be saved. |
| mod_delimiter | Defines how modifications in the peptides are separated from the sequence. Put a comma separated string here, where the first element specifies the start of a modification and the second element defines the end of a modification. |
| delimiter | Defines the delimiter that separates multiple values in one cell in the input evidence file. |
| prot_accession | Defines proteins, for which the core epitopes and landscape should be visualized. |
| plateau_csv | Defines the path to previous computed plateau results. If this parameter is not None the result will not be computed again. Only the defined proteins will get visualized. |

 

#### evidence file
The evidence file is the output file of a search engine. The following file types are supported: csv, tsv, xlsx
#### proteome file
The proteome file should contain the proteome used for the identification of the peptide sequences. The file should follow the fasta format. 


### Output files
There are two csv output files.
#### plateau_result.csv
The csv contains one protein per row. The different columns contain the following information: 
| column | description |
| --- | --- |
| accession | protein accession |
| sequence | sequence of all peptides aligned to that protein |
| start | start positions of all peptides in the protein | 
| end | end positions of all peptides in that protein | 
| grouped petides start | the start positions of all peptides grouped together for each core |
| grouped peptides end | the end position of all peptides grouped together for each core | 
| grouped peptides sequence | peptide sequences that contribute to the same core grouped together |
| sequence group mapping | | 
| landscape | landscape for each group | 
| whole epitopes | whole peptide sequence of each group | 
| core epitopes | core epitope sequence of each group | 
| core epitopes start | start positions of the cores in the protein |
| core epitopes end | end positions of the cores in the protein |

#### evidence_link_group.csv
This csv contains all the information from the initial evidence file. In addition there are two more columns. The column **whole_epitopes** lists the whole sequence of each group associated with the peptide sequence of that row. The column **core_epitopes** lists the core sequences of each group associated with the peptide sequence of that row. 