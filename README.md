# localplateau 
This tool is an adaption from [plateau](https://plateau.bcp.fu-berlin.de/).

## General purpose
The tool can be used to identify and quantify shared consensus epitopes. 

### How to use
1. In your command line go to the directory where this github repository is cloned to.
2. Specify your input as described [here](#Input-files)
3. To compute the consensus epitopes enter the following command:
    ```
    python3 plateau_main.py <PROTEOME_FILE> params.yaml generate-plateau-csv <EVIDENCE_FILE>
    ```
    Replace ```EVIDENCE_FILE``` with the path to your evidence file and ```PROTEOME_FILE``` with the path to the proteome fasta file that was used to generate the evidence file. You can find more detailed information about the input data [here](#input-files).  

   To visualize the landscape of one protein you can use the following command:
    ```
    python3 plateau_main.py <PROTEOME_FILE> params.yaml prot-landscape <PLATEAU_RESULT>
    ```

#### First time set up 
##### Create a conda environment
If miniconda is not installed on your system yet, install miniconda by running the following commands:
```bash
  mkdir -p ~/miniconda3
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
  bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
  rm ~/miniconda3/miniconda.sh
  source .bashrc
```

After a successful installation of miniconda the conda environment can be created and activated by executing the following commands. 
```bash
  conda create --name localplateau python pandas=2.2.3 pyyaml=6.0.2 matplotlib=3.10.0 biopython=1.84
  conda activate localplateau
```

##### Clone the localplateau repository 
```bash 
  git clone git@github.com:AG-Walz/localplateau.git
```

### Input files
#### params.yaml
In the params.yaml file the input parameters can be specified. The file should look like this:
```
parameters:
  seq_column: 
  protacc_column: 
  start_column: 
  end_column: 
  out_dir: 
  mod_pattern: 
  delimiter: 
  prot_accession: 
  plateau_csv: 
  min_overlap: 
  max_step_size: 
  min_epi_length: 
```
The description of each parameter can be found in the table below.
| Parameter | Description |
| --- | --- |
| max_step_size | Defines the maximal step size between two peptides to still be grouped to the same core. If the start positions of two peptides differ by that number the peptides are only grouped together if they overlap by a minimum of min_overlap amino acids.|
| min_overlap | Defines the minimal overlap between two epitopes to still be grouped together if the start position between the epitopes differs more than max step size. |
| min_epi_length | Defines the minimum epitope length. This is the minimal length a core epitope has to have. If for a group the whole sequence is shorter than the minimum epitope length the core will be defined as the whole sequence.| 
| seq_column | Defines the column header in the input evidence file that contains the peptide sequences. |
| protacc_column | Defines the column header in the input evidence file that contains the protein accessions of proteins that contain that peptide. |
| out_dir | Defines the directory in which the results will be saved. |
| mod_pattern | Defines how modifications in the peptides are separated from the sequence. Put a comma separated string here, where the first element specifies the start of a modification and the second element defines the end of a modification. |
| delimiter | Defines the delimiter that separates multiple values in one cell in the input evidence file. |
| [prot_accession] | Defines proteins, for which the core epitopes and landscape should be visualized. Multiple parameters should be separated by commas. |
| [plateau_csv] | Defines the path to previous computed plateau results. If this parameter is not None the result will not be computed again. Only the defined proteins will get visualized. |

Parameters in brackets are optional. An example params.yaml file can be found [here](params.yaml).

#### evidence file
The evidence file is the output file of a search engine. The following file types are supported: csv, tsv, xlsx.

#### proteome file
The proteome file should contain the proteome used for the identification of the peptide sequences. The file should follow the fasta format. 


### Output files
There are two csv output files and optionally pdfs visualizing the landscape and core epitopes of a protein.

#### plateau_result.csv
The csv contains one protein per row. The different columns contain the following information: 
| column | description |
| --- | --- |
| accession | protein accession |
| sequence | sequence of all peptides aligned to that protein |
| start | start positions of the peptides in the protein | 
| end | end positions of the peptides in the protein | 
| grouped peptides start | the start positions of all peptides grouped together for each core |
| grouped peptides end | the end position of all peptides grouped together for each core | 
| grouped peptides sequence | peptide sequences that contribute to the same core grouped together |
| sequence group mapping | the epitope core group of each peptide | 
| landscape | landscape for each group | 
| whole epitopes | whole peptide sequence of each group | 
| core epitopes | core epitope sequence of each group | 
| core epitopes start | start positions of the cores in the protein |
| core epitopes end | end positions of the cores in the protein |

#### evidence_link_group.csv
The evidence_link_group.csv contains all the information from the initial evidence file. In addition there are two more columns. The column **whole_epitopes** lists the whole sequence of each group associated with the peptide sequence of that row. The column **core_epitopes** lists the core sequences of each group associated with the peptide sequence of that row. 

#### landscape visualization
![An example landscape of the protein sp|P62736|ACTA_HUMAN](landscape_example.png)
The script visualizes the landscape and core epitopes of each protein specified in the prot_accession input variable. 

The more intense regions of the histogram represent the core epitopes of the protein. The different colors correspond to distinct core epitopes. Lighter areas correspond to the peptides associated with the respective core epitope.


## Workflow


