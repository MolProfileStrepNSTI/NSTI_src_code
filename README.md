# Description
Source code to reproduce analyses for the paper ***Molecular profiling of tissue biopsies reveals unique signatures associated with streptococcal necrotizing soft tissue infections***

# System Requirements
## Software Dependencies

Code was tested with the following package versions (see requirements.txt)

### Python 
    
    - Python>=3.6
    - biopython==1.72
    - matplotlib==3.0.2
    - networkx==2.2
    - numpy==1.15.4
    - pandas==0.23.4
    - scikit-learn==0.20.1
    - scipy==1.1.0
    - statsmodels==0.9.0
    - urllib3==1.24.1
    - XlsxWriter==1.1.2

### R 
    - R>=3.5
    - Boruta==6.0.0

# Installation

### clone git repository to a local directory
```
    git clone https://github.com/MolProfileStrepNSTI/NSTI_src_code.git
```

### install dependencies (execution time <5min)

**pip**
```
    python -m pip install -r requirements.txt
```

**conda**
```
    conda install --file requirements.txt
```

# Instructions
## Create source data for Figure 3 (execution time <1min)

- script to creates source data table and a GO ontology subgraph in gml file format in the results/figure3 folder
```
    python ./figure3.py
```

### Create source data for Figure 4 (execution time ~2min)

- script to creates source data table for Figure 4B in the results/figure4 folder
```
    python ./Figure4.py
```

### Create source data and subfigures for Figure 7

- scripts to create source data tables and figures in the results/figure7 folder
- hyperparameter tuning for each classifier is located in the notebooks folder
    as a jupyter notebook or in .html format or 


#### Create data and subfigures for Figure 7A + 7C-E (execution time ~25min)
```
    python ./Figure7_Classifier_train_test_val.py
```

***Results***
- **Figure 7A**
<img src="https://github.com/MolProfileStrepNSTI/NSTI_src_code/blob/master/pics/Fig7A_Model_comparison_ROC.svg?sanitize=true">

- **Figure 7C**
<img src="https://github.com/MolProfileStrepNSTI/NSTI_src_code/blob/master/pics/Fig7C_Feat_Comparison_ROC_RF.svg?sanitize=true">

- **Figure 7D**
<img src="https://github.com/MolProfileStrepNSTI/NSTI_src_code/blob/master/pics/Fig7D_FinalModel_ROC.svg?sanitize=true">

#### Create image for Figure 7B (execution time <1min)
```
    Rscript --vanilla ./Figure7B_boruta.R
```

***Results***
- **Figure 7B**
<img src="https://github.com/MolProfileStrepNSTI/NSTI_src_code/blob/master/pics/Fig7B_Boruta.png">

# License
This work is licensed under the MIT license

# Credits
This application uses Open Source components. You can find the source code of their open source projects along with license information below. We acknowledge and are grateful to these developers for their contributions to open source.

## Project: GOENRICH
- https://github.com/jdrudolph/goenrich
- **License** (MIT)
