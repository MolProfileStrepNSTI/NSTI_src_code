# Description
Source code to reproduce analyses for the paper *Molecular profiling of tissue biopsies reveals unique signatures associated with streptococcal necrotizing soft tissue infections*

# System Requirements
## Software Dependencies
- Python 3.6
- pandas

- R 
- Boruta

# Installation
"""
# clone git repository to a local directory
git clone https://....git

# install dependencies
conda install requirements.txt
"""

# Instructions
"""
# Create data for Figure 3 (execution time <1min)

    python ./figure3.py

# Create data for Figure 4 (execution time ~2min)

    python ./Figure4.py

# Create data and images for Figure 7A + 7C-E (execution time ~25min)

    python ./Figure7_Classifier_train_test_val.py

# Create image for Figure 7B (execution time <1min)

    Rscript --vanilla ./Figure7B_boruta.R

"""



# License
This work is licensed under the MIT license

# Credits
This application uses Open Source components. You can find the source code of their open source projects along with license information below. We acknowledge and are grateful to these developers for their contributions to open source.

## Project: GOENRICH
- https://github.com/jdrudolph/goenrich
- License (MIT)