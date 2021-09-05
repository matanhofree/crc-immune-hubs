# CRC Immune Hub expression program analysis 

Principal analysis code used in ["Spatially organized multicellular immune hubs in human colorectal cancer"](https://www.cell.com/cell/fulltext/S0092-8674(21)00945-4).

<br></br>
## Installation 
Clone the repository from github:
```
git clone https://github.com/matanhofree/crc-immune-hubs.git
```
### Prerequisites:
* Matlab (>R2019a)
* python (3.7) & jupyter 

<br></br>
## Data download
Data from used in the analysis can be downloaded from the supplemental website to a data folder.

For example you can use wget:
```
wget -r -np -nH --cut-dirs=2 -R "index.html*" https://portals.broadinstitute.org/crc-immune-hubs/extra/data/colon10x_default/
wget -r -np -nH --cut-dirs=2 -R "index.html*" -R "*.tar.gz" https://portals.broadinstitute.org/crc-immune-hubs/extra/data/matlab/
wget -r -np -nH --cut-dirs=2 -R "index.html*" https://portals.broadinstitute.org/crc-immune-hubs/extra/data/cNMF_tSNE/
```
<br></br>
## Jupyter notebooks
The notebooks folder of this repo contains jupyter notebooks for reproducing the key figures included in our publication. <br>
These notebooks use the Calysto/matlab_kernel to remotely run matlab code in jupyter notebooks instance. See the following for instructions on how to configure this jupyter_kernel: [[1]](https://am111.readthedocs.io/en/latest/jmatlab_install.html), [[2]](https://github.com/Calysto/matlab_kernel).

<br></br>
## Under construction 
(Last update: Sept. 5th)<br>
**Note**: the repositiory will be updated in the coming days to inlclude additional code and examples.<br>
Please check back again soon!


<img src="https://upload.wikimedia.org/wikipedia/commons/1/1a/C%C3%B4ne_orange_-_under_construction.png" width="250" >

