#!/bin/bash
### specify which shell to use
#$ -S /bin/bash
### memory
#$ -l mem_free=8G
### PATH to where STDOUT files should be saved
#$ -o /data/SKDngs01/pg/collaborators/junke/eigenTrait_DEG/Raw
### PATH to where STDERR files should be saved
#$ -e /data/SKDngs01/pg/collaborators/junke/eigenTrait_DEG/Raw
### Specify the working directory if you need to. 
#$ -wd /data/SKDngs01/pg/collaborators/junke/eigenTrait_DEG/Raw
### load env from .bashrc
source ~/.bashrc
 
#### typical R script submission
R -q --vanilla < /data/SKDngs01/pg/collaborators/junke/eigenTrait_DEG/script/DEG_array.R
exit 0

