#!/bin/bash

source activate base
conda activate ezsearch_env

cd /home/hunter/mass_spec/entropies/JOSS/scripts

Rscript get_lcms_library_from_mgf.R \
  --input_path /home/hunter/mass_spec/entropies/JOSS/data/GNPS-SELLECKCHEM-FDA-PART1.mgf \
  --output_path /home/hunter/mass_spec/entropies/JOSS/data/GNPS-SELLECKCHEM-FDA-PART1.csv




