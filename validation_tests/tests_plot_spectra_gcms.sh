#!/bin/bash

source activate base
conda activate ezsearch_env

cd /home/hunter/mass_spec/entropies/JOSS/scripts


echo $'\n\n\n\n\ntest #0'
python plot_spectra_gcms.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --similarity_measure cosine \
  --normalization_method standard \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test1.pdf

echo $'\n\n\n\n\ntest #1'
python plot_spectra_gcms.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test1.pdf

echo $'\n\n\n\n\ntest #2'
python plot_spectra_gcms.py \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test2.pdf

echo $'\n\n\n\n\ntest #3'
python plot_spectra_gcms.py \
  --similarity_measure tsallis \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test3.pdf

echo $'\n\n\n\n\ntest #4'
python plot_spectra_gcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test4.pdf

echo $'\n\n\n\n\ntest #5'
python plot_spectra_gcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test5.pdf

echo $'\n\n\n\n\ntest #6'
python plot_spectra_gcms.py \
  --similarity_measure shannon \
  --normalization_method standard \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test6.pdf

echo $'\n\n\n\n\ntest #7'
python plot_spectra_gcms.py \
  --similarity_measure cosine \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test7.pdf

echo $'\n\n\n\n\ntest #8'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --wf_mz 0.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test8.pdf

echo $'\n\n\n\n\ntest #9'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test9.pdf

echo $'\n\n\n\n\ntest #10'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_max 100 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test10.pdf

echo $'\n\n\n\n\ntest #11'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --int_max 300 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test11.pdf


echo $'\n\n\n\n\ntest #12'
python plot_spectra_gcms.py \
  --similarity_measure cosine \
  --int_min 100 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test12.pdf

echo $'\n\n\n\n\ntest #13'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test13.pdf

echo $'\n\n\n\n\ntest #14'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.4 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test14.pdf

echo $'\n\n\n\n\ntest #15'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test15.pdf

echo $'\n\n\n\n\ntest #16'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test16.pdf

echo $'\n\n\n\n\ntest #17'
python plot_spectra_gcms.py \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/test17.pdf


echo $'\n\nFinished Testing'



