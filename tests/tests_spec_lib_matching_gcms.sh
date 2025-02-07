#!/bin/bash

cd ../src

echo $'\n\n\n\n\ntest #0'
python spec_lib_matching_gcms.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --similarity_measure cosine \
  --normalization_method standard \
  --print_id_results True \

echo $'\n\n\n\n\ntest #1'
python spec_lib_matching_gcms.py \
  --query_data ../data/gcms_query_library.csv \
  --reference_data ../data/gcms_reference_library.csv \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --n_top_matches_to_save 4 \
  --print_id_results True \
  --output_identification ../output_gcms_identification.csv \
  --output_similarity_scores ../output_gcms_all_similarity_scores.csv

echo $'\n\n\n\n\ntest #2'
python spec_lib_matching_gcms.py \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --n_top_matches_to_save 2 \
  --print_id_results True \

echo $'\n\n\n\n\ntest #3'
python spec_lib_matching_gcms.py \
  --similarity_measure tsallis \
  --normalization_method standard \
  --entropy_dimension 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #4'
python spec_lib_matching_gcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #5'
python spec_lib_matching_gcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #6'
python spec_lib_matching_gcms.py \
  --similarity_measure shannon \
  --normalization_method standard \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #7'
python spec_lib_matching_gcms.py \
  --similarity_measure cosine \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #8'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --wf_mz 0.5 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #9'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #10'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_max 100 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #11'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --int_max 300 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #12'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.1 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #13'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.4 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #14'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #15'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --wf_int 1.2 \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #16'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --likely_reference_IDs ../data/likely_gcms_ids.csv \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #17'
python spec_lib_matching_gcms.py \
  --spectrum_preprocessing_order LWF \
  --high_quality_reference_library True \
  --similarity_measure cosine \
  --n_top_matches_to_save 2 \
  --print_id_results True


echo $'\n\nFinished Testing\n'



