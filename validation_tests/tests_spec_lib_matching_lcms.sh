#!/bin/bash

source activate base
conda activate ezsearch_env

cd /home/hunter/mass_spec/entropies/JOSS/scripts

:<< 'END'


END

echo $'\n\n\n\n\ntest #1'
python spec_lib_matching_lcms.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --n_top_matches_to_save 4 \
  --print_id_results True \
  --output_identification /home/hunter/mass_spec/entropies/JOSS/output_lcms_identification.csv \
  --output_similarity_scores /home/hunter/mass_spec/entropies/JOSS/output_lcms_all_similarity_scores.csv

echo $'\n\n\n\n\ntest #2'
python spec_lib_matching_lcms.py \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --n_top_matches_to_save 2 \
  --print_id_results True \

echo $'\n\n\n\n\ntest #3'
python spec_lib_matching_lcms.py \
  --similarity_measure tsallis \
  --normalization_method standard \
  --entropy_dimension 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #4'
python spec_lib_matching_lcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --n_top_matches_to_save 2 \
  --print_id_results True


echo $'\n\n\n\n\ntest #5'
python spec_lib_matching_lcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #6'
python spec_lib_matching_lcms.py \
  --similarity_measure shannon \
  --normalization_method standard \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #7'
python spec_lib_matching_lcms.py \
  --similarity_measure cosine \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #8'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --wf_mz 0.5 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #9'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #10'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --mz_max 300 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #11'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --int_max 400 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #12'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --int_min 80 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #13'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --window_size_centroiding 0.1 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #14'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --window_size_matching 0.1 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #15'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --noise_threshold 0.1 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #16'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #17'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order WCMNLF \
  --similarity_measure cosine \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #18'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order WCNMLF \
  --similarity_measure cosine \
  --window_size_centroiding 0.05 \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #19'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order WCNLMF \
  --similarity_measure shannon \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #20'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order WCNLFM \
  --similarity_measure renyi \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #21'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order WCLNMF \
  --similarity_measure tsallis \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 3 \
  --n_top_matches_to_save 2 \
  --print_id_results True

echo $'\n\n\n\n\ntest #22'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure tsallis \
  --LET_threshold 3 \
  --n_top_matches_to_save 1 \
  --print_id_results True

echo $'\n\n\n\n\ntest #23'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order ML \
  --similarity_measure tsallis \
  --LET_threshold 3 \
  --n_top_matches_to_save 1 \
  --print_id_results True

echo $'\n\n\n\n\ntest #24'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order MNL \
  --similarity_measure tsallis \
  --noise_threshold 0.1 \
  --LET_threshold 3 \
  --n_top_matches_to_save 1 \
  --print_id_results True

echo $'\n\n\n\n\ntest #25'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order MW \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --n_top_matches_to_save 1 \
  --print_id_results True

echo $'\n\n\n\n\ntest #26'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order FWL \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --n_top_matches_to_save 1 \
  --print_id_results True

echo $'\n\n\n\n\ntest #26'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order MNC \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --n_top_matches_to_save 1 \
  --print_id_results True

echo $'\n\n\n\n\ntest #27'
python spec_lib_matching_lcms.py \
  --spectrum_preprocessing_order WM \
  --similarity_measure tsallis \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --n_top_matches_to_save 1 \
  --print_id_results True


echo 'Finished Testing'



