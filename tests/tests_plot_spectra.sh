#!/bin/bash

cd /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1


echo $'\n\n\n\n\ntest #1'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test1.pdf

echo $'\n\n\n\n\ntest #2'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test2.pdf

echo $'\n\n\n\n\ntest #2'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test2.pdf

echo $'\n\n\n\n\ntest #3'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test3.pdf

echo $'\n\n\n\n\ntest #4'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure renyi \
  --chromatography_platform LCMS \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test4.pdf


echo $'\n\n\n\n\ntest #5'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure renyi \
  --chromatography_platform LCMS \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test5.pdf

echo $'\n\n\n\n\ntest #6'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure shannon \
  --chromatography_platform LCMS \
  --normalization_method standard \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test6.pdf

echo $'\n\n\n\n\ntest #7'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test7.pdf

echo $'\n\n\n\n\ntest #8'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --wf_mz 0.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test8.pdf

echo $'\n\n\n\n\ntest #9'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test9.pdf

echo $'\n\n\n\n\ntest #10'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --mz_max 100 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test10.pdf

echo $'\n\n\n\n\ntest #11'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --int_max 300 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test11.pdf

echo $'\n\n\n\n\ntest #12'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --int_min 100 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test12.pdf

echo $'\n\n\n\n\ntest #13'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --window_size_centroiding 0.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test13.pdf

echo $'\n\n\n\n\ntest #14'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --window_size_matching 0.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test14.pdf

echo $'\n\n\n\n\ntest #15'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --noise_threshold 0.4 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test15.pdf

echo $'\n\n\n\n\ntest #16'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --LET_threshold 3 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test16.pdf

echo $'\n\n\n\n\ntest #17'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test17.pdf

echo $'\n\n\n\n\ntest #18'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure cosine \
  --chromatography_platform LCMS \
  --window_size_centroiding 0.05 \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test18.pdf

echo $'\n\n\n\n\ntest #19'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order WCNLMF \
  --similarity_measure shannon \
  --chromatography_platform LCMS \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test19.pdf

echo $'\n\n\n\n\ntest #20'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order WCNLFM \
  --similarity_measure renyi \
  --chromatography_platform LCMS \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test20.pdf

echo $'\n\n\n\n\ntest #21'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order LWCMNF \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test21.pdf

echo $'\n\n\n\n\ntest #22'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test22.pdf

echo $'\n\n\n\n\ntest #23'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order ML \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --LET_threshold 3 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test23.pdf

echo $'\n\n\n\n\ntest #24'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order MNL \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --noise_threshold 0.1 \
  --LET_threshold 3 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test24.pdf

echo $'\n\n\n\n\ntest #25'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order MW \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test25.pdf

echo $'\n\n\n\n\ntest #26'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order WM \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test26.pdf

echo $'\n\n\n\n\ntest #27'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order FWL \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 4 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test27.pdf

echo $'\n\n\n\n\ntest #28'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --spectrum_preprocessing_order MCWL \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 4 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test27.pdf

echo $'\n\n\n\n\ntest #29'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test29.pdf

echo $'\n\n\n\n\ntest #30'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --query_spectrum_ID 212 \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test30.pdf

echo $'\n\n\n\n\ntest #31'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/lcms_reference_library.csv \
  --reference_spectrum_ID "Malyngamide J M+H" \
  --similarity_measure tsallis \
  --chromatography_platform LCMS \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test31.pdf

echo $'\n\n\n\n\ntest #32'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure cosine \
  --normalization_method standard \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test32.pdf

echo $'\n\n\n\n\ntest #33'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test33.pdf

echo $'\n\n\n\n\ntest #34'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test34.pdf

echo $'\n\n\n\n\ntest #35'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure tsallis \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test35.pdf

echo $'\n\n\n\n\ntest #36'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test36.pdf

echo $'\n\n\n\n\ntest #37'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test37.pdf

echo $'\n\n\n\n\ntest #38'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure shannon \
  --normalization_method standard \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test38.pdf

echo $'\n\n\n\n\ntest #39'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure cosine \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test39.pdf

echo $'\n\n\n\n\ntest #40'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --wf_mz 0.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test40.pdf

echo $'\n\n\n\n\ntest #41'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test41.pdf

echo $'\n\n\n\n\ntest #42'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order FNLW\
  --similarity_measure cosine \
  --mz_max 100 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test42.pdf

echo $'\n\n\n\n\ntest #43'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --int_max 300 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test43.pdf


echo $'\n\n\n\n\ntest #44'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --similarity_measure cosine \
  --int_min 100 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test44.pdf

echo $'\n\n\n\n\ntest #45'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.1 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test45.pdf

echo $'\n\n\n\n\ntest #46'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order NFLW \
  --similarity_measure cosine \
  --noise_threshold 0.4 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test46.pdf

echo $'\n\n\n\n\ntest #47'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order FNLW \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test47.pdf

echo $'\n\n\n\n\ntest #48'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test48.pdf

echo $'\n\n\n\n\ntest #49'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test49.pdf

echo $'\n\n\n\n\ntest #50'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --query_spectrum_ID ID_2 \
  --chromatography_platform GCMS \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test50.pdf

echo $'\n\n\n\n\ntest #51'
python plot_spectra.py \
  --query_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_query_library.csv \
  --reference_data /home/hunter/mass_spec/entropies/JOSS/data/gcms_reference_library.csv \
  --chromatography_platform GCMS \
  --reference_spectrum_ID 463514 \
  --spectrum_preprocessing_order LWF \
  --similarity_measure cosine \
  --LET_threshold 1.5 \
  --save_plots /home/hunter/mass_spec/entropies/JOSS/scripts/reviewer1/test51.pdf


echo -e '\n\nFinished Testing\n'


