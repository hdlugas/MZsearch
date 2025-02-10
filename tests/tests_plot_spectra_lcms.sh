#!/bin/bash
cd ../src

echo $'\n\n\n\n\ntest #1'
python plot_spectra_lcms.py \
  --query_data ../data/lcms_query_library.csv \
  --reference_data ../data/lcms_reference_library.csv \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../src/test1.pdf

echo $'\n\n\n\n\ntest #2'
python plot_spectra_lcms.py \
  --similarity_measure tsallis \
  --wf_mz 2 \
  --wf_intensity 0.5 \
  --normalization_method standard \
  --entropy_dimension 0.5 \
  --save_plots ../src/test2.pdf

echo $'\n\n\n\n\ntest #3'
python plot_spectra_lcms.py \
  --similarity_measure tsallis \
  --normalization_method standard \
  --entropy_dimension 2 \
  --save_plots ../src/test3.pdf

echo $'\n\n\n\n\ntest #4'
python plot_spectra_lcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 1.1 \
  --save_plots ../src/test4.pdf


echo $'\n\n\n\n\ntest #5'
python plot_spectra_lcms.py \
  --similarity_measure renyi \
  --normalization_method standard \
  --entropy_dimension 0.9 \
  --save_plots ../src/test5.pdf

echo $'\n\n\n\n\ntest #6'
python plot_spectra_lcms.py \
  --similarity_measure shannon \
  --normalization_method standard \
  --save_plots ../src/test6.pdf

echo $'\n\n\n\n\ntest #7'
python plot_spectra_lcms.py \
  --similarity_measure cosine \
  --save_plots ../src/test7.pdf

echo $'\n\n\n\n\ntest #8'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --wf_mz 0.5 \
  --save_plots ../src/test8.pdf

echo $'\n\n\n\n\ntest #9'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order LFWNCM \
  --similarity_measure cosine \
  --mz_min 200 \
  --mz_max 250 \
  --int_min 50 \
  --int_max 500 \
  --save_plots ../src/test9.pdf

echo $'\n\n\n\n\ntest #10'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --mz_max 100 \
  --save_plots ../src/test10.pdf

echo $'\n\n\n\n\ntest #11'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --int_max 300 \
  --save_plots ../src/test11.pdf

echo $'\n\n\n\n\ntest #12'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --int_min 100 \
  --save_plots ../src/test12.pdf

echo $'\n\n\n\n\ntest #13'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --window_size_centroiding 0.1 \
  --save_plots ../src/test13.pdf

echo $'\n\n\n\n\ntest #14'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --window_size_matching 0.1 \
  --save_plots ../src/test14.pdf

echo $'\n\n\n\n\ntest #15'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --noise_threshold 0.4 \
  --save_plots ../src/test15.pdf

echo $'\n\n\n\n\ntest #16'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWCNLM \
  --similarity_measure cosine \
  --LET_threshold 3 \
  --save_plots ../src/test16.pdf

echo $'\n\n\n\n\ntest #17'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure cosine \
  --LET_threshold 2 \
  --save_plots ../src/test17.pdf

echo $'\n\n\n\n\ntest #18'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure cosine \
  --window_size_centroiding 0.05 \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../src/test18.pdf

echo $'\n\n\n\n\ntest #19'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order WCNLMF \
  --similarity_measure shannon \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../src/test19.pdf

echo $'\n\n\n\n\ntest #20'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order WCNLFM \
  --similarity_measure renyi \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../src/test20.pdf

echo $'\n\n\n\n\ntest #21'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order LWCMNF \
  --similarity_measure tsallis \
  --wf_mz 0.5 \
  --wf_int 1.3 \
  --LET_threshold 2 \
  --save_plots ../src/test21.pdf

echo $'\n\n\n\n\ntest #22'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order CMWNLF \
  --similarity_measure tsallis \
  --LET_threshold 2 \
  --save_plots ../src/test22.pdf

echo $'\n\n\n\n\ntest #23'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order ML \
  --similarity_measure tsallis \
  --LET_threshold 3 \
  --save_plots ../src/test23.pdf

echo $'\n\n\n\n\ntest #24'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order MNL \
  --similarity_measure tsallis \
  --noise_threshold 0.1 \
  --LET_threshold 3 \
  --save_plots ../src/test24.pdf

echo $'\n\n\n\n\ntest #25'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order MW \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --save_plots ../src/test25.pdf

echo $'\n\n\n\n\ntest #26'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order WM \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --save_plots ../src/test26.pdf

echo $'\n\n\n\n\ntest #27'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order FWL \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 4 \
  --save_plots ../src/test27.pdf

echo $'\n\n\n\n\ntest #28'
python plot_spectra_lcms.py \
  --spectrum_preprocessing_order MCWL \
  --similarity_measure tsallis \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 4 \
  --save_plots ../src/test27.pdf


echo $'\n\n\n\n\ntest #29'
python plot_spectra_lcms.py \
  --similarity_measure tsallis \
  --high_quality_reference_library True \
  --wf_mz 0.55 \
  --wf_int 1.1 \
  --LET_threshold 3 \
  --save_plots ../src/test29.pdf

echo 'Finished Testing'