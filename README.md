# MZsearch
Command-line Python tool to perform spectral library matching to identify chemical compounds with host of preprocessing transformations and similarity measures (Cosine and three entropy-based similarity measures). MZsearch is capable of performing spectral library matching with respect to either gas chromatography - mass spectrometry (GC-MS) or liquid chromatography - mass spectrometry (LC-MS) data.

## Table of Contents
- [1. Create conda environment and install dependencies](#create-conda-env)
- [2. Functionality](#functionality)
   - [2.1 Spectrum Preprocessing Transformations](#spec-preprocessing-transformations)
   - [2.2 Similarity Measures](#similarity-measures)
- [3. Usage](#usage)
   - [3.1 Obtain LC-MS or GC-MS library from MGF, mzML, or cdf file](#process-data)
   - [3.2 Run spectral library matching](#run-spec-lib-matching)
   - [3.3 Plot a query spectrum against a reference spectrum before and after spectrum preprocessing transformations](#plotting)
- [4. Bugs/Questions?](#bugs-questions)

<a name="create-conda-env"></a>
## 1. Create conda environment and install dependencies
MZsearch requires the Python dependencies Matplotlib, NumPy, Pandas, and SciPy, Pyteomics, and netCDF4. Specifically, this software was validated with python=3.12.4, matplotlib=3.8.4, numpy=1.26.4, pandas=2.2.2, scipy=1.13.1, pyteomics=4.7.2, and nedCDF4=1.6.5, although it may work with other versions of these tools. For instructions on installing conda on your system, see: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Once conda is installed, clone the GitHub repository, navigate to the directory containing the environment.yml file, and create (and activate) a conda environment with the necessary dependencies with the commands:
```
conda env create -f environment.yml
conda activate mzsearch_env
```

To return to your base environment, you can deactivate the mzsearch_env environment with the command:
```
conda deactivate
```

<a name="functionality"></a>
## 2. Functionality

<a name="spec-preprocessing-transformations"></a>
## 2.1 Spectrum Preprocessing Transformations
Functionality implementing the following spectrum preprocessing
transformations is offered in MZsearch:

-   Filtering: Given user-defined parameters (mz_min,mz_max),
    (int_min,int_max) and spectrum $I$ with m/z values
    $(m_{1},m_{2},...,m_{n})$ and intensities $(x_{1},x_{2},...,x_{n})$,
    the transformed spectrum $I^{\star}$ consists of the peaks
    $(m_{i},x_{i})$ in $I$ such that mz_min $\leq m_{i}\leq$ mz_max and
    int_min $\leq x_{i}\leq$ int_max.
    -   Weight Factor Transformation: Given a pair of user-defined weight
    factor parameters $(\text{a,b})$ and spectrum $I$ with m/z values
    $(m_{1},m_{2},...,m_{n})$ and intensities $(x_{1},x_{2},...,x_{n})$,
    the transformed spectrum $I^{\star}$ has the same m/z values as $I$
    and has intensities given by
    $I^{\star}:=(m_{1}^{\text{a}}\cdot x_{1}^{\text{b}},m_{2}^{\text{a}}\cdot x_{2}^{\text{b}},...,m_{n}^{\text{a}}\cdot x_{n}^{\text{b}})$.

-   Low-Entropy Transformation: Given a user-defined low-entropy
    threshold parameter $T$ and spectrum $I$ with intensities
    $(x_{1},x_{2},...,x_{n})$, $\sum_{i=1}^nx_i = 1$, and Shannon
    entropy $H_{Shannon}(I)=-\sum_{i=1}^{n}x_{i}\cdot ln(x_{i})$, the
    transformed spectrum intensities
    $I^{\star}=(x_{1}^{\star},x_{2}^{\star},...,x_{n}^{\star})$ are such
    that, for all $i\in\{1,2,...,n\}$, $x_{i}^{\star}=x_{i}$ if
    $H_{Shannon}(I)\geq T$ and
    $x_{i}^{\star}=x_{i}^{\frac{1+H_{Shannon}(I)}{1+T}}$ if
    $H_{Shannon}(I)\textless T$.

-   Centroiding (only applicable to LC-MS data): Given a user-defined
    window-size parameter $w_{centroiding}$ and a spectrum $I$ with m/z
    values $(m_{1},m_{2},...,m_{n})$ and intensities
    $(x_{1},x_{2},...,x_{n})$, the transformed spectrum $I^{\star}$
    merges adjacent peaks $(m_{i},x_{i}),(m_{i+1},x_{i+1})$ into the
    peak
    $(\frac{m_{i}\cdot x_{i}+m_{i+1}\cdot x_{i+1}}{x_{i}+x_{i+1}},x_{i}+x_{i+1})$
    if $|m_{i}-m_{i+1}|\textless w_{centroiding}$ for
    $i\in\{1,2,...,n-1\}$. This centroiding procedure generalizes to
    more than two peaks whose m/z values are within distance
    $w_{centroiding}$ of each other.

-   Noise Removal: Given a user-defined noise removal parameter $r$ and
    a spectrum $I$ with intensities $(x_{1},x_{2},...,x_{n})$, noise
    removal removes peaks from $I$ with
    $x_{j}\textless r\cdot\text{max}(\{x_{1},x_{2},...,x_{n}\})$ for
    $j\in\{1,2,...,n\}$.

-   Matching (only applicable to LC-MS data): Given a user-defined
    window-size parameter $w_{matching}$ and two spectra $I$, $J$ with
    m/z ratios $(a_{1},a_{2},...,a_{n}), (b_{1},b_{2},...,b_{m})$ and
    intensities $(x_{1},x_{2},...,x_{n}), (y_{1},y_{2},...,y_{m})$,
    respectively, of which we would like to measure the similarity
    between, the matching procedure outputs two spectra
    $I^{\star},J^{\star}$ containing the same number of peaks with
    $I^{\star}$ and $J^{\star}$ having transformed intensities and
    identical m/z ratios. Specifically, for a given peak $(a_{i},x_{i})$
    of $I$, if there are no peaks $(b_{j},y_{j})$ in $J$ with
    $|a_{i}-b_{j}|\textless w_{matching}$, then the peak $(a_{i},x_{i})$
    remains in $I^{\star}$ and the peak $(a_{i},0)$ is included in
    $J^{\star}$. If there is at least one peak $(b_{j},y_{j})$ with
    $|a_{i}-b_{j}|\textless w_{matching}$, then the peak $(a_{i},x_{i})$
    remains in $I^{\star}$ and the peak
    $(a_{i},\sum_{j\text{ such that }|a_{i}-b_{j}|\textless w_{matching}}b_{j})$
    is included in $J^{\star}$. This procedure is applied when
    transposing the roles of $I$ and $J$ as well.

-   Normalization: Prior to computing entropy - regardless of whether in
    the context of performing the low-entropy transformation or
    computing an entropy-based similarity score - the intensities of
    each spectrum must be normalized to sum to 1 in order to represent a
    probability distribution. To normalize a given spectrum $I$ with
    intensities $(x_{1},x_{2},...,x_{n})$ into spectrum $I^{\star}$ with
    intensities $(x_{1}^{\star},x_{2}^{\star},...,x_{n}^{\star})$ such
    that $\sum_{i=1}^{n}x_{i}^{\star}=1$, two methods are offered:

    \- Standard: $x_{i}^{\star}=\frac{x_{i}}{\sum_{i=1}^{n}x_{i}}$.

    \- Softmax:
    $x_{i}^{\star}=\frac{e^{x_{i}}}{\sum_{i=1}^{n}e^{x_{i}}}$ where
    $e\approx 2.72$ is Euler's constant.

<a name="similarity-measures"></a>
## 2.2 Similarity Measures
Given a pair of processed spectra intensities
$I=(a_{1},a_{2},...,a_{n}), J=(b_{1},b_{2},...,b_{n})\in\mathbb{R}^{n}$
with $0\leq a_{i},b_{i}\leq 1$ for all $i\in\{1,2,...,n\}$ and
$\sum_{i=1}^{n}a_{i}=\sum_{i=1}^{n}b_{i}=1$, MZsearch provides
functionality for computing the following similarity measures:

-   Cosine Similarity Measure:
      $S_{Cosine}(I,J)=\frac{I\cdot J}{|I|_{2}\cdot |J|_{2}}$ where multiplication in the numerator refers to $\frac{1}{2}$ the dot product $I\circ J=a_{1}b_{1}+a_{2}b_{2}+...+a_{n}b_{n}$ of $I$ and $J$ and multiplication in the denominator refers to multiplication of the $L^{2}$-norm of $I$ and $J$, $|I|_{2}=\sqrt{a_{1}^{2}+a_{2}^{2}+...+a_{n}^{2}}, |J|_{2}=\sqrt{b_{1}^{2}+b_{2}^{2}+...+b_{n}^{2}}$.

-   Shannon Entropy Similarity Measure:
      $S_{Shannon}(I,J) = 1-\frac{2\cdot H_{Shannon}\left(\frac{I+J}{2}\right) - H_{Shannon}(I)-H_{Shannon}(J)}{ln(4)}$
      $H_{Shannon}(I)=-\sum_{i=1}^{n}a_{i}\cdot ln(a_{i})$

-    Tsallis Entropy Similarity Measure:
      $S_{Tsallis}(I,J,q)=1-\frac{2\times H_{Tsallis}(I/2+J/2,q)-H_{Tsallis}(I,q)-H_{Tsallis}(J,q)}{N_{Tsallis}(I,J,q)},\\
      N_{Tsallis}(I,J,q):=\frac{\sum_{i=1}^{n}\left(2\left(\frac{a_{i}}{2}\right)^{q}+2\left(\frac{b_{i}}{2}\right)^{q}-a_{i}^{q}-b_{i}^{q}\right)}{1-q},\\
      H_{Tsallis}(I,q)=\frac{\left(\sum_{i=1}^{n}a_{i}^{q}\right)-1}{1-q},\\
      q\neq 1, \ q\textgreater 0$

-   Rényi Entropy Similarity Measure:
      $S_{Renyi}(I,J,q)=1-\frac{2\times H_{Renyi}(I/2+J/2,q)-H_{Renyi}(I,q)-H_{Renyi}(J,q)}{N_{Renyi}(I,J,q)},\\
      N_{Renyi}(I,J,q):=\left(\frac{1}{1-q}\right)\left(2\times ln\left(\sum_{i}(a_{i}/2)^{q}+\sum_{j}(b_{j}/2)^{q}\right)-ln(\sum_{i}a_{i}^{q})-ln(\sum_{i}b_{i}^{q})\right),\\
      H_{Renyi}(I,q)=\frac{1}{1-q}ln(\sum_{i=1}^{n}a_{i}^{q}),\\
      q\neq 1, \ q\textgreater 0$



<a name="usage"></a>
## 3. Usage
This repository has three main capabilities:
1. converting the raw data to the necessary format for spectral library matching
2. running spectral library matching to identify compounds based off of their mass spectrometry data
3. plotting a query spectrum vs a reference spectrum before and after preprocessing transformations.

These tasks are implemented separately for the cases of (i) GC-MS and (ii) LC-MS data due to the different spectrum preprocessing transformations stemming from a different format in the mass to charge ratios in GC-MS vs LC-MS data. To see all parameters for any of the four main scripts (spec_lib_matching_lcms.py, spec_lib_matching_gcms.py, plot_spectra_lcms.py, plot_spectra_gcms.py), run:
```
python build_library.py -h
python spec_lib_matching_.py -h
python spec_lib_matching_.py -h
python plot_spectra_.py -h
python plot_spectra_gcms.py -h
```

<a name="process-data"></a>
### 3.1 Obtain LC-MS or GC-MS library from MGF, mzML, or cdf file
To obtain a CSV file of LC-MS spectra in the format necessary for spectral library matching from raw data in the form of an mgf, mzML, or cdf file, one can run:
```
python build_library.py \
  --input_path path_to_input_mgf_or_mzML_or_cdf_file \
  --output_path path_to_output_csv_file
```

Parameter descriptions are as follows:

--input_path: Path to input file (must be either mgf, mlMZ, or cdf file). Mandatory argument.

--output_path: Path to output CSV file. Default: current working directory.

Some example MGF files one can use to build a LC-MS library can be found from the Global Natural Products Social Molecular Networking databases here: [https://external.gnps2.org/gnpslibrary](https://external.gnps2.org/gnpslibrary). Some example mzML file one can use to build a LC-MS library can be found in this repository: [https://github.com/HUPO-PSI/mzML](https://github.com/HUPO-PSI/mzML).

<a name="run-spec-lib-matching"></a>
### 3.2 Run spectral library matching
To run spectral library matching on LC-MS/GC-MS data, one can use:
```
python spec_lib_matching_.py \
  --query_data path_to_query__CSV_file \
  --reference_data path_to_reference__CSV_file

python spec_lib_matching_gcms.py \
  --query_data path_to_query_gcms_CSV_file \
  --reference_data path_to_reference_gcms_CSV_file
```

Example implementations of these scripts with all parameters specified are:
```
python spec_lib_matching_.py \
  --query_data path_to_query__CSV_file \
  --reference_data path_to_reference__CSV_file \
  --likely_reference_IDs None \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 999999999999 \
  --int_min 0 \
  --int_max 999999999999 \
  --window_size_centroiding 0.5 \
  --window_size_matching 0.5 \
  --noise_threshold 0 \
  --wf_mz 0 \
  --wf_intensity 1 \
  --LET_threshold 0 \
  --entropy_dimension 1.1 \
  --normalization_method standard \
  --n_top_matches_to_save 1 \
  --print_id_results False \
  --output_identification path_to__identification_results_CSV \
  --output_similarity_scores path_to_CSV_of_all__similarity_scores

python spec_lib_matching_gcms.py \
  --query_data path_to_query_gcms_CSV_file \
  --reference_data path_to_reference_gcms_CSV_file \
  --likely_reference_IDs None \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FNLW \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 999999999999 \
  --int_min 0 \
  --int_max 999999999999 \
  --wf_mz 0 \
  --wf_intensity 1 \
  --entropy_dimension 1.1 \
  --normalization_method standard \
  --n_top_matches_to_save 1 \
  --print_id_results False \
  --output_identification path_to_gcms_identification_results_CSV \
  --output_similarity_scores path_to_CSV_of_all_gcms_similarity_scores
```

To implement an example of spectral library matching, you can navigate to the scripts directory and run the following commands which should produce CSV files of the identification results in your current working directory (i.e. the scripts directory) by default:
```
python spec_lib_matching_lcms.py \
  --query_data "$PWD"/../data/lcms_query_library.csv \
  --reference_data "$PWD"/../data/lcms_reference_library.csv

python spec_lib_matching_gcms.py \
  --query_data "$PWD"/../data/gcms_query_library.csv \
  --reference_data "$PWD"/../data/gcms_reference_library.csv
```


Parameter descriptions are as follows:

--query_data: 
  * LC-MS case: 3-column CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a single ion fragment of a mass spectrum, the left-most column should contain an identifier, the middle columns should correspond the mass to charge ratios, and the right-most column should contain the intensities. For example, if spectrum A has 3 ion fragments, then there would be three rows in this CSV file corresponding to spectrum A. Default: LC-MS GNPS library.
  * GC-MS case: CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns contains the intensity with respect to a single mass to charge ratio. Default: GC-MS NIST WebBook library

--reference_data: same format CSV file as query_data except of reference library spectra.

--likely_reference_IDs: CSV file with one column containing the IDs of a subset of all compounds in the reference_data to be used in spectral library matching. Each ID in this file must be an ID in the reference library. Default: None (i.e. default is to use entire reference library)

--similarity_measure: options are 'cosine', 'shannon', 'renyi', and 'tsallis'.

--spectrum_preprocessing_order: The spectrum preprocessing transformations and the order in which they are to be applied. These transformations are applied prior to computing similarity scores. Format must be a string with 2-6 (LC-MS) or 2-4 (GC-MS) characters chosen from F, N, W, C, M, L representing filtering, noise removal, weight-factor-transformation, centroiding, matching, and low-entropy tranformation, respectively. Matching (M) and centroiding (C) are applicable only to LC-MS data. For example, if \'WCM\' is passed, then each (LC-MS) spectrum will undergo a weight factor transformation, then cleaning, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of  data. Default: FCNMWL for LC-MS and FNLW for GC-MS.

--high_quality_reference_library: True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False.

--mz_min: Remove all peaks with mass to charge less than mz_min in each spectrum. Default = 0.

--mz_max: Remove all peaks with mass to charge greater than mz_max in each spectrum. Default = 999999999999.

--int_min: Remove all peaks with intensity less than int_min in each spectrum. Default = 0.

--int_max: Remove all peaks with intensity greater than int_max in each spectrum. Default = 999999999999.

--window_size_centroiding (LC-MS only): Window size parameter used in centroiding a given spectrum. Default = 0.5.

--window_size_matching (LC-MS only): Window size parameter used in matching a query spectrum and a reference library spectrum. Default = 0.5.

--noise_threshold: Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default = 0.

--wf_mz: Mass to charge weight factor parameter. Default = 0.

--wf_intensity: Intensity weight factor parameter. Default = 1.

--LET_threshold: Low-entropy transformation threshold parameter. Spectra with Shannon entropy H less than LET_threshold are transformed according to $\text{intensitiesNew}=\text{intensitiesOriginal}^{\frac{1+S}{1+\text{LETthreshold}}}$. Default = 0.

--entropy_dimension: Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default = 1.1.

--normalization_method: Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: 'standard' and 'softmax'. Default = standard.

--n_top_matches_to_save: The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default = 1.

--print_id_results: Flag indicating whether to print the identification results interactively. Regardless of this flag, the results are saved according to the parameter 'output_identification'. Default = False.
--output_identification: Output CSV file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory (i.e. same directory this script is contained in) with filename 'output_lcms_identification.csv'.

--output_similarity_scores: Output CSV file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this CSV file is written to the current working directory with filename 'output_lcms_all_similarity_scores'.csv.

<a name="plotting"></a>
### 3.3 Plot a query spectrum against a reference spectrum before and after spectrum preprocessing transformations
To plot a query spectrum vs a reference spectrum before and after preprocessing transformations, run:
```
python plot_spectra_lcms.py \
  --query_data path_to_query_lcms_CSV_file \
  --reference_data path_to_reference_lcms_CSV_file \
  --query_spectrum_ID insert_single_ID_from_first_column_of_query_data \
  --reference_spectrum_ID insert_single_ID_from_first_column_of_reference_data \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 999999999999 \
  --int_min 0 \
  --int_max 999999999999 \
  --window_size_centroiding 0.5 \
  --window_size_matching 0.5 \
  --noise_threshold 0 \
  --wf_mz 0 \
  --wf_intensity 1 \
  --LET_threshold 0 \
  --entropy_dimension 1.1 \
  --normalization_method standard \
  --save_plots path_to_output_PDF_file

python plot_spectra_gcms.py \
  --query_data path_to_query_gcms_CSV_file \
  --reference_data path_to_reference_gcms_CSV_file \
  --query_spectrum_ID insert_single_ID_from_first_column_of_query_data \
  --reference_spectrum_ID insert_single_ID_from_first_column_of_reference_data \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FNLW \
  --high_quality_reference_library False \
  --mz_min 0 \
  --mz_max 999999999999 \
  --int_min 0 \
  --int_max 999999999999 \
  --wf_mz 0 \
  --wf_intensity 1 \
  --LET_threshold 0 \
  --entropy_dimension 1.1 \
  --normalization_method standard \
  --save_plots path_to_output_PDF_file
```

To implement an example of plotting spectra, you can navigate to the scripts directory and run the following commands which should produce PDF files of the spectra plots in your current working directory (i.e. the scripts directory) by default:
```
python plot_spectra_lcms.py \
  --query_data "$PWD"/../data/lcms_query_library.csv \
  --reference_data "$PWD"/../data/lcms_reference_library.csv \
  --query_spectrum_ID 100 \
  --reference_spectrum_ID 'Hectochlorin M+H' \
  --noise_threshold 0.05 \
  --wf_mz 0.5 \
  --wf_int 1.2

python plot_spectra_gcms.py \
  --query_data "$PWD"/../data/gcms_query_library.csv \
  --reference_data "$PWD"/../data/gcms_reference_library.csv \
  --query_spectrum_ID ID_1 \
  --reference_spectrum_ID 616386 \
  --similarity_measure tsallis
```


Parameter descriptions are as follows:

--query_data: 
  * LC-MS case: 3-column CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a single ion fragment of a mass spectrum, the left-most column should contain an identifier, the middle columns should correspond the mass to charge ratios, and the right-most column should contain the intensities. For example, if spectrum A has 3 ion fragments, then there would be three rows in this CSV file corresponding to spectrum A. Default: LC-MS GNPS library.
  * GC-MS case: CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns contains the intensity with respect to a single mass to charge ratio. Default: GC-MS NIST WebBook library

--reference_data: Same format CSV file as query_data except of reference library spectra.

--query_spectrum_ID: The identifier of the query spectrum to be plotted. Default: first query spectrum in query_data.

--reference_spectrum_ID: The identifier of the reference spectrum to be plotted. Default: first reference spectrum in reference_data.

--similarity_measure: Options are 'cosine', 'shannon', 'renyi', and 'tsallis'.

--spectrum_preprocessing_order: The LC-MS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-4 characters chosen from W, C, M, L representing weight-factor-transformation, cleaning (i.e. centroiding and noise removal), matching, and low-entropy transformation. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then cleaning, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of LC-MS data. Default: FCNMWL for LC-MS and FNLW for GC-MS .

--high_quality_reference_library: True/False flag indicating whether the reference library is considered to be of high quality. If True, then the spectrum preprocessing transformations of filtering and noise removal are performed only on the query spectrum/spectra. If False, all spectrum preprocessing transformations specified will be applied to both the query and reference spectra. Default: False.

--mz_min: Remove all peaks with mass to charge less than mz_min in each spectrum. Default = 0.

--mz_max: Remove all peaks with mass to charge greater than mz_max in each spectrum. Default = 999999999999.

--int_min: Remove all peaks with intensity less than int_min in each spectrum. Default = 0.

--int_max: Remove all peaks with intensity greater than int_max in each spectrum. Default = 999999999999.

--window_size_centroiding (LC-MS only): Window size parameter used in centroiding a given spectrum. Default = 0.5.

--window_size_matching (LC-MS only): Window size parameter used in matching a query spectrum and a reference library spectrum. Default = 0.5.

--noise_threshold: Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default = 0.

--wf_mz: Mass to charge weight factor parameter. Default = 0.

--wf_intensity: Intensity weight factor parameter. Default = 1.

--LET_threshold: Low-entropy transformation threshold parameter. Spectra with Shannon entropy H less than LET_threshold are transformed according to $\text{intensitiesNew}=\text{intensitiesOriginal}^{\frac{1+S}{1+\text{LETthreshold}}}$. Default = 0.

--entropy_dimension: Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default = 1.1.

--normalization_method: Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: 'standard' and 'softmax'. Default = standard.

--save_plots: Output PDF file containing the plots of the query and reference spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./query_spec_{query_spectrum_ID}_reference_spec_{reference_spectrum_ID}_plot.pdf in the current working directory.

<a name="bugs-questions"></a>
## 4. Bugs/Questions?
If you notice any bugs in this software or have any questions, don't hesitate to reach out to fy7392@wayne.edu.






