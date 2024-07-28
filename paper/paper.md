---
title: 'MZsearch: A Python-Based Compound Idnetification Tool for GC-MS and LC-MS-Based Metabolomics'
tags:
  - Python
  - metabolomics
  - compound identification
  - similarity measure
  - spectral library matching
  - compound identification
  - Cosine correlation
  - Shannon entropy
  - Renyi entropy
  - Tsallis entropy
authors:
  - name: Hunter Dlugas
    orcid: 0000-0002-6819-0045
    equal-contrib: true
    affiliation: "1, 2"
  - name: Xiang Zhang
    orcid: 0000-0003-1102-6313
    affiliation: "3"
  - name: Seongho Kim
    equal-contrib: true 
    affiliation: "1, 2"
affiliations:
 - name: Department of Oncology, School of Medicine, Wayne State University, Detroit, MI, USA
   index: 1
 - name: Biostatistics and Bioinformatics Core, Karmanos Cancer Institute, Wayne State University, Detroit, MI, USA
   index: 2
 - name: Department of Chemistry, University of Louisville, Louisville, KY, USA
   index: 3
date: 25 June 2024
bibliography: paper.bib
editor_options: 
  markdown: 
    wrap: 72
---

# Summary

A primary goal in metabolomics, the qualitative and quantitative study
of metabolites, is to identify chemical compounds in a given sample.
This often involves identifying and quantifying metabolic biomarkers
for use in diagnosing, treating, and stratifying the risk of various
diseases. A central tool for compound identification in metabolomics is
mass spectrometry, which produces the mass spectrum for a given chemical
compound. Each peak in the mass spectrum represents an ion fragment,
with one axis representing mass to charge ratio (m/z) and the other axis
representing relative intensity. A popular method used to identify an
unknown chemical compound based on its mass spectrum is spectral library
matching. This process involves computing a measure of similarity
between the mass spectrum of the unknown chemical compound and each mass
spectrum in a reference library. The unknown chemical compound is then
identified as the one from the reference library whose mass spectrum is
most similar to that of the unknown compound. We present 'MZsearch', a
command-line tool written in Python that implements spectral library
matching. The 'MZsearch' tool offers a range of spectrum preprocessing
transformations and similarity measures, and it is capable of analyzing
data generated from gas
chromatography-mass spectrometry (GC-MS) and liquid chromatography-mass
spectrometry (LC-MS).

# Statement of need

There exist several open-source software tools for working with mass
spectrometry data. The Python package 'matchms' allows one to perform
spectral library matching for both GC-MS and LC-MS and provides an
application programming interface (API) for users to define their own
spectrum preprocessing transformations and similarity measures
[@Huber2020]. The 'pyOpenMS' tool, which is a Python wrapper of the
'OpenMS' tool, is for analyzing LC-MS data and provides common spectrum
preprocessing transformations such as filtering based on m/z and
intensity values [@Rost2014; @Rost2016]. The Python package
'spectrum_utils' provides functionality for preprocessing and
visualizing mass spectrometry data [@Bittremieux2020]. The R package
'metID' is a tool for LC-MS-based compound identification with a focus
on allowing users to combine their in-house database(s) with public
databases [@Shen2022]. The R Shiny package 'ShinyMetID' provides users
six similarity measures (Cosine, Weighted Cosine, Stein & Scott,
Discrete Fourier Transformation, Discrete Wavelet Transformation, and
Semi-Partial Correlation) and a graphical user interface (GUI) to
perform spectral library matching to identify chemical compounds from GC-MS data [@Jeong2023].

In a GC-MS experiment, the resulting mass spectra typically have m/z
values as integers. In contrast, LC-MS
experiments produce mass spectra with high resolution m/z values,
often with many decimal places. This necessitates several spectrum
preprocessing steps unique to LC-MS data analysis, namely centroiding
(i.e., merging peaks that are 'close' with respect to their m/z values)
and matching (i.e., aligning the m/z values so that the query spectrum and reference spectrum have the same length). In
addition to these canonical spectrum preprocessing transformations for
LC-MS data, weight factor transformations and low-entropy
transformations have been proposed to improve the performance of
compound identification for both GC-MS and LC-MS data [@Kim2012;
@Li2021; @Dlugas2024_preprint]. The Shannon Entropy Similarity Measure
has been shown to outperform the Cosine Similarity Measure with respect
to LC-MS-data [@Li2021]. A generalization of
the Shannon Entropy Similarity Measure, the Tsallis Entropy Similarity
Measure, slightly outperforms the
Shannon Entropy Similarity Measure in both GC-MS and LC-MS data
[@Dlugas2024_preprint]. This recent study has further demonstrated that
the order of preprocessing transformations is also critical to achieving
higher accuracy in compound identification. However, to our knowledge,
there is no tool that includes the recently-developed, high-performance
entropy-based similarity measures or considers the order of the spectrum
preprocessing transformations, which are important for effective
metabolomics studies.

The developed 'MZsearch' is a command-line tool for performing spectral
library matching on either GC-MS or LC-MS data. It allows users to
construct their own spectrum preprocessing order using spectrum
preprocessing transformations such as weight factor and low-entropy
transformations. Additionally, the 'MZsearch' tool includes a novel
entropy-based similarity measure, the Rényi Entropy Similarity Measure,
enabling users to choose among four similarity measures: the
commonly-used Cosine Similarity Measure [@Stein1994], the Shannon
Entropy Similarity Measure [@Li2021], the recently-developed Tsallis
Entropy Similarity Measure [@Dlugas2024_preprint], and the novel Rényi
Entropy Similarity Measure.

# Functionality

## Spectrum Preprocessing Transformations

Functionality implementing the following spectrum preprocessing
transformations is offered in 'MZsearch':

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

The flowchart in \autoref{fig:flowchart} depicts the overall workflow of
MZsearch.

![Workflow of
MZsearch.\label{fig:flowchart}](flowchart.png){width="100%,"
height="100%"}

## Similarity Measures

Given a pair of processed spectra intensities
$I=(a_{1},a_{2},...,a_{n}), J=(b_{1},b_{2},...,b_{n})\in\mathbb{R}^{n}$
with $0\leq a_{i},b_{i}\leq 1$ for all $i\in\{1,2,...,n\}$ and
$\sum_{i=1}^{n}a_{i}=\sum_{i=1}^{n}b_{i}=1$, MZsearch provides
functionality for computing the following similarity measures:

-   Cosine Similarity Measure: \begin{equation*}
      S_{Cosine}(I,J)= \frac{I\circ J}{|I|_{2}\cdot |J|_{2}}
    \end{equation*} where multiplication in the numerator refers to the
    dot product $I\circ J=a_{1}b_{1}+a_{2}b_{2}+...+a_{n}b_{n}$ of $I$
    and $J$ and multiplication in the denominator refers to
    multiplication of the $L^{2}$-norm of $I$ and $J$,
    $|I|_{2}=\sqrt{a_{1}^{2}+a_{2}^{2}+...+a_{n}^{2}}, |J|_{2}=\sqrt{b_{1}^{2}+b_{2}^{2}+...+b_{n}^{2}}$.

-   Shannon Entropy Similarity Measure: \begin{gather*}
      S_{Shannon}(I,J) = 1-\frac{2\cdot H_{Shannon}\left(\frac{I+J}{2}\right) - H_{Shannon}(I)-H_{Shannon}(J)}{ln(4)},\\
      H_{Shannon}(I)=-\sum_{i=1}^{n}a_{i}\cdot ln(a_{i})
    \end{gather*}

-   Tsallis Entropy Similarity Measure:
    \begin{gather*}\label{eq:tsallis}
      S_{Tsallis}(I,J,q)=1-\frac{2\times H_{Tsallis}(I/2+J/2,q)-H_{Tsallis}(I,q)-H_{Tsallis}(J,q)}{N_{Tsallis}(I,J,q)},\\
      N_{Tsallis}(I,J,q):=\frac{\sum_{i=1}^{n}\left(2\left(\frac{a_{i}}{2}\right)^{q}+2\left(\frac{b_{i}}{2}\right)^{q}-a_{i}^{q}-b_{i}^{q}\right)}{1-q},\\
      H_{Tsallis}(I,q)=\frac{\left(\sum_{i=1}^{n}a_{i}^{q}\right)-1}{1-q},\\
      q\neq 1, \ q\textgreater 0
    \end{gather*}

-   Rényi Entropy Similarity Measure: \begin{gather*}\label{eq:renyi}
      S_{Renyi}(I,J,q)=1-\frac{2\times H_{Renyi}(I/2+J/2,q)-H_{Renyi}(I,q)-H_{Renyi}(J,q)}{N_{Renyi}(I,J,q)},\\
      N_{Renyi}(I,J,q):=\left(\frac{1}{1-q}\right)\left(2\times ln\left(\sum_{i}(a_{i}/2)^{q}+\sum_{j}(b_{j}/2)^{q}\right)-ln(\sum_{i}a_{i}^{q})-ln(\sum_{i}b_{i}^{q})\right),\\
      H_{Renyi}(I,q)=\frac{1}{1-q}ln(\sum_{i=1}^{n}a_{i}^{q}),\\
      q\neq 1, \ q\textgreater 0
    \end{gather*}

# Usage

This 'MZsearch' has three main capabilities: (i) converting the raw data
to the necessary format for spectral library matching, (ii) running
spectral library matching to identify compounds based on their mass
spectrometry data, and (iii) plotting a query spectrum vs. a reference
spectrum before and after preprocessing transformations.

These tasks are implemented separately for (i) GC-MS and (ii) LC-MS data
due to the different spectrum preprocessing transformations. To see
all parameters for any of the five main scripts (build_library.py,
spec_lib_matching_lcms.py, spec_lib_matching_gcms.py,
plot_spectra_lcms.py, plot_spectra_gcms.py), run:

```         
python build_library.py -h
python spec_lib_matching_lcms.py -h
python spec_lib_matching_lcms.py -h
python plot_spectra_lcms.py -h
python plot_spectra_gcms.py -h
```

The user can use their own libraries, and by default, 'MZsearch'
provides GC-MS and LC-MS reference libraries that have been trimmed. 
The recommended GC-MS and LC-MS reference libraries are available at the Zenodo 
database [@Dlugas2024_database]. The GC-MS library was extracted from 
the National Institute of Standards and Technology WebBook, while the LC-MS 
library was obtained from Global Natural Products Social Molecular Networking.

## Get query/reference library from mgf, mzML, or cdf file

To obtain a query/reference library in the necessary format for spectral
library matching from an mgf, mzML, or cdf file, one can use:

```         
python build_library.py \
  --input_path path_to_input_mgf_or_mzML_or_cdf_file \
  --output_path path_to_output_csv_file 
```

## Run spectral library matching

To run spectral library matching on LC-MS/GC-MS data, one can use:

```         
python spec_lib_matching_lcms.py \
  --query_data path_to_query_lcms_CSV_file \
  --reference_data path_to_reference_lcms_CSV_file \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0\
  --mz_max 999999999999\
  --int_min 0\
  --int_max 999999999999\
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
  --output_identification path_to_lcms_identification_results_CSV \
  --output_similarity_scores path_to_CSV_of_all_lcms_similarity_scores

python spec_lib_matching_gcms.py \
  --query_data path_to_query_gcms_CSV_file \
  --reference_data path_to_reference_gcms_CSV_file \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FNLW \
  --high_quality_reference_library False \
  --mz_min 0\
  --mz_max 999999999999\
  --int_min 0\
  --int_max 999999999999\
  --wf_mz 0 \
  --wf_intensity 1 \
  --entropy_dimension 1.1 \
  --normalization_method standard \
  --n_top_matches_to_save 1 \
  --print_id_results False \
  --output_identification path_to_gcms_identification_results_CSV \
  --output_similarity_scores path_to_CSV_of_all_gcms_similarity_scores
```

## Plot a query spectrum against a reference spectrum before and after spectrum preprocessing transformations

To plot a query spectrum vs. a reference spectrum before and after
preprocessing transformations, run:

```         
python plot_spectra_lcms.py \
  --query_data path_to_query_lcms_CSV_file \
  --reference_data path_to_reference_lcms_CSV_file \
  --likely_reference_IDs path_to_likely_reference_IDs_CSV_file \
  --query_spectrum_ID insert_single_ID_from_first_column_of_query_data \
  --reference_spectrum_ID insert_single_ID_from_first_column_of_reference_data \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FCNMWL \
  --high_quality_reference_library False \
  --mz_min 0\
  --mz_max 999999999999\
  --int_min 0\
  --int_max 999999999999\
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
  --likely_reference_IDs path_to_likely_reference_IDs_CSV_file \
  --query_spectrum_ID insert_single_ID_from_first_column_of_query_data \
  --reference_spectrum_ID insert_single_ID_from_first_column_of_reference_data \
  --similarity_measure cosine \
  --spectrum_preprocessing_order FNLW \
  --high_quality_reference_library False \
  --mz_min 0\
  --mz_max 999999999999\
  --int_min 0\
  --int_max 999999999999\
  --wf_mz 0 \
  --wf_intensity 1 \
  --LET_threshold 0 \
  --entropy_dimension 1.1 \
  --normalization_method standard \
  --save_plots path_to_output_PDF_file
```

An example output figure is shown in (\autoref{fig:spectra_plot}).

![Query spectrum plotted against reference spectrum before/after
spectrum preprocessing transformations. A low-entropy transformation is
applied to the LC-MS spectra, and a weight factor tranformation is
applied to the GC-MS spectra.
\label{fig:spectra_plot}](spectra_figure.png){width="100%,"
height="100%"}

# Acknowledgements

This research was partially funded by NIH grants (National Institutes of
Health R21GM140352), (National Cancer Institute P30CA022453), (National
Institutes of Health P20GM113226) and (National Institutes of Health
P50AA024337). Computations were performed on Wayne State University's
High-Performance Computing Grid.

# References
