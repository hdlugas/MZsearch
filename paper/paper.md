---
title: 'MZsearch: A Python-Based Compound Identification Tool for GC-MS and LC-MS/MS-Based Metabolomics'
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
  - entropy
authors:
  - name: Hunter Dlugas
    orcid: 0000-0002-6819-0045
    affiliation: "1, 2"
  - name: Xiang Zhang
    orcid: 0000-0003-1102-6313
    affiliation: "3"
  - name: Seongho Kim
    orcid: 0000-0003-1120-073X
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

A primary goal in metabolomics, the qualitative and quantitative study of metabolites, is to identify chemical compounds in a given sample. This often involves identifying and quantifying metabolic biomarkers for use in diagnosing, treating, and stratifying the risk of various diseases. A central tool for compound identification in metabolomics is mass spectrometry, which produces the mass spectrum for a given chemical compound. Each peak in the mass spectrum represents an ion fragment, with one axis representing mass-to-charge ratio (m/z) and the other axis representing relative intensity. A popular method used to identify a compound based on its mass spectrum is spectral library matching. This process involves computing a measure of similarity between the mass spectrum of the compound and each mass spectrum in a reference library. The compound is then identified as the one from the reference library whose mass spectrum is most similar to that of the compound. We present MZsearch, a Python-based tool for spectral library matching, available in two versions: a command-line interface and Python modules for integration into custom code. The MZsearch tool offers a range of spectrum preprocessing transformations and similarity measures and can identify compounds from nominal-resolution mass spectra (NR-MS) and high-resolution mass spectra (HR-MS).

# Statement of need

Several open-source software tools have been developed for analyzing mass spectrometry data. The Python package matchms allows one to perform spectral library matching for both GC-MS and LC-MS/MS and provides an application programming interface (API) for users to define their own spectrum preprocessing transformations and similarity measures [@Huber2020]. The pyOpenMS tool, which is a Python wrapper of the OpenMS tool, is for analyzing LC-MS/MS data and provides common spectrum preprocessing transformations such as filtering based on m/z and intensity values [@Rost2014; @Rost2016]. The Python package spectrum_utils provides functionality for preprocessing and visualizing mass spectrometry data [@Bittremieux2020]. The R package metID is a tool for LC-MS/MS-based compound identification with a focus on allowing users to combine their in-house database(s) with public databases [@Shen2022]. The R Shiny package ShinyMetID provides users six similarity measures (Cosine, Weighted Cosine, Stein & Scott, Discrete Fourier Transformation, Discrete Wavelet Transformation, and Semi-Partial Correlation) and a graphical user interface (GUI) to perform spectral library matching to identify chemical compounds from GC-MS data [@Jeong2023].
In many GC-MS experiments, the resulting mass spectra often have m/z values reported as integers due to the use of nominal mass resolution of mass spectrometers. In contrast, LC-MS/MS frequently employs a high-resolution mass spectrometer that provides m/z values with multiple decimal places. This necessitates several spectrum preprocessing steps unique to the high-resolution MS or MS/MS data analysis, namely centroiding (i.e., merging peaks that are 'close' with respect to their m/z values) and matching (i.e., aligning the m/z values so that the query spectrum and reference spectrum have the same length). In addition to these canonical spectrum preprocessing transformations for high-resolution MS/MS data, weight factor transformations and low-entropy transformations have been proposed to improve the performance of compound identification for both GC-MS and LC-MS/MS data [@Kim2012; @Li2021; @Dlugas2024_preprint]. The Shannon Entropy Similarity Measure has been shown to outperform the Cosine Similarity Measure with respect to LC-MS/MS data [@Li2021]. A generalization of the Shannon Entropy Similarity Measure, the Tsallis Entropy Similarity Measure, slightly outperforms the Shannon Entropy Similarity Measure in both GC-MS and LC-MS/MS data [@Dlugas2024_preprint]. This recent study has further demonstrated that the order of preprocessing transformations is also critical to achieving higher accuracy in compound identification. However, to our knowledge, there is no tool that includes the recently-developed, high-performance entropy-based similarity measures or considers the order of the spectrum preprocessing transformations, which are important for effective metabolomics studies.
To address the lack of spectral library matching software that considers both the order of spectrum preprocessing steps and novel entropy-based similarity measures, MZsearch was developed. The developed MZsearch is a Python-based tool for performing spectral library matching on both nominal-resolution mass spectrometry (NR-MS) data (e.g. GC-MS) and high-resolution mass spectrometry  (HR-MS) data (e.g., LC-MS/MS). It is available in two versions: a command-line interface and Python modules for integration into custom code. Users can customize their own spectrum preprocessing order using spectrum preprocessing transformations such as weight factor and low-entropy transformations. Additionally, MZsearch includes a novel entropy-based similarity measure, the Rényi Entropy Similarity Measure, enabling users to choose among four similarity measures: the commonly-used Cosine Similarity Measure [@Stein1994], the Shannon Entropy Similarity Measure [@Li2021], the recently-developed Tsallis Entropy Similarity Measure [@Dlugas2024_preprint], and the novel Rényi Entropy Similarity Measure. **Table 1** compares MZsearch with other similar software.

**Table  1.** Comparison of MZsearch and existing tools.

![\label{comparison}](table.PNG){width="100%," height="100%"}\

# Functionality

## Spectrum Preprocessing Transformations

Functionality implementing the following spectrum preprocessing
transformations is offered in MZsearch (more details found in documentation: <a href="[[url](https://github.com/hdlugas/MZsearch)](https://github.com/hdlugas/MZsearch)"></a>):

-   Filtering based on mz-values and intensities

-   Weight Factor Transformation

-   Low-Entropy Transformation

-   Noise Removal

-   Centroiding (only applicable to LC-MS/MS data)

-   Matching (only applicable to LC-MS/MS data)

-   Normalization (only applicable to entropy-based similarity measures)

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
    dot product $I\circ J=a_{1}b_{1}+a_{2}b_{2}+...+a_{n}b_{n}$ and multiplication in the denominator refers to
    multiplication of the $L^{2}$-norm of $I$ and $J$,
    $|I|_{2}=\sqrt{a_{1}^{2}+a_{2}^{2}+...+a_{n}^{2}}, |J|_{2}=\sqrt{b_{1}^{2}+b_{2}^{2}+...+b_{n}^{2}}$.

-   Shannon Entropy Similarity Measure: \begin{gather*}
      S_{Shannon}(I,J) = 1-\frac{2\cdot H_{Shannon}\left(\frac{I+J}{2}\right) - H_{Shannon}(I)-H_{Shannon}(J)}{ln(4)},\\
      H_{Shannon}(I)=-\sum_{i=1}^{n}a_{i}\cdot ln(a_{i})
    \end{gather*}

-    Tsallis Entropy Similarity Measure [@Tsallis1988; @Havrda1967]:
    \begin{gather*}\label{eq:}
      S_{Tsallis}(I,J,q)=1-\frac{2\times H_{Tsallis}(I/2+J/2,q)-H_{Tsallis}(I,q)-H_{Tsallis}(J,q)}{N_{Tsallis}(I,J,q)},\\
      N_{Tsallis}(I,J,q):=\frac{\sum_{i=1}^{n}\left(2\left(\frac{a_{i}}{2}\right)^{q}+2\left(\frac{b_{i}}{2}\right)^{q}-a_{i}^{q}-b_{i}^{q}\right)}{1-q},\\
      H_{Tsallis}(I,q)=\frac{\left(\sum_{i=1}^{n}a_{i}^{q}\right)-1}{1-q},\\
      q\neq 1, \ q > 0
    \end{gather*}

-   Rényi Entropy Similarity Measure: \begin{gather*}\label{eq:renyi}
      S_{Renyi}(I,J,q)=1-\frac{2\times H_{Renyi}(I/2+J/2,q)-H_{Renyi}(I,q)-H_{Renyi}(J,q)}{N_{Renyi}(I,J,q)},\\
      N_{Renyi}(I,J,q):=\left(\frac{1}{1-q}\right)\left(2\times ln\left(\sum_{i}(a_{i}/2)^{q}+\sum_{j}(b_{j}/2)^{q}\right)-ln(\sum_{i}a_{i}^{q})-ln(\sum_{i}b_{i}^{q})\right),\\
      H_{Renyi}(I,q)=\frac{1}{1-q}ln(\sum_{i=1}^{n}a_{i}^{q}),\\
      q\neq 1, \ q > 0
    \end{gather*}

# Usage

Compound identification occurs during preprocessing or as an independent procedure. In preprocessing, it follows peak detection and should be performed before cross-sample alignment if based on compound names; otherwise, it can be delayed until before pathway analysis. As an independent procedure, it can be conducted ad hoc or post hoc to confirm or identify specific peaks in single or multiple samples. MZsearch performs compound identification and has three main capabilities: (i) converting the raw data to the necessary format for spectral library matching, (ii) running spectral library matching to identify compounds based on their mass spectrometry data, and (iii) plotting a query spectrum vs. a reference spectrum before and after preprocessing transformations. These tasks are implemented separately for (i) NR-MS data such as GC-MS and (ii) HR-MS data such as LC-MS/MS data due to the different spectrum preprocessing transformations. To see all parameters for any of the three main scripts (build_library.py, spec_lib_matching.py, plot_spectra.py), see the documentation at [https://github.com/hdlugas/MZsearch](https://github.com/hdlugas/MZsearch).


# Acknowledgements

This research was partially funded by NIH grants (National Institutes of
Health R21GM140352), (National Cancer Institute P30CA022453), (National
Institutes of Health P20GM113226) and (National Institutes of Health
P50AA024337). Computations were performed on Wayne State University's
High-Performance Computing Grid.


# Author Contributions

S.K. conceived the study and contributed to the overall study design and interpretation. H.D., X.Z., and S.K. collaborated on discussing the design, features, and structures of the developed software package. H.D. implemented the initial idea, further expanded and refined it, and was responsible for coding the developed software package. H.D., X.Z., and S.K. reviewed and edited the manuscript. Both H.D. and S.K. contributed equally to the manuscript preparation, including drafting and revising the manuscript.

# References
