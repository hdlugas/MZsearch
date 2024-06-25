---
title: 'hello'
tags:
  - Python
  - metabolomics
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
  - name: Seongho Kim
    equal-contrib: true 
    affiliation: "1, 3"
affiliations:
 - name: Wayne State University School of Medicine, USA
   index: 1
 - name: Biostatistics and Bioinformatics Core, Karmanos Cancer Institute
   index: 2
 - name: Biostatistics and Bioinformatics Core, Karmanos Cancer Institute/Department of Oncology
   index: 3
date: 25 June 2024
bibliography: paper.bib
```


# Summary

A primary goal of the field of Metabolomics - i.e. the qualitative and quantitative study of metabolites - is to identify chemical compounds in a given sample, oftentimes with the motivation of identifying/quantifying biomarkers for use in diagnosing, treating, and/or stratifying the risk of some disease. A central tool for compound identification in Metabolomics is mass spectrometry which produces a finite set of points - termed the mass spectrum - in the plane $\mathbb{R}^{2}$ for a given chemical compound with each point representing an ion fragment, one axis representing mass/charge, and the other axis representing intensity. The primary method used to identify a given unknown chemical compound based off of its mass spectrum is termed 'spectral library matching' and involves computing a measure of similarity between the mass spectrum of the unknown chemical compound and each mass spectrum of chemical compounds in a predetermined reference library. The unknown chemical compound is then identified as the chemical compound from the reference library whose mass spectrum is most similar to the mass spectrum of the unknown chemical compound. We present (insert package name), a command-line tool written in Python implementing spectral library matching with a host of spectrum preprocessing transformations and similarity measures available in addition to being able to analyze both of the two main types of mass spectrometry data: gas chromatography-mass spectrometry (GC-MS) and liquid chromatography-mass spectrometry (LC-MS) data.


# Statement of need

hello world citation here 

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission by students and experts alike.



# Spectrum Preprocessing Transformations
Functionality implementing the following spectrum preprocessing transformations is offered in (insert package name):

* Weight Factor Transformation: Given a pair of user-defined weight factor parameters $(\text{a,b})$ and spectrum $I$ with mass/charge values $(m_{1},m_{2},...,m_{n})\in\mathbb{R}^{n}$ and intensities $(x_{1},x_{2},...,x_{n})\in\mathbb{R}^{n}$, the transformed spectrum $I^{\star}$ has the same mass/charge values as $I$ and has intensities given by $I^{\star}:=(m_{1}^{\text{a}}\cdot x_{1}^{\text{b}},m_{2}^{\text{a}}\cdot x_{2}^{\text{b}},...,m_{n}^{\text{a}}\cdot x_{n}^{\text{b}}$.

* Low-Entropy Transformation: Given a user-defined low-entropy threshold parameter $T$ and spectrum $I$ with intensities $(x_{1},x_{2},...,x_{n})\in\mathbb{R}^{n}$, then the transformed spectrum intensities $I^{\star}=(x_{1}^{\star},x_{2}^{\star},...,x_{n}^{\star})$ are such that, for all $i\in\{1,2,...,n\}$, $x_{i}^{\star}=x_{i}$ if $H_{Shannon}\geq T$ and $x_{i}^{\star}=x_{i}^{\frac{1+H_{Shannon}(I)}{1+T}}$ if $H_{Shannon}\textless T$.

* Cleaning: Given a user-defined window-size parameter $w_{centroiding}$, a user-defined noise removal parameter $r$, and a spectrum $I$ with mass/charge values $(m_{1},m_{2},...,m_{n})\in\mathbb{R}^{n}$ and intensities $(x_{1},x_{2},...,x_{n})\in\mathbb{R}^{n}$, the transformed spectrum $I^{\star}$ merges adjacent points $(m_{i},x_{i}),(m_{i+1},x_{i+1})$, $i\in\{1,2,...,n-1\}$ if $|m_{i}-m_{i+1}|\textless w_{centroiding}$ into the point $(\frac{m_{i}\cdot x_{i}+m_{i+1}\cdot x_{i+1}}{x_{i}+x_{i+1}},x_{i}+x_{i+1})$. This centroiding procedure generalizes to more than two points whose mass/charge values are within distance $w_{centroiding}$ of each other. If we denote the intensities of $I^{\star}$ as $(y_{1},y_{2},...,y_{m}), noise removal then removes points from $I^{\star}$ with $y_{j}\textless \text{max}(\{y_{1},y_{2},...,y_{m}\})$ for all $j\in\{1,2,...,m\}$.

* Matching (only for LCMS data): Given a user-defined window-size parameter $w_{matching}$ and two spectra $I$, $J$ with mass/charge ratios $(a_{1},a_{2},...,a_{n}), (b_{1},b_{2},...,b_{m})$ and intensities $(x_{1},x_{2},...,x_{n}), (y_{1},y_{2},...,y_{m}$, respectively, of which we would like to measure the similarity between, the matching procedure outputs two spectra $I^{\star},J^{\star}$ containing the same number of points with $I^{\star}$ and $J^{\star}$ having identical mass/charge ratios and transformed intensities. Specifically,for a given point $(a_{i},x_{i})$ of $I$, if there are no points $(b_{j},y_{j})$ in $J$ with $|a_{i}-b_{j}|\textless w_{matching}$, then the point $(a_{i},x_{i})$ remains in $I^{\star}$ and the point $(a_{i},0)$ is included in $J^{\star}$. If there is at least one point $(b_{j},y_{j})$ with $|a_{i}-b_{j}|\textless w_{matching}$, then the point $(a_{i},x_{i})$ remains in $I^{\star}$ and the point $(a_{i},\sum_{j}b_{j})$ is included in $J^{\star}$. This procedure is applied when transposing the roles of $I$ and $J$ as well.

* Normalization: Prior to computing entropy - regardless of whether in the context of performing the low-entropy transformation or computing similarity score - the intensities of each spectrum must be normalized to sum to 1 in order to represent a probability distribution. To normalize a given spectrum $I$ with intensities $(x_{1},x_{2},...,x_{n})$ into spectrum $I^{\star}$ with intensities $(x_{1}^{\star},x_{2}^{\star},...,x_{n}^{\star})$ such that $\sum_{i=1}^{n}x_{i}^{\star}=1$, two methods are offered:

  ** Standard: $x_{i}^{\star}=\frac{x_{i}}{\sum_{i=1}^{n}x_{i}}$
  
  ** Softmax: $x_{i}^{\star}=\frac{e^{x_{i}}}{\sum_{i=1}^{n}e^{x_{i}}}$ where $e\approx 2.72$ is Euler's constant.


# Similarity Measures
Given a pair of processed spectra intensities $I=(a_{1},a_{2},...,a_{n}), J=(b_{1},b_{2},...,b_{n})\in\mathbb{R}^{n}$, with $0\leq a_{i},b_{i}\leq 1$ for all $i\in\{1,2,...,n\}$ and $\sum_{i=1}^{n}a_{i}=\sum_{i=1}^{n}b_{i}=1$, (insert package name) provides functionality for computing the following similarity measures:

* Cosine Similarity Measure:
\begin{equation*}
    S_{Cosine}(I,J)= \frac{I\circ J}{|I|_{2}\cdot |J|_{2}}
\end{equation*}
where multiplication in the numerator refers to the dot product $I\circ J=a_{1}b_{1}+a_{2}b_{2}+...+a_{n}b_{n}$ of $I$ and $J$ and multiplication in the denominator refers to multiplication of the $L^{2}$-norm of $I$ and $J$, $|I|_{2}=\sqrt{a_{1}^{2}+a_{2}^{2}+...+a_{n}^{2}}, |J|_{2}=\sqrt{b_{1}^{2}+b_{2}^{2}+...+b_{n}^{2}}$.

* Shannon Entropy Similarity Measure:
\begin{gather*}
    S_{Shannon}(I,J) = 1-\frac{2\cdot H_{Shannon}\left(\frac{I+J}{2}\right) - H_{Shannon}(I)-H_{Shannon}(J)}{ln(4)},\\
    H_{Shannon}(I)=-\sum_{i=1}^{n}p_{i}\cdot ln(p_{i})
\end{gather*}

* Tsallis Entropy Similarity Measure:
\begin{gather*}\label{eq:tsallis}
    S_{Tsallis}(I_{q},I_{l},q)=1-\frac{2\times H_{Tsallis}(I_{Q}/2+I_{L}/2,q)-H_{Tsallis}(I_{Q},q)-H_{Tsallis}(I_{L},q)}{N_{Tsallis}},\\
    N_{Tsallis}:==\frac{\sum_{i=1}^{n}\left(2\left(\frac{a_{i}}{2}\right)^{q}+2\left(\frac{b_{i}}{2}\right)^{q}-a_{i}^{q}-b_{i}^{q}\right)}{1-q},\\
    H_{Tsallis}(I,q)=\frac{\left(\sum_{i=1}^{n}p_{i}^{q}\right)-1}{1-q},\\
    q\neq 1, \ q\textgreater 0
\end{gather*}

* R\'enyi Entropy Similarity Measure:
\begin{gather*}\label{eq:renyi}
    S_{R\'enyi}(I_{Q}, I_{L})=1-\frac{2\times H_{R\'enyi}(I_{Q}/2+I_{L}/2,q)-H_{R\'enyi}(I_{Q},q)-H_{R\'enyi}(I_{L},q)}{N_{R\'enyi}},\\
    N_{R\'enyi}:=\left(\frac{1}{1-q}\right)\left(2\times ln\left(\sum_{i}(a_{i}/2)^{q}+\sum_{j}(b_{j}/2)^{q}\right)-ln(\sum_{i}a_{i}^{q})-ln(\sum_{i}b_{i}^{q})\right),\\
    H_{R\'enyi}(I,q)=\frac{1}{1-q}ln(\sum_{i=1}^{n}p_{i}^{q}),\\
    q\neq 1, \ q\textgreater 0
\end{gather*}


You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"


# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }


# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.


# References
