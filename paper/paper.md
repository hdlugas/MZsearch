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
---

# Summary

A primary goal of the field of Metabolomics - i.e. the qualitative and quantitative study of metabolites - is to identify chemical compounds in a given sample, oftentimes with the motivation of identifying/quantifying biomarkers for use in diagnosing, treating, and/or stratifying the risk of some disease. A central tool for compound identification in Metabolomics is mass spectrometry which produces a finite set of points - termed the mass spectrum - in the plane $\mathbb{R}^{2}$ for a given chemical compound with each point representing an ion fragment, one axis representing mass/charge, and the other axis representing intensity. The primary method used to identify a given unknown chemical compound based off of its mass spectrum is termed 'spectral library matching' and involves computing a measure of similarity between the mass spectrum of the unknown chemical compound and each mass spectrum of chemical compounds in a predetermined reference library. The unknown chemical compound is then identified as the chemical compound from the reference library whose mass spectrum is most similar to the mass spectrum of the unknown chemical compound. We present (insert package name), a command-line tool written in Python implementing spectral library matching with a host of spectrum preprocessing transformations and similarity measures available in addition to being able to analyze both of the two main types of mass spectrometry data: gas chromatography-mass spectrometry (GC-MS) and liquid chromatography-mass spectrometry (LC-MS) data.

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

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
