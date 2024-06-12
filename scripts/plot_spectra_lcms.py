

from processing import *
from similarity_measures import *
import pandas as pd
import argparse
from pathlib import Path
import sys
import matplotlib.pyplot as plt


# create new ArgumentParser object so we can extract command-line input
parser = argparse.ArgumentParser()

# Add optional command-line arguments
parser.add_argument('--query_data', metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Default = LCMS GNPS Library.')
parser.add_argument('--query_spectrum_ID', metavar='\b', help='The identifier of the query spectrum to be plotted. Default: first query spectrum in query_data.')
parser.add_argument('--reference_spectrum_ID', metavar='\b', help='The identifier of the reference spectrum to be plotted. Default: first reference spectrum in reference_data.')
parser.add_argument('--similarity_measure', metavar='\b', help='Similarity measure: options are \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default = cosine.')
parser.add_argument('--spectrum_preprocessing_order', metavar='\b', help='The LCMS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-4 characters chosen from W, C, M, L representing weight-factor-transformation, cleaning (i.e. centroiding and noise removal), matching, and low-entropy transformation. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then cleaning, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of LCMS data. Default: CMWL')
parser.add_argument('--window_size', metavar='\b', help='Window size parameter used in (i) centroiding and (ii) matching a query spectrum and a reference library spectrum. Default = 0.5')
parser.add_argument('--noise_threshold', metavar='\b', help='Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default = 0')
parser.add_argument('--wf_mz', metavar='\b', help='Mass/charge weight factor parameter. Default = 0.')
parser.add_argument('--wf_intensity', metavar='\b', help='Intensity weight factor parameter. Default = 1.')
parser.add_argument('--LET_threshold', metavar='\b', help='Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default = 0')
parser.add_argument('--entropy_dimension', metavar='\b', help='Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default = 1.1.')
parser.add_argument('--normalization_method', metavar='\b', help='Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default = standard.')
parser.add_argument('--save_plots', metavar='\b', help='Output PDF file containing the plots of the query and reference spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./query_spec_{query_spectrum_ID}_reference_spec_{reference_spectrum_ID}_plot.pdf in the current working directory.')

# parse the user-input arguments
args = parser.parse_args()



# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    df_query = pd.read_csv(f'{Path.cwd()}/../data/lcms_query_library_tmp.csv')
    print('No argument passed to query_data; using default LCMS library')


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('Using default LCMS reference library (from GNPS)\n')
    df_reference = pd.read_csv(f'{Path.cwd()}/../data/lcms_reference_library.csv')


# import the identifier of the query spectrum to be plotted
if args.query_spectrum_ID is not None:
    query_spectrum_ID = str(args.query_spectrum_ID)
else:
    query_spectrum_ID = str(df_query.iloc[0,0])
    print('No argument passed to query_spectrum_ID; using the first spectrum in query_data')


# import the identifier of the reference spectrum to be plotted
if args.reference_spectrum_ID is not None:
    reference_spectrum_ID = str(args.reference_spectrum_ID)
else:
    reference_spectrum_ID = str(df_reference.iloc[0,0])
    print('No argument passed to reference_spectrum_ID; using the first spectrum in reference_data')


# import the identifier of the reference spectrum to be plotted
if args.save_plots is not None:
    path_output = args.save_plots
else:
    path_output = f'{Path.cwd()}/query_spec_{query_spectrum_ID}_reference_spec_{reference_spectrum_ID}_plot.pdf'




# get the spectrum preprocessing order
if args.spectrum_preprocessing_order is not None:
    spectrum_preprocessing_order = list(args.spectrum_preprocessing_order)
else:
    spectrum_preprocessing_order = ['C', 'M', 'W', 'L']


# load the weight factor parameters
if args.wf_mz is not None:
    wf_mz = float(args.wf_mz)
else:
    wf_mz = 0

if args.wf_intensity is not None:
    wf_intensity = float(args.wf_intensity)
else:
    wf_intensity = 1


# load the entropy dimension parameter (if applicable)
if args.similarity_measure == 'renyi' or args.similarity_measure == 'tsallis':
    if args.entropy_dimension is not None:
        q = float(args.entropy_dimension)
    else:
        q = 1.1

# load the window size parameter
if args.window_size is not None:
    window_size = float(args.window_size)
else:
    window_size = 0.5


# load the noise removal parameter
if args.noise_threshold is not None:
    noise_threshold = float(args.noise_threshold)
else:
    noise_threshold = 0



# load the low-entropy transformation threshold
if args.LET_threshold is not None: 
    LET_threshold = float(args.LET_threshold)
else:
    LET_threshold = 0


# set the normalization method
if args.normalization_method is not None:
    normalization_method = args.normalization_method
else:
    normalization_method = 'standard'


# specify the similarity measure to use
if args.similarity_measure is not None:
    similarity_measure = args.similarity_measure
else:
    similarity_measure = 'cosine'



# get unique query/reference library IDs; each query/reference ID corresponds to exactly one query/reference mass spectrum
unique_query_ids = df_query.iloc[:,0].unique()
unique_reference_ids = df_reference.iloc[:,0].unique()

unique_query_ids = np.asarray(list(map(str,unique_query_ids)))
unique_reference_ids = np.asarray(list(map(str,unique_reference_ids)))
df_query.iloc[:,0] = list(map(str,df_query.iloc[:,0].tolist()))
df_reference.iloc[:,0] = list(map(str,df_reference.iloc[:,0].tolist()))

query_idx = np.where(unique_query_ids == query_spectrum_ID)[0][0]
reference_idx = np.where(unique_reference_ids == reference_spectrum_ID)[0][0]

q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[reference_idx])[0]
q_spec = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
r_spec = np.asarray(pd.concat([df_reference.iloc[q_idxs_tmp,1], df_reference.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
#print(q_spec)
#print(r_spec)


fig, axes = plt.subplots(nrows=2, ncols=1)

plt.subplot(2,1,1)
plt.vlines(x=q_spec[:,0], ymin=[0]*q_spec.shape[0], ymax=q_spec[:,1], linewidth=4, color='blue', label=f'Query Spectrum ID: {query_spectrum_ID}')
plt.vlines(x=r_spec[:,0], ymin=[0]*r_spec.shape[0], ymax=r_spec[:,1], linewidth=3, color='red', label=f'Reference Spectrum ID: {reference_spectrum_ID}')
plt.legend(loc='upper right', fontsize=8)
plt.xlabel('Mass:Charge Ratio',fontsize=8)
plt.ylabel('Intensity', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title('Untransformed Query and Reference Spectra', fontsize=12)

for transformation in spectrum_preprocessing_order:
    if transformation == 'C':
        q_spec = clean_spectrum(q_spec, noise_removal=noise_threshold, window_size=window_size) 
        r_spec = clean_spectrum(r_spec, noise_removal=noise_threshold, window_size=window_size) 
    if transformation == 'M':
        m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=window_size)
        q_spec = m_spec[:,0:2]
        r_spec = m_spec[:,[0,2]]
    if transformation == 'W':
        q_spec[:,1] = np.power(q_spec[:,0], wf_mz) * np.power(q_spec[:,1], wf_intensity)
        r_spec[:,1] = np.power(r_spec[:,0], wf_mz) * np.power(r_spec[:,1], wf_intensity)
    if transformation == 'L':
        q_spec[:,1] = transform_int(q_spec[:,1], LET_threshold, normalization_method)
        r_spec[:,1] = transform_int(r_spec[:,1], LET_threshold, normalization_method)


if similarity_measure == 'cosine':
    similarity_score = S_cos(q_spec[:,1], r_spec[:,1])
else:
    q_spec[:,1] = normalize(q_spec[:,1], method = normalization_method)
    r_spec[:,1] = normalize(r_spec[:,1], method = normalization_method)
if similarity_measure == 'shannon':
    similarity_score = S_shannon(q_spec[:,1].to_numpy(), r_spec[:,1].to_numpy())
elif similarity_measure == 'renyi':
    similarity_score = S_renyi(q_spec[:,1].to_numpy(), r_spec[:,1].to_numpy(), q)
elif similarity_measure == 'tsallis':
    similarity_score = S_tsallis(q_spec[:,1].to_numpy(), r_spec[:,1].to_numpy(), q)


plt.subplot(2,1,2)
plt.vlines(x=q_spec[:,0], ymin=[0]*q_spec.shape[0], ymax=q_spec[:,1], linewidth=4, color='blue', label=f'Query Spectrum ID: {query_spectrum_ID}')
plt.vlines(x=r_spec[:,0], ymin=[0]*r_spec.shape[0], ymax=r_spec[:,1], linewidth=3, color='red', label=f'Reference Spectrum ID: {reference_spectrum_ID}')
plt.legend(loc='upper right', fontsize=8)
plt.xlabel('Mass:Charge Ratio', fontsize=8)
plt.ylabel('Intensity', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(f'Transformed Query and Reference Spectra\n Similarity Score: {round(similarity_score,4)}', fontsize=12)

plt.subplots_adjust(hspace=0.7)
plt.savefig(path_output, format='pdf')



