

from processing import *
from similarity_measures import *
import pandas as pd
import argparse
from pathlib import Path
import sys
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser()

# Add optional command-line arguments
parser.add_argument('--query_data', metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Default = GCMS NIST WebBook Library.')
parser.add_argument('--query_spectrum_ID', metavar='\b', help='The identifier of the query spectrum to be plotted. Default: first query spectrum in query_data.')
parser.add_argument('--reference_spectrum_ID', metavar='\b', help='The identifier of the reference spectrum to be plotted. Default: first reference spectrum in reference_data.')
parser.add_argument('--similarity_measure', metavar='\b', help='Similarity measure: options are \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default = cosine.')
parser.add_argument('--wf_mz', metavar='\b', help='Mass/charge weight factor parameter. Default = 0.')
parser.add_argument('--wf_intensity', metavar='\b', help='Intensity weight factor parameter. Default = 1.')
parser.add_argument('--normalization_method', metavar='\b', help='Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default = standard.')
parser.add_argument('--entropy_dimension', metavar='\b', help='Entropy dimension parameter. Note that this only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similairty measure cosine or shannon is chosen. Default = 1.1.')
parser.add_argument('--save_plots', metavar='\b', help='Output PDF file containing the plots of the query and reference spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./query_spec_{query_spectrum_ID}_reference_spec_{reference_spectrum_ID}_plot.pdf in the current working directory.')

# parse the user-input arguments
args = parser.parse_args()



# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    df_query = pd.read_csv(f'{Path.cwd()}/../data/gcms_query_library_tmp.csv')
    print('No argument passed to query_data; using default GCMS NIST WebBook library')


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('Using default GCMS reference library (i.e. NIST WebBook)\n')
    df_reference = pd.read_csv(f'{Path.cwd()}/../data/gcms_reference_library.csv')


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


df_query.iloc[:,0] = list(map(str,df_query.iloc[:,0].tolist()))
df_reference.iloc[:,0] = list(map(str,df_reference.iloc[:,0].tolist()))
query_idx = np.where(df_query.iloc[:,0] == query_spectrum_ID)[0][0]
reference_idx = np.where(df_reference.iloc[:,0] == reference_spectrum_ID)[0][0]
query_spec = df_query.iloc[query_idx,1:df_query.shape[1]].to_numpy()
reference_spec = df_reference.iloc[reference_idx,1:df_reference.shape[1]].to_numpy()

max_mz = max([np.max(np.nonzero(query_spec)), np.max(np.nonzero(reference_spec))])
mzs = list(map(int,np.linspace(1,max_mz,max_mz).tolist()))
query_spec = query_spec[mzs]
reference_spec = reference_spec[mzs]

fig, axes = plt.subplots(nrows=2, ncols=1)
#fig.tight_layout()

plt.subplot(2,1,1)
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=query_spec, linewidth=2, color='blue', label=f'Query Spectrum ID: {query_spectrum_ID}')
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=reference_spec, linewidth=2, color='red', label=f'Reference Spectrum ID: {reference_spectrum_ID}')
plt.legend(loc='upper right', fontsize=8)
plt.xlabel('Mass:Charge Ratio',fontsize=8)
plt.ylabel('Intensity', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title('Untransformed Query and Reference Spectra', fontsize=12)

query_spec = np.power(mzs, wf_mz) * np.power(query_spec, wf_intensity)
reference_spec = np.power(mzs, wf_mz) * np.power(reference_spec, wf_intensity)


if similarity_measure == 'cosine':
    similarity_score = S_cos(query_spec, reference_spec)
else:
    query_spec = normalize(query_spec, method = normalization_method)
    reference_spec = normalize(reference_spec, method = normalization_method)

if similarity_measure == 'shannon':
    similarity_score = S_shannon(query_spec.to_numpy(), reference_spec.to_numpy())
elif similarity_measure == 'renyi':
    similarity_score = S_renyi(query_spec.to_numpy(), reference_spec.to_numpy(), q)
elif similarity_measure == 'tsallis':
    similarity_score = S_tsallis(query_spec.to_numpy(), reference_spec.to_numpy(), q)


plt.subplot(2,1,2)
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=query_spec, linewidth=2, color='blue', label=f'Query Spectrum ID: {query_spectrum_ID}')
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=reference_spec, linewidth=2, color='red', label=f'Reference Spectrum ID: {reference_spectrum_ID}')
plt.legend(loc='upper right', fontsize=8)
plt.xlabel('Mass:Charge Ratio', fontsize=8)
plt.ylabel('Intensity', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(f'Transformed Query and Reference Spectra\n Similarity Score: {round(similarity_score,4)}', fontsize=12)
#plt.show()

plt.subplots_adjust(hspace=0.7)
plt.savefig(path_output, format='pdf')




