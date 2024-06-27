

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
parser.add_argument('--spectrum_preprocessing_order', metavar='\b', help='The GCMS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-4 characters chosen from F, N, L, W representing filtering based on mass/charge and intensity values, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'LW\' is passed, then each spectrum will undergo a low-entropy transformation and then a weight factor transformation. Default: FNLW')
parser.add_argument('--mz_min', metavar='\b', help='Remove all peaks with mass/charge less than mz_min in each spectrum. Default = 0')
parser.add_argument('--mz_max', metavar='\b', help='Remove all peaks with mass/charge greater than mz_max in each spectrum. Default = 999999999999')
parser.add_argument('--int_min', metavar='\b', help='Remove all peaks with intensity less than int_min in each spectrum. Default = 0')
parser.add_argument('--int_max', metavar='\b', help='Remove all peaks with intensity greater than int_max in each spectrum. Default = 999999999999')
parser.add_argument('--noise_threshold', metavar='\b', help='Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default = 0')
parser.add_argument('--wf_mz', metavar='\b', help='Mass/charge weight factor parameter. Default = 0.')
parser.add_argument('--wf_intensity', metavar='\b', help='Intensity weight factor parameter. Default = 1.')
parser.add_argument('--LET_threshold', metavar='\b', help='Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default = 0')
parser.add_argument('--normalization_method', metavar='\b', help='Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default = standard.')
parser.add_argument('--entropy_dimension', metavar='\b', help='Entropy dimension parameter. Note that this only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similairty measure cosine or shannon is chosen. Default = 1.1.')
parser.add_argument('--save_plots', metavar='\b', help='Output PDF file containing the plots of the query and reference spectra before and after preprocessing transformations. If no argument is passed, then the plots will be saved to the PDF ./query_spec_{query_spectrum_ID}_reference_spec_{reference_spectrum_ID}_plot.pdf in the current working directory.')

# parse the user-input arguments
args = parser.parse_args()



# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    df_query = pd.read_csv(f'{Path.cwd()}/../data_all/gcms_query_library_tmp.csv')
    print('No argument passed to query_data; using default GCMS NIST WebBook library')


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('No argument passed to reference_data; using default GCMS reference library (i.e. NIST WebBook)')
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


# get the spectrum preprocessing order
if args.spectrum_preprocessing_order is not None:
    spectrum_preprocessing_order = list(args.spectrum_preprocessing_order)
else:
    spectrum_preprocessing_order = ['F', 'N', 'L', 'W']


# load the filtering parameters
if args.mz_min is not None:
    mz_min = float(args.mz_min)
else: 
    mz_min = 0

if args.mz_max is not None:
    mz_max = float(args.mz_max)
else: 
    mz_max = 999999999999

if args.int_min is not None:
    int_min = float(args.int_min)
else: 
    int_min = 0

if args.int_max is not None:
    int_max = float(args.int_max)
else: 
    int_max = 999999999999


# load the noise removal parameter
if args.noise_threshold is not None:
    noise_threshold = float(args.noise_threshold)
else:
    noise_threshold = 0


# load the weight factor parameters
if args.wf_mz is not None:
    wf_mz = float(args.wf_mz)
else:
    wf_mz = 0

if args.wf_intensity is not None:
    wf_intensity = float(args.wf_intensity)
else:
    wf_intensity = 1


# load the low-entropy transformation threshold
if args.LET_threshold is not None: 
    LET_threshold = float(args.LET_threshold)
else:
    LET_threshold = 0



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


df_query = df_query.astype(object)
df_reference = df_reference.astype(object)
df_query.iloc[:,0] = df_query.iloc[:,0].astype(str)
df_reference.iloc[:,0] = df_reference.iloc[:,0].astype(str)
#df_query.iloc[:,0] = list(map(str,df_query.iloc[:,0].tolist()))
#df_reference.iloc[:,0] = list(map(str,df_reference.iloc[:,0].tolist()))
q_idx = np.where(df_query.iloc[:,0] == query_spectrum_ID)[0][0]
r_idx = np.where(df_reference.iloc[:,0] == reference_spectrum_ID)[0][0]
q_ints = df_query.iloc[q_idx,1:df_query.shape[1]].to_numpy()
r_ints = df_reference.iloc[r_idx,1:df_reference.shape[1]].to_numpy()

max_mz_tmp = max([np.max(np.nonzero(q_ints)), np.max(np.nonzero(r_ints))])
mzs = list(map(int,np.linspace(1,max_mz_tmp,max_mz_tmp).tolist()))
q_ints = q_ints[mzs]
r_ints = r_ints[mzs]

q_spec = np.transpose(np.array([mzs, q_ints]))
r_spec = np.transpose(np.array([mzs, r_ints]))

fig, axes = plt.subplots(nrows=2, ncols=1)
#fig.tight_layout()

plt.subplot(2,1,1)
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=q_spec[:,1]/np.max(q_spec[:,1]), linewidth=2, color='blue', label=f'Query Spectrum ID: {query_spectrum_ID}')
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=r_spec[:,1]/np.max(r_spec[:,1]), linewidth=2, color='red', label=f'Reference Spectrum ID: {reference_spectrum_ID}')
plt.legend(loc='upper right', fontsize=8)
plt.xlabel('Mass:Charge Ratio',fontsize=8)
plt.ylabel('Intensity', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title('Untransformed Query and Reference Spectra', fontsize=12)


for transformation in spectrum_preprocessing_order:
    if transformation == 'W':
        q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
        r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
    if transformation == 'L':
        q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method)
        r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method)
    if transformation == 'N':
        q_spec = remove_noise(q_spec, nr = noise_threshold)
        r_spec = remove_noise(r_spec, nr = noise_threshold)
    if transformation == 'F':
        q_spec = filter_spec_gcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)
        r_spec = filter_spec_gcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max)


if q_spec.shape[0] > 1:
    if similarity_measure == 'cosine':
        similarity_score = S_cos(q_spec[:,1], r_spec[:,1])
    else:
        q_spec[:,1] = normalize(q_spec[:,1], method = normalization_method)
        r_spec[:,1] = normalize(r_spec[:,1], method = normalization_method)

        if similarity_measure == 'shannon':
            similarity_score = S_shannon(q_spec[:,1].astype('float'), r_spec[:,1].astype('float'))
        elif similarity_measure == 'renyi':
            similarity_score = S_renyi(q_spec[:,1], r_spec[:,1], q)
        elif similarity_measure == 'tsallis':
            similarity_score = S_tsallis(q_spec[:,1], r_spec[:,1], q)
else:
    similarity_score = 0


plt.subplot(2,1,2)
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=q_spec[:,1]/np.max(q_spec[:,1]), linewidth=2, color='blue', label=f'Query Spectrum ID: {query_spectrum_ID}')
plt.vlines(x=mzs, ymin=[0]*len(mzs), ymax=r_spec[:,1]/np.max(r_spec[:,1]), linewidth=2, color='red', label=f'Reference Spectrum ID: {reference_spectrum_ID}')
plt.legend(loc='upper right', fontsize=8)
plt.xlabel('Mass:Charge Ratio', fontsize=8)
plt.ylabel('Intensity', fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title(f'Transformed Query and Reference Spectra\n Similarity Score: {round(similarity_score,4)}', fontsize=12)
#plt.show()

plt.subplots_adjust(hspace=0.7)
plt.savefig(path_output, format='pdf')




