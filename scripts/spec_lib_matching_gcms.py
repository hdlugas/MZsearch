

from processing import *
from similarity_measures import *
import pandas as pd
import argparse
from pathlib import Path
import sys


parser = argparse.ArgumentParser()

# Add optional command-line arguments
parser.add_argument('--query_data', metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Default = GCMS NIST WebBook Library.')
parser.add_argument('--similarity_measure', metavar='\b', help='Similarity measure: options are \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default = cosine.')
parser.add_argument('--spectrum_preprocessing_order', metavar='\b', help='The GCMS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-4 characters chosen from F, N, L, W representing filtering based on mass/charge and intensity values, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'LW\' is passed, then each spectrum will undergo a low-entropy transformation and then a weight factor transformation. Default: FNLW')
parser.add_argument('--mz_min', metavar='\b', help='Remove all peaks with mass/charge less than mz_min in each spectrum. Default = 0')
parser.add_argument('--mz_max', metavar='\b', help='Remove all peaks with mass/charge greater than mz_max in each spectrum. Default = 999999999999')
parser.add_argument('--int_min', metavar='\b', help='Remove all peaks with intensity less than int_min in each spectrum. Default = 0')
parser.add_argument('--int_max', metavar='\b', help='Remove all peaks with intensity greater than int_max in each spectrum. Default = 999999999999')
parser.add_argument('--wf_mz', metavar='\b', help='Mass/charge weight factor parameter. Default = 0.')
parser.add_argument('--wf_intensity', metavar='\b', help='Intensity weight factor parameter. Default = 1.')
parser.add_argument('--entropy_dimension', metavar='\b', help='Entropy dimension parameter. Note that this only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similairty measure cosine or shannon is chosen. Default = 1.1.')
parser.add_argument('--normalization_method', metavar='\b', help='Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default = standard.')
parser.add_argument('--n_top_matches_to_save', metavar='\b', help='The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default = 1.')
parser.add_argument('--output_identification', metavar='\b', help='Output CSV file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory (i.e. same directory this script is contained in) with filename \'output_gcms_identification.csv\'.')
parser.add_argument('--output_similarity_scores', metavar='\b', help='Output CSV file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this CSV file is written to the current working directory with filename \'output_gcms_all_similarity_scores\'.csv.')
parser.add_argument('-v', '--verbose', action='store_true')


# parse the user-input arguments
args = parser.parse_args()



# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    df_query = pd.read_csv(f'{Path.cwd()}/../data_all/gcms_query_library_tmp.csv')
    print('No argument passed to query_data; using default GCMS NIST WebBook library')



# begin spectral library matching
print('Performing spectral library matching on GCMS data\n')


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('Using default GCMS reference library (i.e. NIST WebBook)\n')
    df_reference = pd.read_csv(f'{Path.cwd()}/../data/gcms_reference_library.csv')


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


# get the number of most-similar reference library spectra to report for each query spectrum
if args.n_top_matches_to_save is not None:
    n_top_matches_to_save = int(args.n_top_matches_to_save)
else:
    n_top_matches_to_save = 1



# compute the similarity score between each query library spectrum/spectra and all reference library spectra
mzs = np.linspace(1,df_query.shape[1]-1,df_query.shape[1]-1)
all_similarity_scores =  []
for query_idx in range(0,df_query.shape[0]):
    query_spec = df_query.iloc[query_idx,1:df_query.shape[1]]
    query_spec = wf_transform(mzs, query_spec, wf_mz, wf_intensity)

    similarity_scores = []
    for ref_idx in range(0, df_reference.shape[0]):
        if ref_idx % 1000 == 0:
            print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
        ref_spec = df_reference.iloc[ref_idx,1:df_reference.shape[1]]
        ref_spec = wf_transform(mzs, ref_spec, wf_mz, wf_intensity)

        if similarity_measure == 'cosine':
            similarity_score = S_cos(query_spec, ref_spec)
        else:
            query_spec = normalize(query_spec, method = normalization_method)
            ref_spec = normalize(ref_spec, method = normalization_method)

            if similarity_measure == 'shannon':
                similarity_score = S_shannon(query_spec.to_numpy(), ref_spec.to_numpy())
            elif similarity_measure == 'renyi':
                similarity_score = S_renyi(query_spec.to_numpy(), ref_spec.to_numpy(), q)
            elif similarity_measure == 'tsallis':
                similarity_score = S_tsallis(query_spec.to_numpy(), ref_spec.to_numpy(), q)

        similarity_scores.append(similarity_score)
    all_similarity_scores.append(similarity_scores)


df_scores = pd.DataFrame(all_similarity_scores, columns = df_reference.iloc[:,0])
df_scores.index = df_query.iloc[:,0]
df_scores.index.names = ['Query Spectrum ID']


preds = []
scores = []
for i in range(0, df_scores.shape[0]):
    df_scores_tmp = df_scores
    preds_tmp = []
    scores_tmp = []
    for j in range(0, n_top_matches_to_save):
        top_ref_specs_tmp = df_scores_tmp.iloc[i,np.where(df_scores_tmp.iloc[i,:] == np.max(df_scores_tmp.iloc[i,:]))[0]]
        cols_to_keep = np.where(df_scores_tmp.iloc[i,:] != np.max(df_scores_tmp.iloc[i,:]))[0]
        df_scores_tmp = df_scores_tmp.iloc[:,cols_to_keep]

        preds_tmp.append(';'.join(map(str,top_ref_specs_tmp.index.to_list())))
        scores_tmp.append(top_ref_specs_tmp.values[0])
    preds.append(preds_tmp)
    scores.append(scores_tmp)

preds = np.array(preds)
scores = np.array(scores)
out = np.c_[preds,scores]

cnames_preds = []
cnames_scores = []
for i in range(0,n_top_matches_to_save):
    cnames_preds.append(f'RANK.{i+1}.PRED')
    cnames_scores.append(f'RANK.{i+1}.SIMILARITY.SCORE')

df_top_ref_specs = pd.DataFrame(out, columns = [*cnames_preds, *cnames_scores])
df_top_ref_specs.index = df_query.iloc[:,0]
df_top_ref_specs.index.names = ['Query Spectrum ID']



# write spectral library matching results to disk
if args.output_identification is not None:
    df_top_ref_specs.to_csv(args.output_identification)
else:
    df_top_ref_specs.to_csv(f'{Path.cwd()}/output_gcms_identification.csv')


# write all similarity scores to disk
df_scores.columns = ['Reference Spectrum ID ' + col for col in  list(map(str,df_scores.columns.tolist()))]
if args.output_similarity_scores is not None:
    df_scores.to_csv(args.output_similarity_scores)
else:
    df_scores.to_csv(f'{Path.cwd()}/output_gcms_all_similarity_scores.csv')



