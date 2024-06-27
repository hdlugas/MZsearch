

from processing import *
from similarity_measures import *
import pandas as pd
import argparse
from pathlib import Path
import sys


# create new ArgumentParser object so we can extract command-line input
parser = argparse.ArgumentParser()

# Add optional command-line arguments
parser.add_argument('--query_data', metavar='\b', help='CSV file of query mass spectrum/spectra to be identified. Each row should correspond to a mass spectrum, the left-most column should contain an identifier, and each of the other columns should correspond to a single mass/charge ratio. Mandatory argument.')
parser.add_argument('--reference_data', metavar='\b', help='CSV file of the reference mass spectra. Each row should correspond to a mass spectrum, the left-most column should contain in identifier (i.e. the CAS registry number or the compound name), and the remaining column should correspond to a single mass/charge ratio. Default = LCMS GNPS Library.')
parser.add_argument('--similarity_measure', metavar='\b', help='Similarity measure: options are \'cosine\', \'shannon\', \'renyi\', and \'tsallis\'. Default = cosine.')
parser.add_argument('--spectrum_preprocessing_order', metavar='\b', help='The LCMS spectrum preprocessing transformations and the order in which they are to be applied. Note that these transformations are applied prior to computing similarity scores. Format must be a string with 2-6 characters chosen from C, F, M, N, L, W representing centroiding, filtering based on mass/charge and intensity values, matching, noise removal, low-entropy trannsformation, and weight-factor-transformation, respectively. For example, if \'WCM\' is passed, then each spectrum will undergo a weight factor transformation, then centroiding, and then matching. Note that if an argument is passed, then \'M\' must be contained in the argument, since matching is a required preprocessing step in spectral library matching of LCMS data. Default: FCNMWL')
parser.add_argument('--mz_min', metavar='\b', help='Remove all peaks with mass/charge less than mz_min in each spectrum. Default = 0')
parser.add_argument('--mz_max', metavar='\b', help='Remove all peaks with mass/charge greater than mz_max in each spectrum. Default = 999999999999')
parser.add_argument('--int_min', metavar='\b', help='Remove all peaks with intensity less than int_min in each spectrum. Default = 0')
parser.add_argument('--int_max', metavar='\b', help='Remove all peaks with intensity greater than int_max in each spectrum. Default = 999999999999')
parser.add_argument('--window_size_centroiding', metavar='\b', help='Window size parameter used in centroiding a given spectrum. Default = 0.5')
parser.add_argument('--window_size_matching', metavar='\b', help='Window size parameter used in matching a query spectrum and a reference library spectrum. Default = 0.5')
parser.add_argument('--noise_threshold', metavar='\b', help='Ion fragments (i.e. points in a given mass spectrum) with intensity less than max(intensities)*noise_threshold are removed. Default = 0')
parser.add_argument('--wf_mz', metavar='\b', help='Mass/charge weight factor parameter. Default = 0.')
parser.add_argument('--wf_intensity', metavar='\b', help='Intensity weight factor parameter. Default = 1.')
parser.add_argument('--LET_threshold', metavar='\b', help='Low-entropy transformation threshold parameter. Spectra with Shannon entropy less than LET_threshold are transformed according to intensitiesNew=intensitiesOriginal^{(1+S)/(1+LET_threshold)}. Default = 0.')
parser.add_argument('--entropy_dimension', metavar='\b', help='Entropy dimension parameter. Must have positive value other than 1. When the entropy dimension is 1, then Renyi and Tsallis entropy are equivalent to Shannon entropy. Therefore, this parameter only applies to the renyi and tsallis similarity measures. This parameter will be ignored if similarity measure cosine or shannon is chosen. Default = 1.1.')
parser.add_argument('--normalization_method', metavar='\b', help='Method used to normalize the intensities of each spectrum so that the intensities sum to 1. Since the objects entropy quantifies the uncertainy of must be probability distributions, the intensities of a given spectrum must sum to 1 prior to computing the entropy of the given spectrum intensities. Options: \'standard\' and \'softmax\'. Default = standard.')
parser.add_argument('--n_top_matches_to_save', metavar='\b', help='The number of top matches to report. For example, if n_top_matches_to_save=5, then for each query spectrum, the five reference spectra with the largest similarity with the given query spectrum will be reported. Default = 1.')
parser.add_argument('--print_id_results', metavar='\b', help='Flag that prints identification results if True. Default: False')
parser.add_argument('--output_identification', metavar='\b', help='Output CSV file containing the most-similar reference spectra for each query spectrum along with the corresponding similarity scores. Default is to save identification output in current working directory (i.e. same directory this script is contained in) with filename \'output_lcms_identification.csv\'.')
parser.add_argument('--output_similarity_scores', metavar='\b', help='Output CSV file containing similarity scores between all query spectrum/spectra and all reference spectra. Each row corresponds to a query spectrum, the left-most column contains the query spectrum/spectra identifier, and the remaining column contain the similarity scores with respect to all reference library spectra. If no argument passed, then this CSV file is written to the current working directory with filename \'output_lcms_all_similarity_scores\'.csv.')

# parse the user-input arguments
args = parser.parse_args()


print('\nPerforming spectral library matching on LCMS data')

# import the query library
if args.query_data is not None:
    df_query = pd.read_csv(args.query_data)
else:
    df_query = pd.read_csv(f'{Path.cwd()}/../data/lcms_query_library.csv')
    print('No argument passed to query_data; using default LCMS library')


# load the reference library
if args.reference_data is not None:
    df_reference = pd.read_csv(args.reference_data)
else:
    print('No argument passed to reference_data; using default LCMS reference library (from GNPS)\n')
    df_reference = pd.read_csv(f'{Path.cwd()}/../data/lcms_reference_library.csv')


# get the spectrum preprocessing order
preprocessing_error_message = 'Error: \'M\' must be a character in spectrum_preprocessing_order.'
if args.spectrum_preprocessing_order is not None:
    spectrum_preprocessing_order = list(args.spectrum_preprocessing_order)
else:
    spectrum_preprocessing_order = ['F', 'C', 'N', 'M', 'W', 'L']

if 'M' not in spectrum_preprocessing_order:
    print(f'\n{preprocessing_error_message}\n')
    sys.exit()


# load the weight factor parameters
if args.wf_mz is not None:
    wf_mz = float(args.wf_mz)
else:
    wf_mz = 0

if args.wf_intensity is not None:
    wf_intensity = float(args.wf_intensity)
else:
    wf_intensity = 1


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


# load the entropy dimension parameter (if applicable)
if args.similarity_measure == 'renyi' or args.similarity_measure == 'tsallis':
    if args.entropy_dimension is not None:
        q = float(args.entropy_dimension)
    else:
        q = 1.1


# load the centroiding window size parameter
if args.window_size_centroiding is not None:
    window_size_centroiding = float(args.window_size_centroiding)
else:
    window_size_centroiding = 0.5


# load the matching window size parameter
if args.window_size_matching is not None:
    window_size_matching = float(args.window_size_matching)
else:
    window_size_matching = 0.5


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


# get the number of most-similar reference library spectra to report for each query spectrum
if args.n_top_matches_to_save is not None:
    n_top_matches_to_save = int(args.n_top_matches_to_save)
else:
    n_top_matches_to_save = 1


# get the flag to determine whether or not to print identification results
if args.print_id_results is not None:
    print_id_results = str(args.print_id_results)
else:
    print_id_results = 'False'



# get unique query/reference library IDs; each query/reference ID corresponds to exactly one query/reference mass spectrum
unique_query_ids = df_query.iloc[:,0].unique()
unique_reference_ids = df_reference.iloc[:,0].unique()


# compute the similarity score between each query library spectrum/spectra and all reference library spectra
all_similarity_scores =  []
for query_idx in range(0,len(unique_query_ids)):
    q_idxs_tmp = np.where(df_query.iloc[:,0] == unique_query_ids[query_idx])[0]
    q_spec_tmp = np.asarray(pd.concat([df_query.iloc[q_idxs_tmp,1], df_query.iloc[q_idxs_tmp,2]], axis=1).reset_index(drop=True))
    q_spec = q_spec_tmp

    # compute the similarity score between the given query spectrum and all spectra in the reference library
    similarity_scores = []
    for ref_idx in range(0,len(unique_reference_ids)):
        #if ref_idx % 100 == 0:
        #    print(f'Query spectrum #{query_idx} has had its similarity with {ref_idx} reference library spectra computed')
        r_idxs_tmp = np.where(df_reference.iloc[:,0] == unique_reference_ids[ref_idx])[0]
        r_spec = np.asarray(pd.concat([df_reference.iloc[r_idxs_tmp,1], df_reference.iloc[r_idxs_tmp,2]], axis=1).reset_index(drop=True))

        is_matched = False
        for transformation in spectrum_preprocessing_order:
            if transformation == 'C' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                q_spec = centroid_spectrum(q_spec, window_size=window_size_centroiding) 
                r_spec = centroid_spectrum(r_spec, window_size=window_size_centroiding) 
            if transformation == 'M' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                m_spec = match_peaks_in_spectra(spec_a=q_spec, spec_b=r_spec, window_size=window_size_matching)
                q_spec = m_spec[:,0:2]
                r_spec = m_spec[:,[0,2]]
                is_matched = True
            if transformation == 'W' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                q_spec[:,1] = wf_transform(q_spec[:,0], q_spec[:,1], wf_mz, wf_intensity)
                r_spec[:,1] = wf_transform(r_spec[:,0], r_spec[:,1], wf_mz, wf_intensity)
            if transformation == 'L' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                q_spec[:,1] = LE_transform(q_spec[:,1], LET_threshold, normalization_method=normalization_method)
                r_spec[:,1] = LE_transform(r_spec[:,1], LET_threshold, normalization_method=normalization_method)
            if transformation == 'N' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                q_spec = remove_noise(q_spec, nr = noise_threshold)
                r_spec = remove_noise(r_spec, nr = noise_threshold)
            if transformation == 'F' and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
                q_spec = filter_spec_lcms(q_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)
                r_spec = filter_spec_lcms(r_spec, mz_min = mz_min, mz_max = mz_max, int_min = int_min, int_max = int_max, is_matched = is_matched)

        q_ints = q_spec[:,1]
        r_ints = r_spec[:,1]

        if np.sum(q_ints) != 0 and np.sum(r_ints) != 0 and q_spec.shape[0] > 1 and r_spec.shape[1] > 1:
            if similarity_measure == 'cosine':
                similarity_score = S_cos(q_ints, r_ints)
            else:
                q_ints = normalize(q_ints, method = normalization_method)
                r_ints = normalize(r_ints, method = normalization_method)

                if similarity_measure == 'shannon':
                    similarity_score = S_shannon(q_ints, r_ints)
                elif similarity_measure == 'renyi':
                    similarity_score = S_renyi(q_ints, r_ints, q)
                elif similarity_measure == 'tsallis':
                    similarity_score = S_tsallis(q_ints, r_ints, q)
        else:
            similarity_score = 0

        similarity_scores.append(similarity_score)
    all_similarity_scores.append(similarity_scores)


df_scores = pd.DataFrame(all_similarity_scores, columns = unique_reference_ids)
df_scores.index = unique_query_ids
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

        preds_tmp.append(';'.join(top_ref_specs_tmp.index.to_list()))
        if len(top_ref_specs_tmp.values) == 0:
            scores_tmp.append(0)
        else:
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
df_top_ref_specs.index = unique_query_ids
df_top_ref_specs.index.names = ['Query Spectrum ID']

if print_id_results == 'True':
    print(df_top_ref_specs.to_string())


# write spectral library matching results to disk
if args.output_identification is not None:
    df_top_ref_specs.to_csv(args.output_identification)
else:
    df_top_ref_specs.to_csv(f'{Path.cwd()}/output_lcms_identification.csv')


# write all similarity scores to disk
df_scores.columns = ['Reference Spectrum ID: ' + col for col in  list(map(str,df_scores.columns.tolist()))]
if args.output_similarity_scores is not None:
    df_scores.to_csv(args.output_similarity_scores)
else:
    df_scores.to_csv(f'{Path.cwd()}/output_lcms_all_similarity_scores.csv')




