
suppressMessages(library(MSnbase))
suppressMessages(library(optparse))


option_list = list(make_option(c("-i", "--input_path"), type="character", default=NULL, help="path to input MGF file", metavar="character"),
                   make_option(c("-o", "--output_path"), type="character", default=NULL, help="path to output CSV file", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (length(opt) != 3){
  stop('ERROR: User must pass both input_path to MGF file and output_path for CSV file.')
}

mgf_data = readMgfData(opt$input_path)
get_spectra_df = function(mgf_data){
  ids = c()
  mzs = c()
  ints = c()
  for (i in 1:length(mz(mgf_data))){
    print(i)
    mzs_tmp = mz(mgf_data)[i]
    ints_tmp = intensity(mgf_data)[i]
    mzs_tmp = unlist(strsplit(as.character(mzs_tmp[[1]]), ' '))
    ints_tmp = unlist(strsplit(as.character(ints_tmp[[1]]), ' '))
    ids = c(ids, rep(i,length(mzs_tmp)))
    mzs = c(mzs, mzs_tmp)
    ints = c(ints, ints_tmp)
  }
  df_out = data.frame(SPECTRUM.ID=ids, MZ=mzs, INT=ints)
}

spectra_df = get_spectra_df(mgf_data)
write.csv(spectra_df, opt$output_path, row.names=F)


