

df_idrep = read.csv('/wsu/home/fy/fy73/fy7392/mass_spec/data/idrep.csv')
df_rep = read.csv('/wsu/home/fy/fy73/fy7392/mass_spec/data/rep.csv')
df = cbind(df_idrep, df_rep)
write.csv(df, '/wsu/home/groups/tcgabc/2DMS/JOSS/data/gcms_query_library.csv', row.names=F)


df_idnist = read.csv('/wsu/home/fy/fy73/fy7392/mass_spec/data/idnist.csv')
df_nist = read.csv('/wsu/home/fy/fy73/fy7392/mass_spec/data/nist.csv')
df = cbind(df_idnist, df_nist)
write.csv(df, '/wsu/home/groups/tcgabc/2DMS/JOSS/data/gcms_reference_library.csv', row.names=F)


