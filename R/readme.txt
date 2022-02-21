code running steps:
preparing the data and phenophase extraction code: 
tidy_data_fromJiangong.R
  -->tidy_MODIS_VIs_fromWalther.R--extract the VIs from Walther et al., 2021
     -->SSplineFilter_gapfill.R-->filter the outliers of times series 
       -->pheno_extraction_fun.R--prepare the phenophase extraction functions：pheno_extraction_fun1 is the original extraction method; pheno_extraction_fun2 are the updated extraction method after reading Zohner et al., 2022.
         -->Extract_Phenophase.Rmd--extract the phenology from VI(firstly for EVI):Extract_Phenophase_proc1.Rmd is the original extraction code;Extract_Phenophase_proc2.Rmd is the updated extraction code from 2022, Feb.
  

  