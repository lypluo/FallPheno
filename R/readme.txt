code running steps:
preparing the data and phenophase extraction code: 
tidy_data_fromJiangong.R
  -->tidy_MODIS_VIs_fromWalther.R--extract the VIs from Walther et al., 2021
     -->SSplineFilter_gapfill.R-->filter the outliers of times series 
       -->pheno_extraction_fun.R--prepare the phenophase extraction functions.
         -->Extract_Phenophase.Rmd--extract the phenology from VI(firstly for EVI)
  

  