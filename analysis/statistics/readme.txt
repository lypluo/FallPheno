Code Running order:
1_Data_in_growingseason_filter_outliers.R： calculate the GPP/NEP mean and cumulative values (as well as anomaly) during different periods(sos25-eos90,sos25-eos50,sos25-eos25) in growing season.
2_further_evaluate_XEP_vs_eos.R：Apart from analyzing the data in periods of :[sos25,solstice],[sos25,eos90],[sos25,eos50],[sos25,eos25], I also add the analysis for following period on Feb,2022:
[solstice,eos90],[eos90,eos50],[eos50,eos25]
3_slopes_of_XEP_vs_eos.R: conduct the linear regression between XEP and eos, and tidy its slopes.

##
add_0_additional_compare_Lu2022.R：comparing the phenophases extracted from Lu et al., 2022 and this analysis. 