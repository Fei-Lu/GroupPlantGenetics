# Meta-Tissue pipeline

# Prerequisites
1.nominal files :
header must includes
```bash
gene_id	variant_id	tss_distance	pval_nominal	slope	slope_se
```

2.significant pairs files :
header must includes
```bash
variant_id	gene_id	
```

# output files

output/Summary/meta_summary.metasoft.txt.gz
```bash
RSID	#STUDY	PVALUE_FE	BETA_FE	STD_FE	PVALUE_RE	BETA_RE	STD_RE	PVALUE_RE2	STAT1_RE2	STAT2_RE2	PVALUE_BE	I_SQUARE	Q	PVALUE_Q	TAU_SQUARE	pval_Adipose_Subcutaneous_Analysis	pval_Adipose_Visceral_Omentum_Analysis	pval_Adrenal_Gland_Analysis	pval_Artery_Aorta_Analysis	pval_Artery_Coronary_Analysis	pval_Artery_Tibial_Analysis	pval_Brain_Anterior_cingulate_cortex_BA24_Analysis	pval_Brain_Caudate_basal_ganglia_Analysis	pval_Brain_Cerebellar_Hemisphere_Analysis	pval_Brain_Cerebellum_Analysis	pval_Brain_Cortex_Analysis	pval_Brain_Frontal_Cortex_BA9_Analysis	pval_Brain_Hippocampus_Analysis	pval_Brain_Hypothalamus_Analysis	pval_Brain_Nucleus_accumbens_basal_ganglia_Analysis	pval_Brain_Putamen_basal_ganglia_Analysis	pval_Breast_Mammary_Tissue_Analysis	pval_Cells_EBV-transformed_lymphocytes_Analysis	pval_Cells_Transformed_fibroblasts_Analysis	pval_Colon_Sigmoid_Analysis	pval_Colon_Transverse_Analysis	pval_Esophagus_Gastroesophageal_Junction_Analysis	pval_Esophagus_Mucosa_Analysis	pval_Esophagus_Muscularis_Analysis	pval_Heart_Atrial_Appendage_Analysis	pval_Heart_Left_Ventricle_Analysis	pval_Liver_Analysis	pval_Lung_Analysis	pval_Muscle_Skeletal_Analysis	pval_Nerve_Tibial_Analysis	pval_Ovary_Analysis	pval_Pancreas_Analysis	pval_Pituitary_Analysis	pval_Prostate_Analysis	pval_Skin_Not_Sun_Exposed_Suprapubic_Analysis	pval_Skin_Sun_Exposed_Lower_leg_Analysis	pval_Small_Intestine_Terminal_Ileum_Analysis	pval_Spleen_Analysis	pval_Stomach_Analysis	pval_Testis_Analysis	pval_Thyroid_Analysis	pval_Uterus_Analysis	pval_Vagina_Analysis	pval_Whole_Blood_Analysis	mval_Adipose_Subcutaneous_Analysis	mval_Adipose_Visceral_Omentum_Analysis	mval_Adrenal_Gland_Analysis	mval_Artery_Aorta_Analysis	mval_Artery_Coronary_Analysis	mval_Artery_Tibial_Analysis	mval_Brain_Anterior_cingulate_cortex_BA24_Analysis	mval_Brain_Caudate_basal_ganglia_Analysis	mval_Brain_Cerebellar_Hemisphere_Analysis	mval_Brain_Cerebellum_Analysis	mval_Brain_Cortex_Analysis	mval_Brain_Frontal_Cortex_BA9_Analysis	mval_Brain_Hippocampus_Analysis	mval_Brain_Hypothalamus_Analysis	mval_Brain_Nucleus_accumbens_basal_ganglia_Analysis	mval_Brain_Putamen_basal_ganglia_Analysis	mval_Breast_Mammary_Tissue_Analysis	mval_Cells_EBV-transformed_lymphocytes_Analysis	mval_Cells_Transformed_fibroblasts_Analysis	mval_Colon_Sigmoid_Analysis	mval_Colon_Transverse_Analysis	mval_Esophagus_Gastroesophageal_Junction_Analysis	mval_Esophagus_Mucosa_Analysis	mval_Esophagus_Muscularis_Analysis	mval_Heart_Atrial_Appendage_Analysis	mval_Heart_Left_Ventricle_Analysis	mval_Liver_Analysis	mval_Lung_Analysis	mval_Muscle_Skeletal_Analysis	mval_Nerve_Tibial_Analysis	mval_Ovary_Analysis	mval_Pancreas_Analysis	mval_Pituitary_Analysis	mval_Prostate_Analysis	mval_Skin_Not_Sun_Exposed_Suprapubic_Analysis	mval_Skin_Sun_Exposed_Lower_leg_Analysis	mval_Small_Intestine_Terminal_Ileum_Analysis	mval_Spleen_Analysis	mval_Stomach_Analysis	mval_Testis_Analysis	mval_Thyroid_Analysis	mval_Uterus_Analysis	mval_Vagina_Analysis	mval_Whole_Blood_Analysis
1_55164_C_A_b37,ENSG00000224956.5	20	5.26952e-11	0.258737	0.0394229	4.02156e-06	0.265011	0.0574827	6.93307e-12	43.0745	5.86803	NA	51.1915	38.9276	0.00451328	0.0328204	0.0453393	0.525147	NA	0.274242	NA	0.134894	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0.00822152	0.0277308	0.783856	NA	NA	0.0142041	0.238509	0.000420907	0.543647	NA	NA	0.00425359	0.349467	0.420429	NA	NA	NA	NA	5.2195e-05	0.0412067	NA	NA	0.0544247	0.501342	0.0709952	NA	NA	0.151798	0.937	0.643	NA	0.812	NA	0.852	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0.989	0.934	0.539	NA	NA	0.978	0.774	1	0.346	NA	NA	0.988	0.614	0.529	NA	NA	NA	NA	1	0.951	NA	NA	0.923	0.179	0.894	NA	NA	0.025
1_55299_C_T_b37,ENSG00000224956.5	44	8.13752e-18	-0.20736	0.0241182	5.16636e-10	-0.197363	0.0317612	4.00066e-18	73.9193	3.89023	NA	36.2138	67.4127	0.0101001	0.0146772	0.00396313	0.0161918	0.9382	0.00113651	0.00294684	0.066773	0.331159	0.110476	0.961633	0.595057	0.659585	0.134764	0.00514336	0.288588	0.703668	0.403706	0.983961	0.235242	0.13854	0.0136977	0.0182148	0.102395	0.594881	0.0891114	0.309354	0.946334	0.641565	0.0174463	0.0255503	2.02212e-07	0.166701	0.262765	0.522314	0.145819	0.102365	0.244988	0.298329	0.495241	0.29353	0.048439	0.007991	0.758177	0.631149	0.893613	0.991	0.973	0.632	0.989	0.987	0.931	0.556	0.94	0.631	0.555	0.667	0.918	0.126	0.859	0.667	0.64	0.591	0.459	0.895	0.971	0.977	0.949	0.319	0.933	0.86	0.506	0.712	0.97	0.982	1	0.88	0.801	0.82	0.897	0.932	0.751	0.845	0.7780.82	0.956	0.991	0.751	0.75	0.337
1_55326_T_C_b37,ENSG00000187583.6	8	0.0155583	-0.15894	0.0657019	0.114243	-0.218739	0.138494	0.000299749	5.85212	8.17671	NA	74.1142	27.0418	0.000327555	0.10732	0.00587521	NA	NA	NA	NA	0.0564325	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0.730297	NA	NA	0.150424	NA	NA	NA	NA	NA	0.697602	0.667234	NA	NA	NA	NA	0.40102	NA	NA	2.6e-05	0.92	NA	NA	NA	NA	0.768	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	0.232	NA	NA	0.073	NA	NA	NA	NA	NA	0.065	0.068	NA	NA	NA	NA	0.328	NA	NA	0.987
```

# parameters
1.list all significant pairs files like
```bash
/data2/xiaohan/metasoft/v6_sig_egenes/GTEx_Analysis_v6p_eQTL/Brain_Anterior_cingulate_cortex_BA24_Analysis.v6p.signif_snpgene_pairs.txt.gz
/data2/xiaohan/metasoft/v6_sig_egenes/GTEx_Analysis_v6p_eQTL/Uterus_Analysis.v6p.signif_snpgene_pairs.txt.gz
```

2.list all nominal pairs files like
```bash
/data2/xiaohan/metasoft/v6_sig_egenes/GTEx_Analysis_v6p_eQTL/Brain_Anterior_cingulate_cortex_BA24_Analysis.v6p.signif_snpgene_pairs.txt.gz
/data2/xiaohan/metasoft/v6_sig_egenes/GTEx_Analysis_v6p_eQTL/Uterus_Analysis.v6p.signif_snpgene_pairs.txt.gz
```

3.parameter files like
contains 4 rows including this Dir Path, sig_list files, nominal_list files, and your output
```bash
#This is where the Java files Dir is
/data2/xiaohan/jar/metasoft
# This is a file contains where sig files are (file in parameters 1)
/data2/xiaohan/metasoft/testoutput/sigpairs.txt
# This is a file contains where nominal files are (file in parameters 2)
/data2/xiaohan/metasoft/testoutput/allpairs.txt
# This is your STUDY's output Dis
/data2/xiaohan/metasoft/testoutput
```


# usage
nohup java -jar MetaTissue.jar Meta_par.txt > log.txt & 
