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
java -jar Meta_par.txt 