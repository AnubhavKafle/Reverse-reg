#!/bin/bash

EXPR=$1

python /home/anubhavk/Desktop/Internship_Qstat/codebase/transeqtl.py --genotype /home/anubhavk/Desktop/Internship_Qstat/data/GTEx.VCF/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr1.gz  --sample /home/anubhavk/Desktop/Internship_Qstat/data/GTEx.VCF/donor_ids.fam --expression ${EXPR} --output pp --start 0 --end 100000  --transgeno /home/anubhavk/Desktop/Internship_Qstat/data/TEJAAS_Gtex_lmcorr_empirical/best_Qscores_lmcorr_Gtex_empirical/Genotype_filtered/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_filteredTEJAAS_transeqtls.gz --optimize 0 --sigbeta 0.0138262336815
