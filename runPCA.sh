#!/bin/bash

source activate py36

for j in {1..15}; 
do
    #COVAR="/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.expr/GTEx_WB_Covar_PC"${j};
    EXPR="/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.expr/GTEx_normalized_WB_expr_NoPEER_lmcorrected_gene_genetic_PCA_"${j}".txt";
    #GENETIC_3="/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.expr/3genetic_PC.txt";
    #PC="/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.expr/GTEx_wholeblood_PC${j}.txt";
    #python /home/anubhavk/Desktop/Internship_Qstat/codebase/PCA_Gene_expr.py --components ${j};
    #cat ${PC} ${GENETIC_3} > ${COVAR};
    #Rscript /home/anubhavk/Desktop/Internship_Qstat/codebase/lib/Gene_expr_residual.R ${COVAR} ${j};
    bash /home/anubhavk/Desktop/Internship_Qstat/codebase/run_transeqtl.sh ${EXPR};
done
