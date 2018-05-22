from sklearn.decomposition import PCA
import subprocess
import numpy as np
from main.iotools.dosageparser import DosageParser
import argparse

def parse_args():

	parser = argparse.ArgumentParser(description='Take out PC from gene expression ')

	parser.add_argument('--components',
						type=int,
						dest='pc_num',
						help='No. of Principle components')

	opts = parser.parse_args()
	return opts

arguess = parse_args()

NC = arguess.pc_num


expression_filename = "/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.expr/gtex_wholeblood_normalized.expression.txt"
genotype_filename   = "/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.VCF/GTEx_450Indiv_genot_imput_info04_maf01_HWEp1E6_dbSNP135IDs_donorIDs_dosage_chr1.gz"
sample_filename     = "/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.VCF/donor_ids.fam"
startsnp = 1
endsnp = 10  

#Functions required
def read_expression(filename):
    gene_names = list()
    expression = list()
    with open(filename, 'r') as genefile:
        header = genefile.readline()
        donorids = header.strip().split("\t")[1:]
        for line in genefile:
            linesplit = line.strip().split("\t")
            expression.append(np.array(linesplit[1:], dtype=float))
            gene_names.append(linesplit[0])
    expression = np.array(expression)
    return donorids, expression, gene_names

#Lazy ass code to match the sample IDs already with genotype for no problems later. So reading few genotypes
ds = DosageParser(genotype_filename, sample_filename, startsnp, endsnp)
dosage = ds.dosage 
snpinfo = ds.snpinfo
donorids = ds.sample_id
nsnps = ds.nsnps
nsample = ds.nsample

#gene expression 
sampleids, expression, gene_names = read_expression(expression_filename)

#quality check, matching 
choose_ids = [x for x in sampleids if x in donorids]
exprsn_indices = [i for i, x in enumerate(sampleids) if x in choose_ids]

#commented out to later check it 
#geno = (geno - np.mean(geno,axis=1).reshape(-1,1)) / np.std(geno,axis = 1).reshape(-1,1) # scaling and normalizing
expr = expression[:, exprsn_indices]


# PCA gene expression 
pca = PCA(n_components=NC)
pca.fit(np.transpose(expr))
expr_pca = pca.transform(np.transpose(expr))
PC = expr_pca.T

Outputfile  = "/home/anubhavk/Desktop/Internship_Qstat/data/GTEx.expr/GTEx_wholeblood_PC{}.txt".format(NC)

with open(Outputfile, "w") as outfile:
    ids = "\t".join([str(i) for i in sampleids])
    header = "ID" "\t" + ids + "\n"
    outfile.write(header)
    for i in range(PC.shape[0]):
        a = "\t".join([str(j) for j in PC[i,:]])
        b = "PCA{}".format(i) + "\t" + a + "\n"
        outfile.write(b)

print("written PCA covar for", Outputfile)
