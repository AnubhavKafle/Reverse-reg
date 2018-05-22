
#! usr/bin/python

import numpy as np
from scipy.stats import chi2
from statsmodels.sandbox.stats.multicomp import fdrcorrection_twostage as fdr
import time
import argparse
import matplotlib.pyplot as plt
from main.iotools.dosageparser import DosageParser
#from main.iotools.readVCF import ReadVCF 
from main.inference.regulizer_optimizer import OptimizeRegularizer


def parse_args():

    parser = argparse.ArgumentParser(description='Trans-Eqtls from Joint Association AnalysiS (TEJAAS)')

    parser.add_argument('--genotype',
                        type=str,
                        dest='genotype_filename',
                        metavar='FILE',
                        help='input genotype file')

    parser.add_argument('--transgeno',
                        type=str,
                        dest='trans_genofile',
                        metavar='FILE',
                        default = " ",
                        help='input potential trans-eQTLs genotype filename for optimization.')

    parser.add_argument('--optimize',
                        type=int,
                        dest='optimize_sigbeta',
                        default = 0,
                        help='1/0 option to optimize sigmabeta')

    parser.add_argument('--sigbeta',
                        type=float,
                        dest='sigbeta',
                        default = 0.006,
                        help='User given sigbeta')
    
    parser.add_argument('--sample',
                        type=str,
                        dest='sample_filename',
                        metavar='FILE',
                        help='input fam file')

    parser.add_argument('--expression',
                        type=str,
                        dest='expression_filename',
                        metavar='FILE',
                        help='input expression file')

    parser.add_argument('--output',
                        type=str,
                        dest='output_fileprefix',
                        metavar='FILE',
                        help='output file prefix')

    parser.add_argument('--start',
                        type=int,
                        dest='startsnp',
                        help='starting SNP index')

    parser.add_argument('--end',
                        type=int,
                        dest='endsnp',
                        help='ending SNP index')

    opts = parser.parse_args()
    return opts


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

def norm_binom(genotype, f):
    genotype = (genotype - (2 * f)) / np.sqrt(2 * f * (1 - f))
    return genotype

    #Read genotype here
def parse_geno (genotype_filename, sample_filename, startsnp, endsnp):
	ds        = DosageParser(genotype_filename, sample_filename, startsnp, endsnp)
	dosage    = ds.dosage
	snpinfo   = ds.snpinfo
	donorids  = ds.sample_id
	nsnps     = ds.nsnps
	nsample   = ds.nsample
	return dosage, snpinfo, donorids, nsnps, nsample

#quality check, matching
def quality_check ( sampleids, donorids):
	choose_ids     = [x for x in sampleids if x in donorids]
	dosage_indices = [i for i, x in enumerate(donorids)  if x in choose_ids]
	exprsn_indices = [i for i, x in enumerate(sampleids) if x in choose_ids]
	return dosage_indices, exprsn_indices

#Read python variables here

opts                = parse_args()
out_fileprefix      = opts.output_fileprefix
genotype_filename   = opts.genotype_filename
sample_filename     = opts.sample_filename
expression_filename = opts.expression_filename
startsnp            = opts.startsnp
endsnp              = opts.endsnp
optim               = opts.optimize_sigbeta
user_sigbeta        = opts.sigbeta
transgeno           = opts.trans_genofile

#Read gene expression here 
sampleids, expression, gene_names = read_expression(expression_filename)

'''
       Read dosage maf > 0.1 and < 0.9, monoalleleic SNPs. \  
       Check the Dosageparser file before doing this as there \ 
       is dependencies on MAF.  
'''
#Check if optimization is demanded.
if optim:
    print ("\nSigma beta optimization started. Reading from the provided trans-genotype file.")
    tic = time.time()
    try:
        trans_dosage, trans_snpinfo, trans_donorids, trans_nsnps, trans_nsample = parse_geno (transgeno, sample_filename, 0, 50000)
        dosage_indices, exprsn_indices = quality_check (sampleids , trans_donorids)
        conv_dosage = np.round(trans_dosage)  #Best performance if dosages are rounded off
        geno        = conv_dosage[:, dosage_indices].T 
          # for the matrix to be of NXI size
        freq        = np.array([x.freq for x in trans_snpinfo])
        expr        = expression[:, exprsn_indices]
        geno        = norm_binom(geno, freq)
        geno        = geno.T

        optimize_sigbeta   = OptimizeRegularizer(geno, expr)
        optimize_sigbeta.update()
        optimize_sigbeta.niter
        sigbeta            = optimize_sigbeta.sigmareg
        toc = time.time()
        print ("Sigma beta optimization completed in :", toc - tic , "seconds\n")
        print ("Optimized sigma beta value is:" , sigbeta,"\n")
        del geno, conv_dosage, trans_dosage
    except OSError as e:
        print("\n",e, ". Trans-eQTLs file not provided for Optimization\n")
        raise SystemExit


else:
    print("\n=============================================================")
    print("\nsigma beta optimization not requested")
    if user_sigbeta == 0.006:
        print ("\nsigma beta not provided. Using default value of 0.006")
    else:
        print("\nsigma beta value provided as : " , user_sigbeta)
    print("=============================================================")
    sigbeta         = user_sigbeta


tic = time.time()

print ("\nReading Genotype")
dosage, snpinfo, donorids, nsnps, nsample = parse_geno ( genotype_filename, sample_filename, startsnp, endsnp)
dosage_indices, exprsn_indices = quality_check (sampleids , donorids)
conv_dosage = np.round(dosage)  #Best performance if dosages are rounded off
geno        = conv_dosage[:, dosage_indices].T   # for the matrix to be of NXI size
freq        = np.array([x.freq for x in snpinfo])

# Scaling and normalization
geno  = norm_binom(geno, freq)  #commented out to later check it 
expr  = expression[:, exprsn_indices]
geno  = geno.T  # SXN
print ("Completed data loading and processing\n")


randoms = np.random.choice(geno.shape[0], 30000)
genos = geno[randoms,:]

#Calculate Qscores here
sigma_geno = 1
U, S, V_t = np.linalg.svd(np.transpose(expr),full_matrices=False)
fact = sigma_geno**2 / sigbeta**2
diagw = np.square(S) / (np.square(S) + fact)
Wsvd = np.dot(U, np.dot(np.diag(diagw), U.T))
Qsvd  = np.diag(np.dot(genos, np.dot(Wsvd,  genos.T)))

#Inference parameters
# alpha_numer = 0
# alpha_denom = 0

# for i in range(S.shape[0]):
#     alpha_numer += S[i]**4/(S[i]**2 + sigma_geno**2/sigbeta**2)**2
#     alpha_denom += S[i]**2/(S[i]**2 + sigma_geno**2/sigbeta**2)

#alpha_numer = np.sum(S**4 / (S**2 + fact)**2 )
#alpha_denom = np.sum(S**2 / (S**2 + fact))

#alpha = alpha_numer/alpha_denom

#f_degree = (U.shape[0]+2)/(U.shape[0]-1) * (alpha_denom**2/ alpha_numer)

alpha = np.square(np.std(Qsvd)) / (2* np.mean(Qsvd))
f_degree = 2* np.square(np.mean(Qsvd)) / np.square(np.std(Qsvd))

Q_pval = 1 - chi2.cdf(Qsvd/alpha,f_degree)

toc = time.time()

print( "\nTime taken for the Qscores calculation:", toc-tic, "seconds\n")
#############################################################################################################################
# Plotting here

Qnull = np.random.chisquare(f_degree, Qsvd.shape[0])
Qnull_sort = Qnull[np.argsort(Qnull)]
Qsvd_sort = Qsvd[np.argsort(Qsvd)]/alpha

NC = expression_filename.split("_")[-1].split(".")[0]

output_real = "/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qsvd_PC{}.npy".format(NC)
output_null = "/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qnull_PC{}.npy".format(NC)
np.save(output_real, Qsvd_sort)
np.save(output_null, Qnull_sort)

print("+++++++++++++++++++++++++++++++++++++++++")

print("\nWritten Qvalues for PC analysis:{}".format(NC))

print("\n++++++++++++++++++++++++++++++++++++++++")
