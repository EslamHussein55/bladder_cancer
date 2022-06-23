'''
We may check for significance between expression levels by using the p value from a t test.  There are several problems with this methodology. T test assumes normality, and the genes do not have a normal distribution. Also, since we are testing so many genes, some p value will be small by chance. In some cases, genes only have a few nonzero values for males but all zeros for females, whi, 
Also, p value may be small because of limited sampling
To fix this,  use resampling. see if p values could be obtained with sample of males.  Compare with several male samples throw out p values that are not smallest, find distribution of remaining p values.
Compare this with result of using male sample as “female”. Out of these, plot max of cdf of p values.

Given a p value, what is the likelihood that this gene is significantly different?

Unfortunately, some of the genes that may have clinically significant differences may have only slight differences in level, so the differences may not be statistically significant. 
'''

eqVar = False
nSamples = 10
nCompares = 6000

pMaxFor_cdf = 0.01
cdfN = 20

# Remove all genes expressed fewer than this many times
minExpress = 5


altMF = 'less'  # Male is less
# altMF = 'greater' # Male is greater
# altMF = 'two-sided'

# switch sexes and test the other way
switchSex = False
if switchSex and (altMF == 'less'):
    altMF = 'greater'
if switchSex and (altMF == 'greater'):
    altMF = 'less'


import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
#!pip install scipy --upgrade # Needed for permutations

# For random shuffling
rng = np.random.default_rng()

# Computes CDF with evenly spaced points from data
def distFn(X,xMin,xMax, N):
    delta = (xMax-xMin)/N #size of intervals, [Xmin+j*delta,  Xmin+(j+1)*delta]
    X = np.maximum(xMin,np.minimum(X,xMax))
    Y = np.floor((X-xMin)/delta).astype(int) # 
    unique_elements, counts_elements = np.unique(Y, return_counts=True)
    pVec = np.zeros(N+1)
    pVec[unique_elements] = counts_elements
    Fvec = np.cumsum(pVec) #empirical cumulative distribution function
    xVec = xMin + delta*np.arange(N+1)
    pVec = pVec/delta
    Fvec = Fvec/Fvec[-1]
    return xVec,Fvec,pVec


Dataset = pd.read_csv(r'data_clinical_patient.txt', skiprows = 4, sep='\t')
ClinicalDF = pd.DataFrame(Dataset)
ClinicalDF = ClinicalDF.drop(['SUBTYPE','CANCER_TYPE_ACRONYM','OTHER_PATIENT_ID','DAYS_TO_BIRTH','ETHNICITY'], axis = 1)
fieldNames = list(ClinicalDF)
ClinicalDF = ClinicalDF.set_index('PATIENT_ID')
ClinicalPatientID = list(ClinicalDF.index.values)

Dataset = pd.read_csv(r'data_mrna_seq_v2_rsem.txt', sep='\t')
Dataset.Hugo_Symbol.fillna(Dataset.Entrez_Gene_Id, inplace = True)
del Dataset['Entrez_Gene_Id']
RSEMDF = Dataset.set_index('Hugo_Symbol')
RSEMDF.columns = RSEMDF.columns.str.replace('-01','')
geneIDs = np.array(RSEMDF.index.values)
RSEMPatientID = list(RSEMDF)

# Reconcile patient IDs between two datasets
ClinicalDF = ClinicalDF.drop(list(set(ClinicalPatientID) - set(RSEMPatientID)), axis = 0)
RSEMDF = RSEMDF.drop(list(set(RSEMPatientID) - set(ClinicalPatientID)), axis = 1)
ClinicalPatientID = list(ClinicalDF.index.values)
RSEMPatientID = list(RSEMDF)

geneVals = RSEMDF.to_numpy()

# remove genes that don't express enough times
geneValsSelect = (np.sum(geneVals>0,axis=1)>=minExpress)
geneVals = geneVals[geneValsSelect]
geneIDs = geneIDs[geneValsSelect]

ng = len(geneVals)

sexIx = ClinicalDF['SEX'].to_numpy()
if switchSex:
    maleIx = (sexIx=='Female')
else:
    maleIx = (sexIx=='Male')

femaleIx = (maleIx == False)
nMale = np.sum(maleIx)
nFemale = np.sum(maleIx)

geneValsMale = geneVals[:,maleIx].T
geneValsFemale = geneVals[:,femaleIx].T

# Store remaining number of genes
nRemainGenes = np.zeros((nSamples,nCompares))

# Store maximum cdf
cdfMax = np.zeros(cdfN+1)


# Store when gene is eliminated
geneRank = np.zeros(ng) + nCompares

# Do resampling to compare male with female.
for k in range(nSamples):
    print(k)
    geneVals0 = geneValsMale
    geneIx = np.arange(ng)

    # On last iteration, do with female data
    if k==(nSamples-1):
        geneVals2 = geneValsFemale
    #Otherwise, make random selection from male data
    else:
        tmp1 = rng.choice(nMale,nFemale)
        geneVals2 = geneValsMale[tmp1,:]

    # Compute "female" p-values
    tValsC,pValsC = stats.ttest_ind(geneVals0, geneVals2, alternative = altMF, equal_var = eqVar)

    for j in range(nCompares):
        if j%50==0:
            print(j)
        
        # Choose 'H0' sample for comparison
        tmp = rng.choice(nMale,nFemale)
        # Compute 'H0' p-values
        tValsO,pValsO = stats.ttest_ind(geneVals0, geneVals0[tmp,:],alternative = altMF, equal_var = eqVar)

        # Remove genes for which H0 p-values are smaller
        removeIx = np.logical_or(np.isnan(pValsO),  (pValsC >= pValsO))
        selectIx = np.logical_not(removeIx)

        # Reduce genes that are compared    
        geneVals0=geneVals0[:,selectIx]
        geneVals2 = geneVals2[:,selectIx]
        pValsC = pValsC[selectIx]

        # Record how long the gene lasted
        if k==(nSamples-1):
            geneRank[geneIx[removeIx]] = j
        
        
        geneIx = geneIx[selectIx]
        nRemainGenes[k,j] = np.sum(selectIx.astype(int)) 
    
        xVec,Ftmp,pTmp = distFn(pValsC,0,pMaxFor_cdf,cdfN)
        Ftmp = Ftmp*len(pValsC)
        if k<nSamples-1:
            cdfMax = np.maximum(Ftmp,cdfMax)
        

FemGenes = (nRemainGenes[-1,:].T).reshape((-1,1))
plt.plot( nRemainGenes[:-1,:].T) 
plt.plot(FemGenes,'k--')
plt.title("Number of genes remaining by iteration")
plt.yscale("log")
plt.show()

plt.plot(np.sort(nRemainGenes[:-1,-1].T));
plt.plot(FemGenes[-1] + 0*nRemainGenes[:-1,-1].T,'k--')
plt.title("Distribution of number of remaining genes")
plt.show()

plt.plot(xVec[:-1],cdfMax[:-1]);
plt.plot(xVec[:-1],Ftmp[:-1],'k--')
plt.title("Worst-case null versus observed cdf")
plt.show()


