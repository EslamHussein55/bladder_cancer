'''
Test 1:  reserve some male genes as "female", and see if the algorithm falsely identifies differences
'''

def testAlgFn(shuffleData, impersonateWithMale,sampleWithoutReplacement,eqVar,nSamples,nCompares,pMaxFor_cdf,cdfN,altMF,switchSex,thisSeed, ClinicalDF,geneVals):

    if switchSex and (altMF == 'less'):
        altMF = 'greater'
    if switchSex and (altMF == 'greater'):
        altMF = 'less'
    
    ng = len(geneVals)

    import numpy as np
    from scipy import stats
    import warnings
    warnings.filterwarnings("ignore")
    #!pip install scipy --upgrade # Needed for permutations
    
    # For random shuffling
    rng = np.random.default_rng(seed=thisSeed)
    
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
    
    
    
    sexIx = ClinicalDF['SEX'].to_numpy()
    if switchSex:
        maleIx = (sexIx=='Female')
    else:
        maleIx = (sexIx=='Male')
    
    femaleIx = (maleIx == False)
    nMale = np.sum(maleIx)
    nFemale = np.sum(femaleIx)
    
    geneValsMale = geneVals[:,maleIx].T

    
    # Do with real female data
    geneValsFemale = geneVals[:,femaleIx].T
    
    # Shuffle because there appears to be ordering within males
    if shuffleData: 
        rng.shuffle(geneValsMale,axis=0)
        
    # Option 1:  Select male as simulated "female" and run 
    if impersonateWithMale:
        geneValsFemale = geneValsMale[:nFemale,:]
    
    # Option 2: Reduce males so it is without replacement
    if sampleWithoutReplacement:
        geneValsMale = geneValsMale[nFemale:,:]
    
    nMale = np.shape(geneValsMale)[0]
    
    # Store remaining number of genes
    nRemainGenes = np.zeros((nSamples,nCompares))
    
    # Store maximum cdf
    cdfMax = np.zeros(cdfN+1)
    
    
    # Store when gene is eliminated
    geneRank = np.zeros(ng) + nCompares
    
    # Do resampling to compare male with female.
    for k in range(nSamples):
        if nSamples > 1:
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
            if j%6000==0 and j>0:
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
    
    return nRemainGenes,FemGenes,xVec,cdfMax,Ftmp