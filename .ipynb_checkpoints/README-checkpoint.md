# Analysis of genetic data for male and female bladder cancer patients

We will be looking at genetic data from bladder cancer patients. The idea is to see if there is any relation between the expression levels of specific genes and patients' ability to recover from cancer.

### Downloading the data

The data to be analyzed is obtainable from the cBioPortal website (https://www.cbioportal.org/). The specific dataset we will be using is from this study: https://www.cbioportal.org/study/summary?id=blca_tcga_pan_can_atlas_2018.  

1. Download all data by clicking the `download` icon next to the blue underlined title `Bladder Urothelial Carcinoma (TCGA, PanCancer Atlas)` (see upper left of page). 

To unzip the data, you may need to install `7-zip` on your computer (https://www.7-zip.org/download.html). You actually have to extract twice, because the folder has been compressed twice:  first with `gzip` and then with `tar`.  

The files that you will need are:

* `data_clinical_patient.txt`  (in the main folder)
* `data_mrna_seq_v2_rsem.txt`  (in the main folder)

RSEM stands for "RNA-Seq by Expectation Maximization". The two RSEM files gives measured gene expression levels for all of the 20,000+ genes on the human genome.  For background information about RSEM, see https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-323.  The RSEM file has data from badder cells from about 400 bladder cancer patients. 
