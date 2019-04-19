# JPR-201712_MS2-MS3_PSM-peptide-protein

This is publically available data from this publication:
> D’Angelo, G., Chaerkady, R., Yu, W., Hizal, D.B., Hess, S., Zhao, W., Lekstrom, K., Guo, X., White, W.I., Roskos, L. and Bowen, M.A., 2017. Statistical models for the analysis of isobaric tags multiplexed quantitative proteomics. Journal of proteome research, 16(9), pp.3124-3136.

## First analysis (June 15, 2018): 

### [Direct link](https://pwilmart.github.io/TMT_analysis_examples/MS2MS3_peptides_proteins.html) to HTML view of the notebook.

**Click on the Jupyter notebook file (_MS2MS3_peptides_proteins.ipynb_) to see the notebook in your browser. It may take a minute for the page to render and load, so please be patient.**


In shotgun (bottom-up) proteomics experiments there can be a lot of data redundancy. There can be multiple MS2 scans acquired for the same analyte. There can be multiple PSM (the sequences matched to the MS2 scans) to the same sequence form (a peptide sequence in a given modification state). There can be multiple peptide forms of the same peptide sequence. There can be (and usually are) multiple peptides from the same protein.

This leads to many questions. Is it better to have lots of lower quality data points? Or is it better to try and aggregate the data into fewer, higher-quality measures? What level of aggregations are better? What kind of aggregation operation to use?
This notebook will explore the simplest aggregation operation of summation. It will look at some data distribution properties at three levels of aggregation: PSMs, peptides (a peptide sequence in a given modification state), and proteins (using the razor peptide approach to shared and unique peptides).

One fun, extra possibility with the dataset from the publication below is to compare the same sample run on a more traditional MS2 reporter ion method, and the newer synchronous precursor scan MS3 method. The dataset does not lend itself to an analysis of accuracy and dynamic range, but we can see if the duty cycle or signal intensities are very different between the two ways of acquiring reporter ions signals.
