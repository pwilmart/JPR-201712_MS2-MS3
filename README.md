# JPR-201712_MS2-MS3

This is publicly available data from this publication:
> D’Angelo, G., Chaerkady, R., Yu, W., Hizal, D.B., Hess, S., Zhao, W., Lekstrom, K., Guo, X., White, W.I., Roskos, L. and Bowen, M.A., 2017. Statistical models for the analysis of isobaric tags multiplexed quantitative proteomics. Journal of proteome research, 16(9), pp.3124-3136.

Part of that study was 13 human proteins spiked into a common E. coli background. Samples were labeled with 10-plex TMT. The samples were run on a Thermo Fusion instrument in two experiments. One experiment produced the reporter ion signals in the MS2 scans in a more conventional strategy. The other experiment used the newer SPS MS3 method available on the Tribrid platforms. If the human proteins are excluded, this is essentially two sets of 10 technical replicates of E. coli lysates.

An experiment comparing commercial serum systemic lupus erythematosus (SLE) data to normal serum using 10-plex TMT labeling analyzed on a Q-Exactive platform was also part of the study. The RAW files were downloaded and processed with Comet/PAW.

## First analysis (June 15, 2018)

### [Direct link](https://pwilmart.github.io/TMT_analysis_examples/MS2MS3_peptides_proteins.html) to an HTML view of the notebook.

### Folder: `first_analysis_20180615`

The first analysis is described in the README.md file in the above folder. It explores how similar the identical data is at the PSM, peptide, or proteins level. The analysis was done with MaxQuant.

## Second analysis (April 19, 2019)

### [Direct link](https://pwilmart.github.io/TMT_analysis_examples/JPR-2017_E-coli_MS2-MS3.html) to an HTML view of the notebook.

### Folder: 'second_analysis_20190419'

The second analysis uses the PAW pipeline to directly compare the MS2 data to the MS3 data.

## Serum analysis (April 28, 2019)

### [Direct link](https://pwilmart.github.io/TMT_analysis_examples/JPR-2017_serum.html) to an HTML view of the notebook.

### Folder: 'serum_20190428'

Added analysis of the Q-Exactive SLE serum experiments using the PAW pipeline. Depleted serum is a little different kind of sample compared to cell lysates. Choice of search parameter needs more thought. Contaminant designations have to be expanded to include the depleted protein families. Some of the practical aspects of biological variability in human samples and smaller datasets are illustrated.
