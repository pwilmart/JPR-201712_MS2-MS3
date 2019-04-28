# JPR-201712_serum

This is publicly available data from this publication:
> D’Angelo, G., Chaerkady, R., Yu, W., Hizal, D.B., Hess, S., Zhao, W., Lekstrom, K., Guo, X., White, W.I., Roskos, L. and Bowen, M.A., 2017. Statistical models for the analysis of isobaric tags multiplexed quantitative proteomics. Journal of proteome research, 16(9), pp.3124-3136.

An experiment comparing commercial serum systemic lupus erythematosus (SLE) data to normal serum using 10-plex TMT labeling analyzed on a Q-Exactive platform was also part of the study. The RAW files were downloaded and processed with Comet/PAW.

## Serum analysis (April 28, 2019)

### [Direct link](https://pwilmart.github.io/TMT_analysis_examples/JPR-2017_serum.html) to an HTML view of the notebook.

### Folder: 'serum_20190428'

Added analysis of the Q-Exactive SLE serum experiments using the PAW pipeline. Depleted serum is a little different kind of sample compared to cell lysates. Choice of search parameter needs more thought. Contaminant designations have to be expanded to include the depleted protein families. Some of the practical aspects of biological variability in human samples and smaller datasets are illustrated.

The protein summary report (and its peptide companion) has this progression:

1. protein_summary_9.txt - redundant protein summary
1. grouped_protein_summary_9.txt - non-redundant, grouped protein summary
1. grouped_protein_summary_TMT_9.txt - with added reporter ion intensities
1. ave_labeled_grouped_protein_summary_TMT_9.txt - added sample labels and mock reference
1. ave_labeled_grouped_protein_summary_TMT_9_IRS_normalized.txt - after IRS norm script
1. JPR-2017_serum_IRS_edgeR_annotated.xlsx - with added edgeR results and annotations
