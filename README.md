# qgrs-Rwrapper
A R wrapper for qgrs-cpp. The latter is a C++ implementation of QGRS mapping algorithm published by freezer333 on https://github.com/freezer333/qgrs-cpp. This project support a R wrapper that allow analysis possible G4 structures in transcript from a specific gene list. It also provide a executable file of qgrs-cpp on win64.

The original papers:

- Frees, S., Crum, M., Menendez, C., Bagga, P.S. (2014) ["QGRS-Conserve: A Computational Method for Discovering Evolutionarily Conserved G-quadruplex Motifs"](https://humgenomics.biomedcentral.com/articles/10.1186/1479-7364-8-8). Human Genetics, BioMed Central, Vol 8 (8)
- Menendez, C., Frees, S., and Bagga, P. (2012) [QGRS-H Predictor: A Web Server for Predicting Homologous Quadruplex forming G-Rich Sequence Motifs in Nucleotide Sequences](https://academic.oup.com/nar/article/40/W1/W96/1074452/QGRS-H-Predictor-a-web-server-for-predicting). Nucleic Acids Res. doi: 10.1093/nar/gks42240: W96-W103.

## Usage

```
Rscript G4_analysis.R [regionType] [gene.csv] [qgrs_exec_filepath]
```

The code hasbeen tested on Windows.

## Arguments

The following command line options are supported:

```
regionType: fiveUTR, threeUTR, cds
gene.csv: input file, example in gene.csv
qgrs_exec_filepath: qgrs executable file path, default is The built-in qgrs executable, please specific it if it not work
```

## Output
By default, output is in res folder, temporary files is in tmp folder. The implication of column are showed on column titles, the last 8 columns are the output for qgrs-cpp, please refer to https://github.com/freezer333/qgrs-cpp for more informations.


## Requirements
```
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
```

