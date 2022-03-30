message('Will install necessary for G4_analysis.R, please wait a moment')
if(!require('BiocManager')){
    message('BiocManager is not installed, install it')
    install.packages('BiocManager')
}


options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

require_packages <- c('GenomicRanges', 'BSgenome.Hsapiens.UCSC.hg38', 
                      'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db', 
                      'dplyr')
message(paste("check and install packages if need ", require_packages, collapse = ', '))
BiocManager::install(require_packages)
