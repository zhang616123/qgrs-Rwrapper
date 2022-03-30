#author: zhang616123
#date: 20220330
#run qgrs, a tools for G4 predict, with specific gene list

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  warning('no specific command arguments, will run with default, set region type as fiveUTR \n
          command usage is Rscript G4_analysis.R [regionType] [gene.csv] [qgrs_exec_filepath]')
  regionType <- 'fiveUTR'
  Sys.sleep(3)
}else {
  regionType <- args[1]
}

message('run qgrs, a tools for G4 predict, with specific gene list\n 
      author: zhang616123, more informations about qgrs output, please see https://github.com/freezer333/qgrs-cpp')
Sys.sleep(3)


if(length(args) == 2){
  gene_list1 <- args[2]
} else {
  gene_list1 <- 'gene.csv'
}

if(Sys.info()['sysname'] == "Windows" && length(args) < 3){
  qgrsExecPath = "./qgrs-cpp-master/qgrs-cpp-master/qgrs.exe"
} else {
  warning('operating system is not Windows, will try qgrs executable file for Linux (Ubuntu), 
          please specific executable file path if it not work')
  qgrsExecPath = "./qgrs-cpp-master/qgrs-cpp-master/qgrs"
}

if(length(args) == 3){
  qgrsExecPath <- args[3]
  
}


geneList <- read.csv(gene_list1, 
                     header = F)
colnames(geneList) <- c('ENSEMBL', 'SYMBOL')

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
if(regionType == "fiveUTR"){
  trs <- fiveUTRsByTranscript(txdb, use.names=TRUE)
  filePath1 <- "res/fiveUtrQgrs.tab"
  logPath1 <- "res/fiveUtrQGRS.log"
} else if(regionType == "threeUTR"){
  trs <- threeUTRsByTranscript(txdb, use.names=TRUE)
  filePath1 <- "res/threeUtrQgrs.tab"
  logPath1 = "res/threeUtrQGRS.log"
} else if(regionType == "cds") {
  trs <- cdsBy(txdb, use.names=TRUE)
  filePath1 <- "res/cdsQgrs.tab"
  logPath1 <- "res/cdsQGRS.log"
} else {
  stop(paste0("inputed argument ", regionType, "is not accepected, only accepecte fiveUTR, threeUTR and cds"))
}

genes <- genes(txdb)

trans_id <- AnnotationDbi::select(org.Hs.eg.db, geneList$ENSEMBL, 'ENSEMBLTRANS', 
                                  'ENSEMBL') %>% 
  na.omit() %>% 
  left_join(geneList)

#names(trs) <- gsub("\\..*$", '', names(trs))

gene_trs <- trs[which(gsub("\\..*$", '', names(trs)) %in% trans_id$ENSEMBLTRANS)]

#trans_grs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gene_trs)

if(!dir.exists("tmp"))dir.create("tmp")
if(!dir.exists("res"))dir.create("res")

write_file <- function(strs, filePath, append = TRUE){
  #fileConn <- file(filePath)
  write(strs, filePath, append=append)
  #close(fileConn)
}

get_geneid <- function(transId, orgDb = org.Hs.eg.db, 
                       logPath = NULL){
  trans_id <- gsub("\\..*$", '', transId)
  trans_id <- suppressMessages(AnnotationDbi::select(orgDb, trans_id, 
                                    c('ENSEMBL', 'SYMBOL'), 'ENSEMBLTRANS'))
  if(dim(trans_id)[1] > 1){
    warning(paste0(transId, " return 1:many ENSEMBL or SYMBOL mapped, 
                                      Its not be expected\n", trans_id))
    warning(paste0(transId, " will use the first row of mapped res"))
    if(!is.null(logPath)){
      write_file(paste0(transId, " return 1:many ENSEMBL or SYMBOL mapped, 
                                      Its not be expected, will use first row of mapped res\n", trans_id), 
               file = logPath)
    }
    trans_id <- trans_id[1, ]
  }
  #print(trans_id)
  paste(trans_id, collapse = "\t")
}


write_file("ENSEMBLETRANSWITHVERSION\tENSEMBLTRANS\tENSEMBL\tSYMBOL\tTARGETSEQ\tID\tT1\tT2\tT3\tT4\tTS\tGS\tSEQ", 
           filePath1, 
           append = FALSE)
write_file("", 
           logPath1, 
           append = FALSE)

run_qgrs <- function(gr, trans_id, 
                     qgrsExecPath = "./qgrs-cpp-master/qgrs-cpp-master/qgrs.exe", 
                     grepExecPath = "./qgrs-cpp-master/gnugrep/grep.exe", 
                     outputPath = filePath1, 
                     logPath = logPath1
                     ){
  gr  <- gr[trans_id]
  trans_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  trans_seqFull <- lapply(trans_seq, as.character) %>% 
    unlist() %>% 
    paste(sep = "", collapse = "")
  tmpFile <- paste0("tmp/", trans_id, "_seq.txt")
  fileConn<-file(paste0("tmp/", trans_id, "_seq.txt"))
  writeLines(trans_seqFull, fileConn)
  close(fileConn)
  print(trans_id)
  #cmd <- paste0("powershell", qgrsExecPath, " -i ", tmpFile, " | ", grepExecPath, " -v ID")
  #print(cmd)
  cmd_res <- system(paste0(qgrsExecPath, " -i ", tmpFile), 
         intern = TRUE)%>%
    as.character() %>% 
    strsplit("\n") %>% 
    unlist()
  print(cmd_res)
  cmd_res <- cmd_res[c(-1, -2)]
  #print(cmd_res)
  if(grepl('No QGRS found', cmd_res[1])){
    warning(paste0(trans_id, 'No QGRS found'))
    return(NA)
  } else {
    #print(cmd_res)
    sapply(cmd_res, function(res, id = trans_id, seq = trans_seqFull, 
                             filePath = outputPath){
      gene_id <- get_geneid(id, logPath = logPath)
      res <- gsub("\\s{1,10}", "\t", res)
      res_str <- paste(trans_id, gene_id, seq, res, sep = "\t", collapse = "\t")
      write_file(res_str, filePath = outputPath)
    })
  }
}

running_res <- sapply(names(gene_trs), run_qgrs, gr = gene_trs)
message(paste0("analysis for ", regionType, " sucessful, results is in res folder, ", 
               "tmp is temporary files and can be delete safely"))
