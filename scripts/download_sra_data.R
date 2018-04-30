#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1){
    write('To Run: Rscript download_sra_data.R <file>. Where file is path to a file containing SRP ids one per line', stderr())
    quit(status=1)

}
suppressMessages(library(SRAdb))
suppressMessages(library(DBI))
suppressMessages(library(GEOmetadb))


###### User defined #########
ascpKeyPath <- '~/.aspera/connect/etc/asperaweb_id_dsa.openssh'
download.root.location <- '/staging/as/wenzhenl/re-ribo-data' #'/home/cmb-06/as/skchoudh/dna/re-ribo-data'
######              #########


getInfo <- function(accession, sra_con, sraType='sra') {
  return (getSRAinfo(in_acc=accession, sra_con = sra_con, sraType=sraType)[,1])
}

getSampleAlias <- function(srs, sra_con) {
  return (unlist(dbGetQuery(sra_con, paste("SELECT sample_alias, scientific_name",
                                           " FROM sample where sample_accession='", srs, "' ", sep="") ) ) )
}

getExperimentAlias <- function(srx, sra_con) {
    return (unlist(dbGetQuery(sra_con, paste("SELECT experiment_alias",
                                             " FROM experiment WHERE experiment_accession='", srx, "'", sep="") ) ) )
}




################### Do not edit below this line #########


SRAsqlfile <- 'SRAmetadb.sqlite'
GEOsqlfile <- 'GEOmetadb.sqlite'
ascpCMD <- paste0('ascp -k 1 -QT -l 2000m -i ', ascpKeyPath)
if(!file.exists(SRAsqlfile)) SRAsqlfile <<- getSRAdbFile()
if(!file.exists(GEOsqlfile)) GEOsqlfile <<- getSQLiteFile()

input_list <- scan(args[1], what="", sep='\n')

sra_con <- dbConnect(SQLite(), SRAsqlfile)
geo_con <- dbConnect(SQLite(), GEOsqlfile)

get_SRP_from_GSE <- function(gse){
    if (grepl('SRP', gse)){
        return (gse)
    }
    results <- as.list(dbGetQuery(geo_con, paste0("SELECT * FROM gse WHERE gse = '", gse, "'")))
    supplementary_file <- strsplit(results$supplementary_file, ';')[[1]]
    srp_match <- grep('SRP', supplementary_file, value=TRUE)
    return (sapply(strsplit(srp_match, '/'), tail, 1))

}

download_metadata <- function(srp, destination_dir){
    url <- paste0('https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=', srp, '&get=csv_all')
    cmd <- paste0('wget -O ', file.path(destination_dir, 'SraRunTable.txt'), ' -c ', url)
    #print(cmd)
    #system(cmd)
    return (url)
}
getMetadata <- function(accession, sra_con, cols){
    cols <- paste(cols,  collapse=", ")
    return (dbGetQuery(sra_con, paste("SELECT ", cols,
                                      " FROM sra WHERE run_accession='", run_accession, "'", sep="")))
}


download <- function(x){
  srx_id <- x[2]
  srp_id <- x[1]
  print(srp_id)
  #sra_id <- x[2]
  #url <- x[6]
  #filename <- paste(sra_id, 'sra', sep='.')
  dir.create(file.path(download.root.location), showWarnings = F)
  dir.create(file.path(download.root.location, srp_id), showWarnings = F)
  #download_metadata(srp_id, file.path(download.root.location, srp_id))
  dir.create(file.path(download.root.location, srp_id, srx_id), showWarnings = F)
  dir_name <- file.path(download.root.location, srp_id, srx_id)
  getSRAfile(srx_id, sra_con, fileType = 'sra',
             srcType = 'fasp', destDir= file.path(download.root.location, srp_id, srx_id),
             ascpCMD = ascpCMD )
}

srp_list <- c()
for (element in input_list) {
    if (grepl('GSE*', element)){
        srp_list <- c(srp_list, get_SRP_from_GSE(element))
    }
    else{
        srp_list <- c(srp_list, element)
    }

}

for (srp in srp_list){
    print(srp)
    conversion <- sraConvert(srp, sra_con = sra_con)
    conversion$ftp <- apply(matrix(conversion$run, ncol=1), MAR=1, FUN = getInfo, sra_con = sra_con)
    ## Correct ftp link
    conversion$ftp <- gsub("ByExp", "ByRun", conversion$ftp)

    #conversion$experiment_alias <- cbind(conversion, t(apply(matrix(conversion$experiment, ncol=1), MAR=1, FUN=getExperimentAlias, sra_con=sra_con)))
    conversion <- cbind(conversion, t(apply(matrix(conversion$sample, ncol=1), MAR=1, FUN=getSampleAlias, sra_con=sra_con)))
    #srx_list <- apply(conversion, 2, unique)$experiment

    #sra_df <- as.data.frame(listSRAfile(srx_list, sra_con))
    cols <- c("library_name", "library_strategy", "library_source", "library_selection", "library_layout", "experiment_attribute", "experiment_alias", "experiment_title",
                         "sample_name", "read_spec", "spots", "sample_attribute", "study_title", "study_accession", "common_name", "taxon_id")
    metadata <- as.data.frame(matrix(NA, nrow(conversion), length(cols)+2))
    for (i in 1:nrow(metadata)) {
    run_accession <- conversion$run[i]
    metadata[i, 1:length(cols) ] <- getMetadata(run_accession, sra_con, cols)
    metadata[i, length(cols)+1] <- conversion$sample_alias[i]
    metadata[i, length(cols)+2] <- conversion$scientific_name[i]

    }
    names(metadata) <- c(cols, 'sample_alias', 'scientific_name')
    metadata <- metadata[!duplicated(metadata), ]
    dir.create(file.path(download.root.location, srp), showWarnings = F)
    write.table(metadata, file.path(download.root.location, srp, 'metadata.tsv'), row.names=FALSE, sep='\t')
    sra_df <- as.data.frame(conversion)
    sra_df <- sra_df[!duplicated(sra_df),]
    srp_srx = sra_df[, c(1,4)]
    sra_srx <- srp_srx[!duplicated(srp_srx),]
    apply(sra_srx, 1, download)
}
dbDisconnect(sra_con)
dbDisconnect(geo_con)

