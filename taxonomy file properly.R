tax_file <- scan(file="trainset14_032015_rmdup.tax", what="", sep="\n", quiet=TRUE)

accession <- gsub("^(\\S*).*", "\\1", tax_file) #some are separated by tabs or spaces or both

taxonomy <- gsub(".*(Root.*)", "\\1", tax_file)
taxonomy <- gsub(" ", "_", taxonomy)    #remove spaces and replace with '_'
taxonomy <- gsub("\t", "", taxonomy)    #remove extra tab characters
taxonomy <- gsub("[^;]*_incertae_sedis$", "", taxonomy)
taxonomy <- gsub('\"', '', taxonomy) #remove quote marks

levels <- read.table(file="trainset14_db_taxid.txt", sep="*", stringsAsFactors=FALSE)
subs <- levels[grep("sub", levels$V5),]
sub.names <- subs$V2

tax.split <- strsplit(taxonomy, split=";")

remove.subs <- function(tax.vector){
  return(tax.vector[which(!tax.vector %in% sub.names)])
}

no.subs <- lapply(tax.split, remove.subs)
no.subs.str <- unlist(lapply(no.subs, paste, collapse=";"))
no.subs.str <- gsub("^Root;(.*)$", "\\1;", no.subs.str)

write.table(cbind(as.character(accession), no.subs.str), "trainset14_032015.rdp.tax", row.names=F, col.names=F, quote=F, sep="\t")

