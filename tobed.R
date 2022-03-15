# Loading packages
foo <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE , quietly = TRUE) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE , quietly = TRUE)
    }
  }
}
foo( c("argparser" , "roperators", "dplyr", "foreach","doParallel") )

p <- arg_parser("MeH result file to .bedGraph")

# Add command line arguments
p <- add_argument("-m","--Meh", help="input Meh resultes csv file", type="character")
p <- add_argument("-o", "--output", help="output to bedGraph file")
p <- add_argument("-r", "--reverse", default="all" ,help="reverse strand as negative MeH")
# Parse the command line arguments
args <- parse_args(p)

CG <- read.csv(args$Meh,header=TRUE)
CG=CG[which(apply(CG,1,function(x) sum(is.na(x)))==0),]

if( args$reverse == "all") {
  for (i in 1:dim(CG)[2]){
  if (!colnames(CG)[i] %in% c("chrom","bin","strand")){
    write.table(x = cbind(CG$chrom,format(CG$bin+(CG$strand=="r"), scientific = FALSE),
                          format(CG$bin+1+(CG$strand=="r"), scientific = FALSE),CG[,i]), 
                file= gsub(" ","",paste("PW_",colnames(CG)[i],".bedGraph")), row.names = FALSE, sep = " ",col.names = FALSE,quote = FALSE)
  }}       
} if (args$reverse == "n") {
  for (i in 1:dim(CG)[2]){
  if (!colnames(CG)[i] %in% c("chrom","bin","strand")){
  write.table(x = cbind(CG$chrom,format(CG$bin, scientific = FALSE),
                      format(CG$bin+1, scientific = FALSE),CG[,i]*(2*(CG$strand=="f")-1)), 
            file= gsub(" ","",paste("PW_",colnames(CG)[i],".bedGraph")), row.names = FALSE, sep = " ",col.names = FALSE,quote = FALSE)
  }}
} 