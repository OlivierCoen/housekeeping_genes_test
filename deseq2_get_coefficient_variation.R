library(DESeq2)
library(argparser)
library(tibble)

get_args <- function() {
    argp <- arg_parser(description="")
    argp <- add_argument(
            parser=argp,
            arg="--count-file",
            help="file containing counts for all samples")
    argp <- add_argument(
            parser=argp,
            arg="--outfile",
            help="outfile")
    argv <- parse_args(
            parser=argp,
            argv=commandArgs(trailingOnly = TRUE))
    return(argv)
}

arg <- get_args()

data <- read.csv(arg$count_file, header=TRUE)
gene_ids <- rownames(data)

# we do not need an experiment design since we just want coefficient of variations
# making fake experiment design
colData <- DataFrame(rep("single", ncol(data)))

dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = ~ 1)

# estimating size factors for normalization and dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
gene_disp <- dispersions(dds)

# calculating the square root of the dispersion as a CV-like measure
cv <- sqrt(gene_disp)
# concatenate original columns with normalized counts

cv_df <- cbind(gene_ids, data.frame(cv))
# set gene ids as index

cv_df <- column_to_rownames(cv_df, var = "gene_ids")

write.table(cv_df, sep="\t", file=arg$outfile, row.names = TRUE, quote = FALSE)
