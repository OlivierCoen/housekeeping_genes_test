library(edgeR)
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
gene_counts <- data.matrix(data)

dge <- DGEList(counts=gene_counts, genes=gene_ids)

#normalization
dge <- calcNormFactors(dge, method="TMM")
norm_counts <- cpm(dge, lof=FALSE)

# estimating dispersion
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# getting the coefficient of variation as the sqrt of the tagwise dispersion
cv <- sqrt(dge$tagwise.dispersion)

#design.mat <- model.matrix(~ 0 + d$samples$group)
#colnames(design.mat) <- levels(d$samples$group)
# dge <- estimateGLMCommonDisp(dge,design.mat)
# dge <- estimateGLMTrendedDisp(dge,design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
# dge <- estimateGLMTagwiseDisp(dge,design.mat)

# concatenate original columns with normalized counts
df <- cbind(gene_ids, data.frame(cv))
# set gene ids as index
df <- column_to_rownames(df, var = "gene_ids")
print(head(df))

write.table(df, sep="\t", file=arg$outfile, row.names = TRUE, quote = FALSE)