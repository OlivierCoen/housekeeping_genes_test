FROM rocker/verse

RUN Rscript -e "install.packages('BiocManager', repos='http://cran.us.r-project.org', dependencies=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install(c('ArrayExpress', 'affy', 'edgeR', 'DESeq2', 'limma'), dependencies=TRUE, ask=FALSE)"