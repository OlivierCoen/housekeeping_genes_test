library("ExpressionAtlas")

cwd <- getwd()
folder_name <- "output"

# setting outdir
outdir <- file.path(cwd, folder_name)

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

properties <- c('salt')
species <- 'arabidopsis'

atlasRes <- searchAtlasExperiments( properties = properties, species = species )
if (length(atlasRes) == 0) {
    properties_string <- paste(properties, collapse = ", ")
    msg <- paste('No Expression Atlas experiment found with keywords', properties_string)
    stop(msg)
}

atlasData <- getAtlasData( atlasRes$Accession )
#atlasData <- getAtlasData( "E-GEOD-10496" )

# looping through experiments fetched
for (experiment_id in names(atlasData)) {

  eset <- atlasData[[ experiment_id ]]

  # looping through each data type (ex: 'rnaseq') in the experiment
  for (data_type in names(eset)) {

    data <- eset[[ data_type ]]

    skip_iteration <- FALSE
    # getting count dataframe
    tryCatch({
        df <- assays(data)$counts
    }, error = function(e) {
        print(paste("Caught an error: ", e$message))
        print(paste('ERROR: Could not get assay data for experiment ID', experiment_id, 'and data type', data_type))
        skip_iteration <- TRUE
    })

    # If an error occurred, skip to the next iteration
    if (skip_iteration) {
        next
    }

    outfilename <- paste0(experiment_id, '_', data_type, '.csv')
    outfile <- file.path(outdir, outfilename)
    # exporting to CSV file
    # index represents gene names
    write.table(df, outfile, sep=',')
  }
  
}