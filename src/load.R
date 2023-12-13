#' Available RAM in kB
CheckRAM <- function() as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))


#' Available cores
AvailableCores <- function(prop2use = .9) {
    max(1, floor(parallel::detectCores(logical = FALSE)*prop2use))
}


#' Get number of cores to fit RAM needs
Cores4RAM <- function(need) max(1, min(AvailableCores(), floor(CheckRAM() / need)))


#' Accessors to colData
CData <- function(name, vr, list = NULL) {
    if (!is.null(list)) {
        pluck(list, name, 'colData', 'listData', vr)
    }else{
        pluck(name, 'colData', 'listData', vr)
        }
}


#' Accessors to rowData
GData <- function(name, vr, list = NULL) {
    if (!is.null(list)) {
        pluck(.GlobalEnv, list, name, 'rowRanges', 'elementMetadata', 'listData', vr)
    }else{
        pluck(.GlobalEnv, name, 'rowRanges', 'elementMetadata', 'listData', vr)
    }
}


#' Read 10x data
#'
#' Read 10x data into a SingleCellExperiment object
#'
#' @param path Directory containing Cell Ranger counts matrix
#' @param dataset Name of the dataset
#' @param verbose Whether to print progress messages?
#'
#' @details
#' Based on DropletUtils::read10xCounts and adapted for Cell Ranger 3 output
#'
#' @return SingleCellExperiment object containing 10x data
read10x <- function(sample, path = here::here("data/cellranger"), verbose = FALSE) {

    path <- file.path(path, sample)
    
    if (verbose) {message("Reading gene info...")}
    gene_info <- readr::read_tsv(file.path(path, "genes.tsv.gz"),
                                 col_names = c("ID", "Symbol"),
                                 col_types = readr::cols(
                                     .default = readr::col_character()
                                 ))
    gene_info <- S4Vectors::DataFrame(gene_info)
    gene_info$Name <- scater::uniquifyFeatureNames(gene_info$ID,
                                                   gene_info$Symbol)
    rownames(gene_info) <- gene_info$Name

    if (verbose) {message("Reading cell info...")}
    cell_names <- readr::read_lines(file.path(path, "barcodes.tsv.gz"))
    n_cells <- length(cell_names)
    cell_ids <- str_c(sample, 1:n_cells, sep = "-")

    cell_info <- S4Vectors::DataFrame(Cell      = cell_ids,
                                      Sample    = rep(sample, n_cells),
                                      Barcode   = stringr::str_sub(cell_names, end =  -3),
                                      Batch     = rep(stringr::str_sub(sample, end =  -4), n_cells),
                                      Stage     = rep(stringr::str_sub(sample, start = 4), n_cells),
                                      row.names = NULL)

    if (verbose) {message("Reading expression matrix...")}
    mat <- as(Matrix::readMM(file.path(path, "matrix.mtx.gz")), "dgCMatrix")
    colnames(mat) <- cell_info$Cell

    if (verbose) {message("Creating SingleCellExperiment...")}
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = mat),
        rowData = gene_info,
        colData = cell_info
    )
    
    if (verbose) {message("Mark barcodes filtered by CellRanger in SingleCellExperiment object...")}
    filt_barcodes <- readr::read_lines(file.path(path, "filtered_barcodes.tsv.gz")) %>% 
        stringr::str_sub(string = ., end = -3)
    colData(sce)$CellRangerFilt <- colData(sce)$Barcode %in% filt_barcodes

    if (verbose) {message("Done!")}
    invisible(gc())
    return(sce)
}



#' Read dropEst
#'
#' Load data from dropEst corrected with Bayesian estimation
#' 
readDropEst <- function(sample, path = here::here("data/dropest"), verbose = FALSE) {
    holder <- read_rds(file.path(sprintf("%s/%s.rds", path, sample)))
    if (length(holder$reads_per_umi_per_cell$reads_per_umi[[1]][[1]]) != 2)
        stop("Quality must be provided")
    
    verbosity <- ifelse(verbose, 2, 0) 
    
    # Derive matrix corrected by Bayesian approach
    do_correction <- !file.exists(file.path(sprintf("%s/%s_corrected_Bayesian_cms.rds", path, sample)))
    if (do_correction) {
        reads_per_umi_per_cell <- holder$reads_per_umi_per_cell
        invisible(gc())
        
        umi_distribution <- dropestr::GetUmisDistribution(reads_per_umi_per_cell$reads_per_umi)
        umi_probs <- umi_distribution / sum(umi_distribution)
        collisions_info <- dropestr::FillCollisionsAdjustmentInfo(umi_probs, max(holder$cm))
        
        # Run Bayesian approach
        do_correction <- !file.exists(file.path(sprintf("%s/%s_corrected_Bayesian.rds", path, sample)))
        if (do_correction) {
            corrected_reads <- list()
            corrected_reads$Bayesian <- dropestr::CorrectUmiSequenceErrors(
                reads.per.umi.per.cb.info = reads_per_umi_per_cell,
                method = 'Bayesian',
                return = 'reads',
                collisions.info = collisions_info,
                umi.probabilities = umi_probs,
                verbosity.level = verbosity,
                mc.cores = if (exists("bpparam")) bpparam else Cores4RAM(27777777)
            )
            write_rds(corrected_reads, 
                      file.path(sprintf("%s/%s_corrected_Bayesian.rds", path, sample)))
        }else{
            corrected_reads <- read_rds(file.path(sprintf("%s/%s_corrected_Bayesian.rds", path, sample)))
        }
        
        BuildCountMatrixFromReads <- function(filt.rpus, reads.per.umi.per.cb.info, collisions.info=NULL) {
            filt.umis.per.gene <- sapply(filt.rpus, length)
            if (!is.null(collisions.info)) {
                filt.umis.per.gene <- collisions.info[filt.umis.per.gene]
            }
            reads.per.umi.per.cb.info$umis_per_gene <- filt.umis.per.gene
            return(dropestr::BuildCountMatrix(reads.per.umi.per.cb.info))
        }
        
        # Corrected Matrix
        corrected_cms <- lapply(
            corrected_reads,
            BuildCountMatrixFromReads,
            reads.per.umi.per.cb.info = holder$reads_per_umi_per_cell,
            collisions.info = collisions_info
        )
        
        corrected_cms <- lapply(corrected_cms,
                                function(cm) cm[grep("^[^;]+$", rownames(cm)), ])
        if (length(corrected_cms) == 1) {
            corrected_cms <- unlist(corrected_cms)
        }
        
        write_rds(corrected_cms, 
                  file.path(sprintf("%s/%s_corrected_Bayesian_cms.rds", path, sample)))
    }else{
        corrected_cms <- read_rds(file.path(sprintf("%s/%s_corrected_Bayesian_cms.rds", path, sample)))
    }
    
    if (typeof(corrected_cms) == "list") {
        if (length(corrected_cms) == 1) {
            corrected_cms <- unlist(corrected_cms)
        } else{
            corrected_cms <- corrected_cms[["Bayesian"]]
        }
    }
    
    holder$cm <- corrected_cms[grep("^[^;]+$", rownames(corrected_cms)),]
    holder$cm_raw <- holder$cm_raw[grep("^[^;]+$", rownames(holder$cm_raw)),]
    invisible(gc())
    return(holder)
}


prepRefGraph <- function(tax, iversion, outdir, biomart) {
    # Functions:
    ## Need to modify this function to work with version 17:
    get_irefindex <- function(tax_id = "All", iref_version = "current", data_folder = getwd()) {
    if (data_folder == "data") {
        datafolder = system.file("data", package = "iRefR")
    } else if (data_folder == "home") {
        datafolder = R.home()
    } else {
        datafolder = data_folder
    }
    if (iref_version == "current") {
        iref_version = "17.0"
        release_date = "27062020"
        url = paste("https://irefindex.vib.be/download/irefindex/data/archive/release_17.0/psi_mitab/MITAB2.6/",
                    tax_id, ".mitab.", release_date, ".txt.zip", sep = "")
    } else {
        if (iref_version == "17.0") {
        release_date = "27062020"
        url = paste("https://irefindex.vib.be/download/irefindex/data/archive/release_",
                    iref_version, "/psi_mitab/MITAB2.6/",
                    tax_id, ".mitab.", release_date, ".txt.zip", sep = "")
        }
    }
    file_location = paste(datafolder, "/", tax_id, ".mitab.",
                            release_date, ".txt", sep = "")
    if (file.exists(file_location) == TRUE) {
        cat("Reading available iRefIndex file...\n")
        irefindex_tab = unique(read.table(file_location, header = TRUE,
                                        comment.char = "", sep = "\t", quote = ""))
    } else {
        cat("Downloading iRefIndex file...\n")
        zipfile = paste(datafolder, "/", tax_id, ".mitab.", release_date,
                        ".txt.zip", sep = "")
        download.file(url, destfile = zipfile)
        unzip(zipfile, exdir = datafolder)
        file.remove(zipfile)
        cat("Reading downloaded file...\n")
        txtfile = paste(datafolder, "/", tax_id, ".mitab.", release_date,
                        ".txt", sep = "")
        irefindex_tab = unique(read.table(txtfile, header = TRUE,
                                        comment.char = "", sep = "\t", quote = ""))
        cat("File has been saved as:\n")
        cat(paste(txtfile, "\n"))
        save(file = paste(datafolder, "/", tax_id, ".mitab.",
                        release_date, ".RData", sep = ""), list = "irefindex_tab")
    }
    irefindex_tab
    }


    graph.stats <- function(g, suffix = "") {
        nodes <- V(g)
        n <- length(nodes)
        cat("Degrees...\n")
        x.deg <- igraph::degree(g)[names(nodes)]
        cat("Hub scores...\n")
        x.hub <- igraph::hub_score(g, scale = TRUE)$vector[names(nodes)]
        #cat("Authority scores...\n") #Same as hub scores here
        #x.aut <- igraph::authority_score(g, scale = TRUE)$vector[names(nodes)]
        cat("Betweenness...\n")
        x.bet <- igraph::betweenness(g)[names(nodes)]
        cat("Closeness...\n")
        x.clo <- igraph::closeness(g)[names(nodes)]
        cat("Transitivity (clustering coefficient)...\n")
        x.tra <- igraph::transitivity(g, type = "local", vids = nodes, isolates = "zero")
        df <- data.frame(degree = x.deg, hub_score = x.hub,
                        betweenness = x.bet, closeness = x.clo, transitivity = x.tra)
        names(df) <- paste(names(df), suffix, sep = "")
        df$icrogid <- names(nodes)
        return(df)
    }


    convert.ids <- function(stats, idtype = "uniprotkb", cname = "uniprotkb", id_conversion_table = id_conversion_table_mus) {
        cat(idtype, "IDs...\n")
        pb <- txtProgressBar(0, nrow(stats), style = 3)
        for (i in 1:nrow(stats)) {
            setTxtProgressBar(pb, i)
            uid <- convert_protein_ID("icrogid", stats[i, "icrogid"], idtype, id_conversion_table)
            stats[i, cname] <- ifelse(length(uid) == 1, uid, NA)
        }
        return(stats)
    }


    index  <- get_irefindex(tax_id = tax,
                            iref_version = iversion,
                            data_folder = outdir)

    # Get interaction types:
    index.binary  = select_interaction_type("binary", index)
    index.complex = select_interaction_type("complex", index)
    index.polymer = select_interaction_type("polymer", index)
    index.inter   = rbind(index.binary, index.complex)

    # Convert to edge lists:
    edgeList.binary  <- convert_MITAB_to_edgeList(index.binary)
    edgeList.complex <- convert_MITAB_to_edgeList(index.complex)
    edgeList.polymer <- convert_MITAB_to_edgeList(index.polymer)
    edgeList.inter   <- convert_MITAB_to_edgeList(index.inter)

    # Convert edge lists to graphs:
    graph.binary  <- convert_edgeList_to_graph(edgeList.binary)
    graph.complex <- convert_edgeList_to_graph(edgeList.complex)
    graph.polymer <- convert_edgeList_to_graph(edgeList.polymer)
    graph.inter   <- convert_edgeList_to_graph(edgeList.inter)

    # Degree distribution:
    plot(degree_distribution(graph.inter), log = "xy")

    # get all nodes and convert to uniprot:
    id_conversion_table_mus <- create_id_conversion_table(index,
                                                          tempdir(),
                                                          sprintf("id_convertable_%s", tax))

    stats <- graph.stats(graph.inter, ".ppi")
    stats <- convert.ids(stats, "uniprotkb", "uniprotkb")
    stats <- convert.ids(stats, "entrezgene/locuslink", "entrez")
    stats <- convert.ids(stats, "refseq", "refseq")
    cat("\nConversion was successful for:\n")
    cat(apply(stats, 2, function(x) sum(!is.na(x))))

    # Tell if genes are involved in complex interactions:
    cpx <- data.frame(id = names(V(graph.complex)), complex.ppi = TRUE)
    stats <- left_join(stats, cpx,
                        by = c("icrogid" = "id"))
    stats %<>% replace_na(list(complex.ppi = FALSE))
    cat("\nNumber of complex interactions:\n")
    cat(table(stats$complex.ppi))

    # Tell if genes are members of polymerisation processes
    pol <- data.frame(id = names(V(graph.polymer)), polymer.ppi = TRUE)
    stats <- left_join(stats, pol,
                        by = c("icrogid" = "id"))
    stats %<>% replace_na(list(polymer.ppi = FALSE))
    cat("\nNumber of polymer interactions:\n")
    cat(table(stats$polyner.ppi))

    # Add ensembl IDs and common gene names
    tln <- getBM(attributes = c("ensembl_gene_id",
                                "entrezgene_id",
                                "external_gene_name"),
                filters = "entrezgene_id",
                values = stats$entrez,
                mart = biomart,
                verbose = FALSE)
    tln$entrezgene_id <- as.character(tln$entrezgene_id)
    tln$external_gene_name <- as.character(tln$external_gene_name)
    tln %<>%
        distinct(ensembl_gene_id, .keep_all = TRUE) %>%
        distinct(entrezgene_id, .keep_all = TRUE) %>%
        distinct(external_gene_name, .keep_all = TRUE)
    stats <- left_join(stats, tln,
                       by = c("entrez" = "entrezgene_id"))
    skimr::skim(stats)

    path_stats <- here(outdir,
                       sprintf("ppi_centrality_statistics_%s_v%s.tsv",
                               tax, iversion))
    write_tsv(stats,
              file = path_stats)
    cat(sprintf("\nPPI Centrality statistics has been written to %s\n",
                    path_stats))

    path_graph <- here(outdir,
                       sprintf("ppi_graph_%s_iref_v%s.rds",
                               tax, iversion))
    write_rds(graph.inter,
              file = path_graph)
    cat(sprintf("PPI graph for taxon %s has been written to %s",
                    tax, path_graph))

    return(list(graph.inter, stats))
}
