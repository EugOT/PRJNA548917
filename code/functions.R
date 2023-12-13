#' Get download markdown link
#'
#' Convert an output file name and location to a URL that can be used to
#' download the file.
#'
#' @param file name of the output file
#' @param folder name of the directory in the output directory containing the
#' output file
#'
#' @return Markdown URL link to the file
get_download_link <- function(file, folder = NULL) {
  remote <- workflowr::wflow_git_remote(verbose = FALSE)["origin"]

  url <- gsub(":", "/", remote)
  url <- gsub("git@", "http://", url)
  url <- gsub(".git", "", url, fixed = TRUE)
  url <- paste(url, "raw/master/output", sep = "/")

  if (is.null(folder)) {
    url <- paste(url, file, sep = "/")
  } else {
    url <- paste(url, folder, file, sep = "/")
  }

  link <- glue::glue("[{file}]({url})")

  return(link)
}


#' Helper function to write number of cells per file
#'
#' From the SRA id of dataset check number of cells in filtered seurat object
#'
#' @param sra id of dataset
#'
#' @return list with number of cells parameter, value and description
n_cells_per_file <- function(sra, srt = combined_srt) {
  dataset <- list(
    Parameter = sprintf("n_cells-%s", sra),
    Value = sum(srt$orig.ident == sra),
    Description = sprintf("Number of cells in the filtered %s dataset", sra)
  )
  return(dataset)
}

#' Write gene table
#'
#' Save a gene table to a file
#'
#' @param gene_table gene table to save
#' @param path file path to save location
#'
#' @details
#' data.frame objects will be saved as a (zipped) CSV file. List objects will
#' be saved in XLSX format.
write_gene_table <- function(gene_table, path) {
  if (is.data.frame(gene_table)) {
    zip_path <- paste0(path, ".zip")
    if (file.exists(zip_path)) {
      file.remove(zip_path)
    }
    readr::write_csv(gene_table, path, na = "")
    zip(zip_path, path, flags = "-q -j")
    invisible(file.remove(path))
  } else {
    writexl::write_xlsx(gene_table, path)
  }
}

#' Save plot
save_my_plot <- function(name,
                         plt,
                         type,
                         h = 12,
                         asp = 1.618,
                         path = plots_dir,
                         format = ".pdf") {
  cowplot::save_plot(
    filename = here::here(
      path,
      stringr::str_glue(type,
        as.character(name),
        format,
        .sep = "_"
      )
    ),
    plot = plt,
    base_height = h,
    base_asp = asp,
    limitsize = FALSE
  )
}


#' Read in Scrublet results
read_scrublet <- function(org, fpr) {
  read_tsv(here(
    "scrublet", org,
    sprintf("%s_scrublet_calls_FPR_%s.tsv", org, fpr)
  )) %>%
    mutate(
      origin = org,
      cell_name = str_c(org, barcode, sep = "_")
    )
}


#' Select mrtree resolution
select_resolution <- function(df) {
  #' Get distance between two resolutions with top ARI score
  get_top_res_diff <- function(dat) {
    tmp.ari <-
      dat |>
      top_n(n = 2, wt = ari) |>
      purrr::pluck(2)
    tmp.res <-
      dat |>
      top_n(n = 2, wt = ari) |>
      purrr::pluck(1)
    tmp.ari <- tmp.ari[1] - tmp.ari[2]
    tmp.res <- tmp.res[1] - tmp.res[2]
    return(c(tmp.ari, tmp.res))
  }

  #' Pick one of two top resolutions with parameters
  pick_res_param <- function(dat,
                             ari.dif, res.dif,
                             ari.thd = .05, res.thd = 0) {
    if (ari.dif < ari.thd & res.dif < res.thd) {
      res <-
        dat |>
        top_n(n = 2, wt = ari) |>
        purrr::pluck(1)
      res <- res[2]
    } else {
      res <-
        dat |>
        top_n(n = 1, wt = ari) |>
        purrr::pluck(1)
    }
    return(res)
  }


  df %<>% as_tibble()
  ein.check <-
    df |>
    top_n(n = 2, wt = ari) |>
    purrr::pluck(2) |>
    purrr::map_lgl(~ .x == 1)
  if (all(ein.check)) {
    df %<>% arrange(-resolution) |> distinct(ari, .keep_all = TRUE)
    c(tmp.ari, tmp.res) %<-% get_top_res_diff(df)
    resK <- pick_res_param(df, ari.dif = tmp.ari, res.dif = tmp.res)
  } else {
    df %<>%
      filter(ari != 1) |>
      arrange(-resolution) |>
      distinct(ari, .keep_all = TRUE)
    c(tmp.ari, tmp.res) %<-% get_top_res_diff(df)
    resK <- pick_res_param(df, ari.dif = tmp.ari, res.dif = tmp.res)
  }
  return(resK)
}


pc_score <- function(object = srt, PCs = 1:5, score.thresh = 1e-05) {
  pAll <- object[["pca"]]@jackstraw$empirical.p.values
  pAll <- pAll[, PCs, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- reshape2::melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  score.df <- NULL
  for (i in PCs) {
    pc.score <-
      suppressWarnings(
        prop.test(
          x = c(
            length(
              x = which(x = pAll[
                ,
                i
              ] <= score.thresh)
            ),
            floor(x = nrow(x = pAll) *
              score.thresh)
          ), n = c(nrow(pAll), nrow(pAll))
        )$p.val
      )
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    } else {
      score.df <- rbind(score.df, data.frame(PC = paste0(
        "PC",
        i
      ), Score = pc.score))
    }
  }
  return(score.df)
}


#' Derive MRTree clustering of Seurat object
derive_k_tree <- function(
    srt,
    n.pcs = n_pcs,
    vseed = reseed,
    n.cores = n_cores) {
  srt <- NormalizeData(srt)
  srt <-
    FindVariableFeatures(
      srt,
      selection.method = "vst",
      nfeatures = 3000
    )
  all.genes <- rownames(srt)
  hvg <- VariableFeatures(srt)
  var_regex <- "^Hla-|^Ig[hjkl]|^Rna|^mt-|^Rp[sl]|^Hb[^(p)]|^Gm"
  hvg <- hvg[str_detect(pattern = var_regex, string = hvg, negate = TRUE)]
  srt <-
    ScaleData(srt,
      features = all.genes,
      vars.to.regress = c(
        "var_regex", "log10GenesPerUMI",
        "S.Score", "G2M.Score"
      )
    )
  srt <-
    RunPCA(srt,
      features = hvg,
      npcs = n.pcs,
      seed.use = vseed,
      verbose = TRUE
    )
  srt <-
    JackStraw(
      object = srt,
      assay = "RNA",
      reduction = "pca",
      dims = n.pcs,
      num.replicate = 100,
      prop.freq = 0.01,
      maxit = 1000
    )
  srt <-
    ScoreJackStraw(srt,
      dims = seq_along(srt[["pca"]]@stdev)
    )
  test_pc <-
    pc_score(
      object = srt,
      PCs = seq_along(srt[["pca"]]@stdev),
      score.thresh = 1e-05
    )
  selected_pcs <-
    seq_along(
      srt[["pca"]]@stdev
    )[test_pc$Score <= 1e-03 &
      srt[["pca"]]@stdev > quantile(srt[["pca"]]@stdev, .25)]
  srt <-
    srt |>
    FindNeighbors(
      dims = selected_pcs,
      k.param = 15,
      annoy.metric = "euclidean",
      n.trees = 100,
      verbose = FALSE
    ) |>
    RunUMAP(
      dims = selected_pcs,
      reduction.name = "umap",
      reduction.key = "UMAP_",
      return.model = FALSE,
      umap.method = "umap-learn",
      densmap = TRUE,
      dens.lambda = 1L,
      dens.frac = 0.3,
      n.epochs = 1000L,
      n.neighbors = 15L,
      min.dist = 0.01,
      spread = 2L,
      metric = "correlation",
      init = "pca",
      seed.use = vseed,
      verbose = FALSE
    )
  resolutions <-
    modularity_event_sampling(
      A = srt@graphs$RNA_snn,
      n.res = 30,
      gamma.min = 0.2,
      gamma.max = 2.50001
    ) # sample based on the similarity matrix
  srt <- FindClusters(
    srt,
    algorithm = 4, method = "igraph",
    resolution = resolutions, random.seed = vseed,
    verbose = FALSE
  )
  out <- mrtree(
    srt,
    prefix = "RNA_snn_res.",
    n.cores = n.cores,
    consensus = FALSE,
    sample.weighted = TRUE,
    augment.path = FALSE,
    verbose = FALSE
  )
  # Adjusted Multiresolution Rand Index (AMRI)
  ks.flat <- apply(
    out$labelmat.flat,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  ks.mrtree <- apply(
    out$labelmat.mrtree,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  amri.flat <- sapply(seq_len(ncol(out$labelmat.flat)), function(i) {
    AMRI(out$labelmat.flat[, i], srt$seurat_clusters)$amri
  })
  amri.flat <- aggregate(amri.flat, by = list(k = ks.flat), FUN = mean)
  amri.recon <- sapply(seq_len(ncol(out$labelmat.mrtree)), function(i) {
    AMRI(out$labelmat.mrtree[, i], srt$seurat_clusters)$amri
  })
  df <- rbind(
    data.frame(
      k = amri.flat$k,
      amri = amri.flat$x,
      method = "Seurat flat"
    ),
    data.frame(k = ks.mrtree, amri = amri.recon, method = "MRtree")
  )
  stab.out <- stability_plot(out)
  resK <- SelectResolution(stab.out$df)
  srt$k_tree <- out$labelmat.mrtree[, which.min(
    abs(as.integer(
      str_remove(dimnames(
        out$labelmat.mrtree
      )[[2]], "K")
    ) - resK)
  )] %>%
    as.numeric() %>%
    as.factor()
  Idents(srt) <- "k_tree"
  if (length(unique(srt$k_tree)) > 1) {
    srt.markers.lr <-
      FindAllMarkers(
        srt,
        assay = "RNA",
        verbose = FALSE,
        random.seed = reseed,
        latent.vars = c(
          "var_regex", "log10GenesPerUMI",
          "S.Score", "G2M.Score"
        ),
        only.pos = TRUE,
        min.pct = 0.1,
        base = 10,
        logfc.threshold = 0.2,
        test.use = "LR"
      )
    if (length(unique(srt.markers.lr$cluster)) > 1) {
      write_csv(
        srt.markers.lr,
        here(
          tables_dir,
          sprintf(
            "%s_all-mrk_logreg.csv",
            unique(srt$orig.ident)
          )
        )
      )
    }
    srt.markers.mast <-
      FindAllMarkers(
        srt,
        assay = "RNA",
        verbose = FALSE,
        random.seed = reseed,
        latent.vars = c(
          "var_regex", "log10GenesPerUMI",
          "S.Score", "G2M.Score"
        ),
        only.pos = TRUE,
        min.pct = 0.1,
        base = 10,
        logfc.threshold = 0.2,
        test.use = "MAST"
      )
    if (length(unique(srt.markers.lr$cluster)) > 1) {
      write_csv(
        srt.markers.lr,
        here(
          tables_dir,
          sprintf(
            "%s_all-mrk_mast.csv",
            unique(srt$orig.ident)
          )
        )
      )
    }
  }
  return(list(srt, srt.markers.lr, srt.markers.mast))
}

###### COLOURS: ################################################################
qc_palette <- c(
  "Doublet" = "#5050FFFF",
  "Pass" = "#ff3700f1",
  "High_MT" = "#6BD76BFF",
  "High_Hgb" = "#C75127FF",
  "Low_Complexity" = "#3B1B53FF",
  "High_Ribo" = "#E7C76FFF",
  "Low_nFeature" = "#A9A9A9FF",
  "High_UMIs" = "#1A0099FF",
  "High_Hgb,High_Ribo,Low_nFeature" = "#660099FF",
  "High_UMIs,Doublet" = "#990080FF",
  "High_Hgb,High_Ribo,High_MT,Low_nFeature,Doublet" = "#FF1463FF",
  "High_Ribo,Low_nFeature,Doublet" = "#00D68FFF",
  "High_Hgb,High_Ribo,Low_nFeature,Doublet" = "#D60047FF",
  "High_Ribo,High_MT" = "#00CC33FF",
  "High_Hgb,High_Ribo" = "#4775FFFF",
  "High_Ribo,High_MT,Low_Complexity" = "#991A00FF",
  "High_Hgb,Low_nFeature,Low_Complexity" = "#00991AFF",
  "High_Hgb,Low_nFeature" = "#CE3D32FF",
  "High_Hgb,High_MT" = "#BA6338FF",
  "High_MT,Low_nFeature" = "#D595A7FF",
  "High_Ribo,High_MT,Doublet" = "#D58F5CFF",
  "High_Hgb,High_Ribo,High_MT,Low_nFeature" = "#CDDEB7FF",
  "High_MT,Low_Complexity" = "#5A655EFF",
  "High_Ribo,High_MT,Low_nFeature" = "#CC9900FF",
  "High_MT,Low_nFeature,Low_Complexity" = "#00CC99FF",
  "High_MT,Doublet" = "#FFC20AFF",
  "High_Hgb,High_Ribo,High_MT" = "#996600FF",
  "Low_nFeature,Low_Complexity" = "#009966FF",
  "High_Ribo,Doublet" = "#749B58FF",
  "High_Hgb,Doublet" = "#5DB1DDFF",
  "Low_nFeature,Doublet" = "#924822FF",
  "High_Hgb,High_MT,Low_nFeature" = "#7A65A5FF",
  "High_Hgb,High_Ribo,Low_Complexity" = "#612A79FF",
  "High_MT,Low_nFeature,Doublet" = "#CC9900FF",
  "High_Ribo,High_MT,Low_nFeature,Low_Complexity" = "#99CC00FF",
  "High_UMIs,Low_Complexity" = "#0099CCFF",
  "High_Ribo,High_UMIs,Low_Complexity" = "#FFD147FF",
  "High_UMIs,Low_Complexity,Doublet" = "#809900FF",
  "High_Hgb,High_Ribo,High_MT,Doublet" = "#008099FF",
  "High_Hgb,Low_Complexity" = "#F0E685FF",
  "High_Ribo,Low_nFeature" = "#802268FF",
  "High_UMIs,High_MT,Low_Complexity" = "#837B8DFF",
  "High_Hgb,Low_Complexity,Doublet" = "#E4AF69FF",
  "High_Hgb,High_Ribo,Doublet" = "#AE1F63FF",
  "High_Hgb,High_MT,Low_nFeature,Low_Complexity" = "#99CC00FF",
  "High_Hgb,High_Ribo,High_MT,Low_Complexity,Doublet" = "#33CC00FF",
  "High_Hgb,High_Ribo,High_UMIs,Low_Complexity" = "#0A47FFFF",
  "Low_Complexity,Doublet" = "#990033FF",
  "High_Ribo,Low_Complexity" = "#339900FF",
  "High_Hgb,High_Ribo,Low_nFeature,Low_Complexity" = "#003399FF"
)
