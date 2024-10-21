#' @title run init analysis of scrublet
#' @description See preprint: Scrublet: computational identification of cell doublets in single-cell transcriptomic data
#' Samuel L Wolock, Romain Lopez, Allon M Klein.  bioRxiv 357368; doi: https://doi.org/10.1101/357368
#' @param seurat_obj the `Seurat` object to perform Scrublet
#' @param python_home The python home directory where Scrublet is installed
#' @param min_counts, int (optional, default=2), See scrublet reference
#' @param min_cells, int (optional, default=3), See scrublet reference
#' @param expected_doublet_rate, float (optional, default=0.06), See scrublet reference - expected_doublet_rate: the
#' fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter. For this example, the expected doublet rate comes from the Chromium User Guide: https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC
#' @param min_gene_variability_pctl, int (optional, default=85), See scrublet reference
#' @param n_prin_comps, int (optional, default=50), See scrublet reference  (Number of principal components to use)
#' @param sim_doublet_ratio, int (optional, default=2),  the number of doublets to simulate, relative to the number of observed transcriptomes. This should be high enough that all doublet states are well-represented by simulated doublets. Setting it too high is computationally expensive. The default value is 2, though values as low as 0.5 give very similar results for the datasets that have been tested.
#' @param n_neighbors, int (optional) n_neighbors: Number of neighbors used to construct the KNN classifier of observed transcriptomes and simulated doublets. The default value of round(0.5*sqrt(n_cells)) generally works well.
#' Return only a list containing scrublet output
#' @return return a python scrublet object, it will be used for further analysis.
#'
#' @export
#'
get_init_scrublet <- function(seurat_obj, python_home = Sys.which("python"),
                              min_counts=2, min_cells=3,
                              expected_doublet_rate=0.06,
                              min_gene_variability_pctl=85,
                              n_prin_comps=50, sim_doublet_ratio=2,
                              n_neighbors=NULL) {
  ## source py
  source_py(python_home = python_home,
            py_script = paste(system.file(package = "scrubletR"), "scrublet.py", sep = "/"))

  ## prepare the data
  count_mat <- Seurat::GetAssayData(seurat_obj, layer = "counts", assay = 'RNA')
  X <- as(Matrix::t(count_mat), "TsparseMatrix")
  i <- as.integer(X@i)
  j <- as.integer(X@j)
  val <- X@x
  dim <- as.integer(X@Dim)
  if(is.null(n_neighbors)){
    n_neighbors<-round(0.5*sqrt(nrow(X)))
  }

  ## do the scrublet analysis
  scrublet_py_args<-c(list(i=i, j=j, val=val, dim=dim,
                           expected_doublet_rate=expected_doublet_rate, min_counts=min_counts,
                           min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=n_prin_comps,
                           sim_doublet_ratio=sim_doublet_ratio, n_neighbors=n_neighbors))
  scrublet_res <- do.call(scrublet_py, scrublet_py_args)

  return(scrublet_res)
}



#' @title plot the doublet score histogram
#' @description plot the histogram of doublet scores for observed transcriptomes and simulated doublets
#' @param scrublet_obj The result from \code{get_init_scrublet} or \code{call_doublets}
#' @param threshold the threshold to select doublet. If NULL, use the threshold in \code{scrublet_obj}
#'
#' @import ggplot2
#' @importFrom scales trans_breaks trans_format math_format
#' @importFrom patchwork wrap_plots
#'
#' @return return a `ggplot2` object
#'
#' @export
#'
plot_histogram <- function(scrublet_obj, threshold = NULL) {

  if (is.null(threshold)) {
    thres <- scrublet_obj$threshold_
  } else {
    thres <- threshold
  }

  w <- as.data.frame(scrublet_obj$doublet_scores_obs_)
  names(w) <- 'x'

  p1 <- ggplot(w) +
    aes(x) +
    geom_histogram(bins = 50) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(limits = c(0, 1)) +
    annotation_logticks(sides="l") +
    geom_vline(xintercept = thres, linetype = 'solid', color = 'black') +
    labs(x = 'Doublet score', y = 'Prob. density', title = 'Observed transriptomes') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = 'black'))

  k = as.data.frame(scrublet_obj$doublet_scores_sim_)
  names(k) <- 'x'

  p2 <- ggplot(k) +
    aes(x) +
    geom_histogram(bins = 50) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_continuous(limits = c(0, 1)) +
    annotation_logticks(sides="l") +
    geom_vline(xintercept = thres, linetype = 'solid', color = 'black') +
    labs(x = 'Doublet score', y = 'Prob. density', title = 'Observed transriptomes') +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(color = 'black'))

  p_list <- patchwork::wrap_plots(list(p1, p2), ncol = 2)

  return(p_list)
}


#' @title call doublets by new threshold
#' @description Call trancriptomes as doublets or singlets
#' @param scrublet_obj The result from \code{get_init_scrublet}
#' @param threshold Doublet score threshold for calling a transcriptome a doublet.
#'        If `None`, this is set automatically.
#' @param return_obj If FALSE, only return a logical vector show doublet, or a python scrublet object.
#' @param python_home The python home directory where Scrublet is installed
#'
#' @return return a logical vector or a python scrublet object.
#'
#' @export
#'
call_doublets <- function(scrublet_obj, threshold = NULL, python_home = Sys.which("python"), return_obj = TRUE) {
  if (is.null(threshold)) {
    thres <- reticulate::py_none()
  } else {
    thres <- threshold
  }

  ## source py
  source_py(python_home = python_home,
            py_script = paste(system.file(package = "scrubletR"), "scrublet.py", sep = "/"))

  scrub_obj_new = call_scrublet_doublets(scrublet_obj, threshold = thres)

  if (return_obj) {
    return(scrub_obj_new)
  } else {
    scrub_obj_new$predicted_doublets_
  }
}


#' @title scrublet_R
#' @description See preprint: Scrublet: computational identification of cell doublets in single-cell transcriptomic data
#' Samuel L Wolock, Romain Lopez, Allon M Klein.  bioRxiv 357368; doi: https://doi.org/10.1101/357368
#' @param seurat_obj the `Seurat` object to perform Scrublet
#' @param python_home The python home directory where Scrublet is installed
#' @param return_results_only bool (optional, default False)
#' @param min_counts, int (optional, default=2), See scrublet reference
#' @param min_cells, int (optional, default=3), See scrublet reference
#' @param expected_doublet_rate, float (optional, default=0.06), See scrublet reference - expected_doublet_rate: the
#' fraction of transcriptomes that are doublets, typically 0.05-0.1. Results are not particularly sensitive to this parameter. For this example, the expected doublet rate comes from the Chromium User Guide: https://support.10xgenomics.com/permalink/3vzDu3zQjY0o2AqkkkI4CC
#' @param min_gene_variability_pctl, int (optional, default=85), See scrublet reference
#' @param n_prin_comps, int (optional, default=50), See scrublet reference  (Number of principal components to use)
#' @param sim_doublet_ratio, int (optional, default=2),  the number of doublets to simulate, relative to the number of observed transcriptomes. This should be high enough that all doublet states are well-represented by simulated doublets. Setting it too high is computationally expensive. The default value is 2, though values as low as 0.5 give very similar results for the datasets that have been tested.
#' @param n_neighbors, int (optional) n_neighbors: Number of neighbors used to construct the KNN classifier of observed transcriptomes and simulated doublets. The default value of round(0.5*sqrt(n_cells)) generally works well.
#' Return only a list containing scrublet output
#' @param threshold Doublet score threshold for calling a transcriptome a doublet. if set NULL, set this threshold automatically.
#' @return return a `Seurat` object with an additional column added to `meta.data` with both the doublet_score output from scrublet, and predicted_doublets
#' @importFrom reticulate use_python source_python py_module_available
#' @importFrom Matrix t
#' @importFrom methods as
#' @import Seurat
#' @export
scrublet_R <- function(seurat_obj, python_home = Sys.which("python"),
                       return_results_only = FALSE, min_counts=2,
                       min_cells=3, expected_doublet_rate=0.06,
                       min_gene_variability_pctl=85,
                       n_prin_comps=50, sim_doublet_ratio=2, n_neighbors=NULL,
                       threshold = NULL) {
  ## source py
  source_py(python_home = python_home,
            py_script = paste(system.file(package = "scrubletR"), "scrublet.py", sep = "/"))

  ## do the scrublet analysis
  scrublet_obj <- get_init_scrublet(seurat_obj = seurat_obj, python_home = python_home,
                                    min_counts=min_counts, min_cells=min_cells,
                                    expected_doublet_rate=expected_doublet_rate,
                                    min_gene_variability_pctl=min_gene_variability_pctl,
                                    n_prin_comps=n_prin_comps, sim_doublet_ratio=sim_doublet_ratio,
                                    n_neighbors=n_neighbors)
  if (!is.null(threshold)) {
    scrublet_obj <- call_doublets(call_doublets, threshold = threshold)
  }
  scrublet_res <- list(doublet_scores = scrublet_obj$doublet_scores_obs_,
                       predicted_doublets = scrublet_obj$predicted_doublets_)

  names(scrublet_res)<-c("doublet_scores", "predicted_doublets")

  ## get the result
  if (return_results_only) {
    return(scrublet_res)
  } else {
    seurat_obj[["doublet_scores"]] <- scrublet_res$doublet_scores
    seurat_obj[["predicted_doublets"]] <- scrublet_res$predicted_doublets
    return(seurat_obj)
  }
}


source_py <- function(python_home = Sys.which("python"),
                      py_script = paste(system.file(package = "scrubletR"), "scrublet.py", sep = "/")) {
  ## use the python
  reticulate::use_python(python_home)

  ## test whether have the modules
  if (!reticulate::py_module_available("scrublet")) {
    stop("python module scrublet does not seem to be installed; - try running 'py_config()'")
  }

  ## source .py file
  reticulate::source_python(py_script)

}
