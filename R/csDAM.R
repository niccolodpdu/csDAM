#### Compute cosine profile curves for each gene (with power transformation) ####
#' @title Compute cosine profile curves for each gene
#'
#' @description
#' Computes cumulative cosine profile curves for each gene, optionally applying
#' a mean-ratio power transformation to emphasize genes with few high and many
#' low expression values. The resulting matrix shows, for each gene, the cosine
#' similarity between its descending-sorted expression pattern and reference
#' vectors of the form (1, ..., 1, 0, ..., 0) for k = 1..n_samples.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item (Optional) Apply a mean-ratio power transform  
#'         \deqn{x'_i = (x_i / \mu)^\gamma}
#'         where \eqn{\mu} is the mean or median of the gene's values.  
#'         This transformation accentuates front-loaded profiles while
#'         suppressing long tails of small nonzero values.
#'   \item Sort each gene’s expression values in descending order.
#'   \item Compute cumulative cosine similarities with partial reference vectors
#'         of length k:
#'         \deqn{res[k] = \frac{\sum_{j \le k} x_{(j)}}
#'                        {\sqrt{k} \; \|x\|_2}}
#' }
#'
#' If \code{auto_gamma = TRUE}, the power exponent \eqn{\gamma} is automatically
#' scaled with sample size (number of columns) as
#' \deqn{\gamma_\mathrm{eff} = \gamma \times
#'       \frac{\log_2(10 + n_\mathrm{smp})}{\log_2(60)}}
#' which slightly increases nonlinearity for longer vectors to compensate for
#' dilution by many small values.
#'
#' @param data A numeric matrix with samples in columns and genes in rows.
#' @param transform Character string specifying the transformation to apply:
#'   \itemize{
#'     \item \code{"none"} – no transformation (default);
#'     \item \code{"mean_power"} – mean-ratio power transform.
#'   }
#' @param gamma Numeric scalar. Base power exponent controlling nonlinearity
#'   (recommended range 1.2–1.8). Ignored if \code{transform = "none"}.
#' @param auto_gamma Logical. If \code{TRUE} (default), \eqn{\gamma} is adjusted
#'   automatically with sample size.
#' @param use_median Logical. If \code{TRUE} (default), use the median instead
#'   of the mean as the central value \eqn{\mu} for robustness.
#'
#' @return
#' A numeric matrix of dimension \code{[n_genes × n_samples]}, where each row
#' represents one gene’s cumulative cosine profile curve.
#' Columns correspond to k = 1..n_samples (i.e., cumulative fraction of samples).
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rpois(200, lambda = 5), nrow = 50, ncol = 4)
#'
#' # No transformation
#' cos_mat0 <- cos_iDEG_fast(data, transform = "none")
#'
#' # Mean-ratio power transform with automatic gamma adjustment
#' cos_mat1 <- cos_iDEG_fast(data, transform = "mean_power")
#'
#' # Fixed gamma = 1.8, using mean instead of median
#' cos_mat2 <- cos_iDEG_fast(data, transform = "mean_power",
#'                            gamma = 1.8, auto_gamma = FALSE,
#'                            use_median = FALSE)
#'
#' @export
cos_iDEG_fast <- function(data,
                          transform = c("none", "mean_power"),
                          gamma = 1.5,
                          auto_gamma = TRUE,
                          use_median = TRUE) {
  
  transform <- match.arg(transform)
  n_genes <- nrow(data)
  n_smp   <- ncol(data)
  
  # ---------- Helper: mean-ratio power transform ----------
  mean_ratio_power <- function(v, gamma, use_median = TRUE) {
    mu <- if (use_median) median(v) else mean(v)
    (pmax(v, 0) / (mu + 1e-8))^gamma
  }
  
  # ---------- Optional transformation ----------
  if (transform == "mean_power") {
    # auto-adjust gamma based on sample size
    if (auto_gamma) {
      gamma_eff <- gamma * log2(10 + n_smp) / log2(60)
    } else {
      gamma_eff <- gamma
    }
    
    data <- t(apply(data, 1, mean_ratio_power,
                    gamma = gamma_eff, use_median = use_median))
  }
  
  # ---------- Core cosine accumulation ----------
  res <- matrix(NA_real_, nrow = n_genes, ncol = n_smp)
  
  for (i in 1:n_genes) {
    this_gene <- sort(data[i, ], decreasing = TRUE)
    norm_val  <- sqrt(sum(this_gene^2))
    if (norm_val == 0) {
      res[i, ] <- 0
      next
    }
    cumsums <- cumsum(this_gene)
    res[i, ] <- cumsums / (sqrt(1:n_smp) * norm_val)
  }
  
  rownames(res) <- rownames(data)
  colnames(res) <- paste0("k", 1:n_smp)
  return(res)
}


#### Select genes by best-k threshold from cosine profiles ####
#' Select genes by k threshold
#'
#' This function selects genes based on the maximum cosine similarity index 
#' and orders their expression across samples. 
#' 
#' Optionally it rules out features that are unlikely to be the markers, e.g., 
#' those features that resemble 1,1,...,1,1 than any other reference vectors.
#' 
#' @title Select features by k threshold
#' @param data A gene expression matrix (rows = genes, columns = samples).
#' @param cos_mat A cosine similarity matrix (rows = genes, columns = samples).
#' @param k_thr Threshold of k for selecting genes. If NULL (default), it will  
#'     be set to the number of samples (i.e., ncol(data)).
#' @return A list containing:
#'   \item{selected_genes}{Names of selected genes.}
#'   \item{best_k}{Vector of best-k indices for selected genes.}
#'   \item{pi}{A list of permutations (descending order of expression across samples).}
#' @examples
#' \dontrun{
#' res <- select_genes_by_k(data, cos_mat)
#' }
#' 
#' @export
select_genes_by_k <- function(data, cos_mat, k_thr = NULL) {
  stopifnot(is.matrix(data), is.matrix(cos_mat))
  stopifnot(nrow(data) == nrow(cos_mat))  # genes must match
  
  # Set default threshold to number of samples
  if (is.null(k_thr)) {
    k_thr <- ncol(data)
  }
  
  # 1. For each gene (row), find the best-k
  best_k <- apply(cos_mat, 1, which.max)
  
  # 2. Select genes with best_k < threshold
  sel_idx <- which(best_k < k_thr)
  
  # 3. Compute pi (descending order of expression across samples)
  pi_list <- lapply(sel_idx, function(g) order(data[g, ], decreasing = TRUE))
  names(pi_list) <- rownames(data)[sel_idx]
  
  # 4. Return results
  list(
    selected_genes = rownames(data)[sel_idx],
    best_k = best_k[sel_idx],
    pi = pi_list
  )
}


#### Group genes by their top-n π indices, with optional merging ####
#'
#' @title Group genes by π ranking
#' @description
#' This function groups genes based on the top-n indices of their π ranking.
#' Genes sharing the same top-n set (unordered) are grouped together. Groups
#' can optionally be merged if they overlap in their top-n indices by more
#' than a user-defined threshold. The order in which groups are processed can
#' be controlled (\code{"size"} vs \code{"alphabet"}). When merging with
#' \code{key_mode = "union"}, the key expansion can be limited with
#' \code{union_max}.
#'
#' @param cos_mat A numeric matrix of cosine values, with genes as columns
#'   and possible k-values as rows (output of \code{cos_iDEG_fast}).
#' @param best_k An integer vector giving the best k-value for each gene.
#'   Names of this vector must match the column names of \code{cos_mat}.
#' @param pi_list A list of integer vectors, one per gene, representing the
#'   ranking of samples (π). Names must correspond to genes.
#' @param top_n Integer, how many of the top π indices to use for initial grouping.
#' @param merge_k Integer, minimum number of overlapping indices (minus 1)
#'   required to merge groups. If set to 0, no merging is performed.
#' @param key_mode Character, either \code{"largest"} or \code{"union"}.
#'   \code{"largest"} keeps the key of the larger group when merging.
#'   \code{"union"} expands the key to the union of both groups (snowball logic).
#' @param init_order Character, either \code{"size"} or \code{"alphabet"}.
#'   Controls the initial order of groups: by group size (descending) or
#'   alphabetical order of group keys.
#' @param union_max Integer, maximum size of a key when using
#'   \code{key_mode = "union"}. Once this size is reached, further merges do
#'   not expand the key but still merge genes. Default is \code{Inf} (no limit).
#'
#' @return A list of groups, where each element is a character vector of genes
#'   belonging to that group. Each group is internally sorted by best_k
#'   (ascending) and then cosine value (descending). An attribute
#'   \code{"merged_history"} records details of merges performed, including
#'   the original keys and the final key.
#'
#' @examples
#' \dontrun{
#' sorted_groups <- group_genes_by_pi(cos_mat, best_k, pi_list,
#'                                    top_n = 3, merge_k = 1,
#'                                    key_mode = "union",
#'                                    init_order = "size",
#'                                    union_max = 5)
#' }
#'

#' @export
group_genes_by_pi <- function(cos_mat, best_k, pi_list,
                              top_n = 3, merge_k = 0,
                              key_mode = c("largest", "union"),
                              init_order = c("size", "alphabet"),
                              union_max = Inf) {
  key_mode <- match.arg(key_mode)
  init_order <- match.arg(init_order)
  merged_history <- list()
  
  genes <- names(pi_list)
  
  # 1. Build grouping key = first top_n of π (unordered set)
  key_list <- lapply(pi_list, function(pi) sort(pi[1:min(top_n, length(pi))]))
  key_str  <- vapply(key_list, paste, collapse = "-", FUN.VALUE = character(1))
  
  # 2. Split genes into initial groups
  groups <- split(genes, key_str)
  group_keys <- split(key_list, key_str)
  
  # Order groups
  if (init_order == "size") {
    group_sizes <- sapply(groups, length)
    groups <- groups[order(group_sizes, decreasing = TRUE)]
    group_keys <- group_keys[names(groups)]
  } else if (init_order == "alphabet") {
    groups <- groups[order(names(groups))]
    group_keys <- group_keys[names(groups)]
  }
  
  # 3. Post-hoc merging if merge_k > 0
  if (merge_k > 0 && length(groups) > 1) {
    changed <- TRUE
    while (changed) {
      changed <- FALSE
      new_groups <- list()
      new_keys   <- list()
      used <- rep(FALSE, length(groups))
      
      for (i in seq_along(groups)) {
        if (used[i]) next
        current_genes <- groups[[i]]
        current_key   <- group_keys[[i]][[1]]
        merged_keys   <- names(groups)[i]
        
        for (j in seq_along(groups)) {
          if (i == j || used[j]) next
          overlap <- length(intersect(group_keys[[i]][[1]], group_keys[[j]][[1]]))
          if (overlap >= merge_k + 1) {
            # merge genes
            current_genes <- c(current_genes, groups[[j]])
            
            if (key_mode == "union") {
              # expand only if current length < union_max
              new_union <- union(current_key, group_keys[[j]][[1]])
              if (length(new_union) <= union_max) {
                current_key <- new_union
              }
              # else: keep current_key unchanged (only merge genes)
            } else if (key_mode == "largest") {
              if (length(groups[[j]]) > length(groups[[i]])) {
                current_key <- group_keys[[j]][[1]]
              }
            }
            
            merged_keys <- c(merged_keys, names(groups)[j])
            used[j] <- TRUE
            changed <- TRUE
          }
        }
        
        new_name <- paste(sort(current_key), collapse = "-")
        new_groups[[new_name]] <- unique(current_genes)
        new_keys[[new_name]]   <- list(sort(current_key))
        merged_history[[new_name]] <- list(
          original_keys = unique(merged_keys),
          final_key = new_name
        )
        used[i] <- TRUE
      }
      groups <- new_groups
      group_keys <- new_keys
    }
  }
  
  # 4. Sort genes inside each group (by best_k, then cosine)
  sorted_groups <- lapply(groups, function(glist) {
    ks <- as.numeric(best_k[glist])
    cos_vals <- vapply(glist, function(g) {
      k <- as.integer(best_k[g])
      cos_mat[g, k]
    }, numeric(1))
    
    ord <- order(ks, -cos_vals)
    glist[ord]
  })
  
  attr(sorted_groups, "merged_history") <- merged_history
  return(sorted_groups)
}




#### Select representative sources from grouped genes ####
#'
#' @title Select representative sources
#'
#' @description
#' Selects representative sources (groups of genes) from the output of
#' \code{group_genes_by_pi()} using either distance-based or size-based
#' strategies.  
#' 
#' Depending on the chosen \code{method}, sources can be prioritized by
#' centroid dissimilarity (cosine or Spearman distance) or by the number of
#' genes contained.  
#' 
#' - When \code{method = "cosine"} (default), group centroids are compared
#'   by cosine distance, and the most distinct groups are iteratively
#'   selected.  
#' - When \code{method = "spearman"}, the same logic is applied using
#'   Spearman rank correlation instead of cosine similarity.  
#' - When \code{method = "size_only"}, all distance-based calculations are
#'   skipped and only the groups with the largest gene counts are returned.
#'
#' The default configuration (\code{method = "cosine"}, \code{mode = "distance"},
#' \code{cutoff = 0.2}, \code{first = "maxdist"}) iteratively adds sources
#' that are at least 0.2 distance apart from all previously selected ones,
#' starting from the most isolated centroid.
#'
#' @param groups A list of groups (e.g., output of \code{group_genes_by_pi}),
#'   where each element is a character vector of gene names.
#' @param data A numeric matrix with samples in columns and genes in rows.
#'   Used to compute group centroids if \code{method = "cosine"} or
#'   \code{method = "spearman"}.
#' @param method Character, one of \code{"cosine"}, \code{"spearman"},
#'   or \code{"size_only"}.
#'   - \code{"cosine"} (default): select sources based on centroid cosine distances.  
#'   - \code{"spearman"}: select sources based on centroid Spearman correlation.  
#'   - \code{"size_only"}: skip distance-based logic, return groups with the
#'     largest number of markers (size). Ignores \code{mode}, \code{cutoff},
#'     and \code{first}.
#' @param mode Character, one of \code{"distance"} or \code{"number"}.
#'   Used only when \code{method != "size_only"}.
#'   - \code{"distance"} (default): iteratively select sources whose minimum
#'     distance to already selected ones exceeds the threshold in
#'     \code{cutoff}.  
#'   - \code{"number"}: select exactly \code{n_sources} sources by maximizing
#'     centroid separation.
#' @param cutoff Numeric, distance threshold used when
#'   \code{mode = "distance"} (default = 0.2). 
#' @param n_sources Integer, number of sources to return when
#'   \code{mode = "number"} or \code{method = "size_only"} (default = 5).
#' @param per Integer, optional (default = 1). If not NULL, limit each source to the first
#'   \code{per} genes (according to their internal sorting).
#' @param centroid_genes Integer, optional. Number of top genes to use per
#'   source when computing centroids (default = use all genes). Ignored when
#'   \code{method = "size_only"}.
#' @param first Character, either \code{"maxdist"} or \code{"maxgenes"}.
#'   Determines how the first source is selected when \code{mode = "number"}
#'   and \code{method != "size_only"}. Default = \code{"maxdist"} (the most
#'   isolated centroid).
#'
#' @return
#' A list of selected sources (subset of \code{groups}).  
#' Each source is represented as a character vector of genes.  
#' If \code{per} is not NULL, each group is truncated accordingly.
#'
#' @examples
#' \dontrun{
#' # --- Default behavior (cosine distance, distance mode) ---
#' sel <- select_sources(groups, df_raw)
#'
#' # --- Select top 10 sources by maximal centroid distance ---
#' sel <- select_sources(groups, df_raw, method = "cosine",
#'                       mode = "number", n_sources = 10,
#'                       centroid_genes = 5, per = 1000)
#'
#' # --- Select 30 largest groups by gene count ---
#' sel <- select_sources(groups, df_raw, method = "size_only", n_sources = 30)
#' }
#'
#' @details
#' **Important note for proportion estimation:**  
#' The downstream function \code{proportion_rescaling()} requires each selected
#' source to contain exactly one representative marker gene.  
#' Therefore, when csCAM is used for proportion estimation, \strong{per must be set to 1}
#' so that each source contributes only its top-ranked gene.  
#' Using \code{per > 1} will produce groups containing multiple genes per source,
#' which makes proportion estimation invalid.
#' 
#' @export
select_sources <- function(groups, data,
                           method = c("cosine", "spearman", "size_only"),
                           mode = c("distance", "number"),
                           cutoff = 0.2, n_sources = 5,
                           per = 1, centroid_genes = NULL,
                           first = c("maxdist", "maxgenes")) {
  method <- match.arg(method)
  mode   <- match.arg(mode)
  first  <- match.arg(first)
  
  # Special case: method = "size_only" → return sources with most genes
  if (method == "size_only") {
    gene_counts <- sapply(groups, length)
    top_sources <- names(sort(gene_counts, decreasing = TRUE))[1:min(n_sources, length(groups))]
    result <- groups[top_sources]
    
    if (!is.null(per)) {
      result <- lapply(result, function(g) head(g, per))
    }
    return(result)
  }
  
  # 1. Compute centroids
  centroids <- t(sapply(groups, function(g) {
    if (!is.null(centroid_genes)) {
      g <- head(g, centroid_genes)
    }
    colMeans(data[g, , drop = FALSE])
  }))
  
  # 2. Compute distance matrix
  sim_mat <- switch(method,
                    cosine = {
                      num <- tcrossprod(centroids)
                      denom <- sqrt(rowSums(centroids^2))
                      denom_mat <- outer(denom, denom)
                      num / denom_mat
                    },
                    spearman = cor(t(centroids), method = "spearman", use = "pairwise.complete.obs")
  )
  dist_mat <- 1 - sim_mat
  
  # 3. Select sources
  if (mode == "distance") {
    selected <- character(0)
    for (src in rownames(centroids)) {
      if (length(selected) == 0) {
        selected <- c(selected, src)
      } else {
        dmin <- min(dist_mat[src, selected])
        if (dmin >= cutoff) {
          selected <- c(selected, src)
        }
      }
    }
  } else if (mode == "number") {
    remaining <- rownames(centroids)
    if (first == "maxgenes") {
      selected <- names(which.max(sapply(groups, length)))
    } else {
      selected <- remaining[which.max(rowSums(dist_mat))]  # default: most isolated
    }
    
    while (length(selected) < n_sources && length(setdiff(remaining, selected)) > 0) {
      candidates <- setdiff(remaining, selected)
      scores <- sapply(candidates, function(c) {
        min(dist_mat[c, selected])
      })
      selected <- c(selected, names(which.max(scores)))
    }
  }
  
  # 4. Subset groups
  result <- groups[selected]
  
  # 5. Limit per-source genes if requested
  if (!is.null(per)) {
    result <- lapply(result, function(g) head(g, per))
  }
  
  return(result)
}

#### Rank sources by maximal distinctness ####
#'
#' @title Rank sources by maximal distinctness
#'
#' @description
#' Iteratively ranks gene groups (candidate cell types) by maximizing
#' pairwise distance between their centroids. The algorithm selects the
#' first source as the most isolated group (or the largest one, depending
#' on \code{first}), and then successively adds groups that are maximally
#' distinct from all previously selected ones.  
#' 
#' This function can be viewed as an ordered version of
#' \code{select_sources()}, returning not only the final subset but also
#' the selection order, distance-based scores, and gene counts per group.
#' 
#' The optional parameter \code{max_n} allows early stopping after a fixed
#' number of cell types have been ranked, avoiding unnecessary computation
#' on large sets of groups.
#'
#' @param groups A list of gene groups, typically the output of
#'   \code{group_genes_by_pi()}, where each element is a character vector
#'   of gene names.
#' @param data A numeric expression matrix with samples in columns and
#'   genes in rows. Used to compute group centroids.
#' @param method Character, one of \code{"cosine"} or \code{"spearman"}.
#'   Determines how pairwise similarities between group centroids are
#'   calculated. Default = \code{"cosine"}.
#' @param first Character, one of \code{"maxdist"} or \code{"maxgenes"}.
#'   When \code{"maxdist"}, the first selected source is the most isolated
#'   centroid (with largest sum of pairwise distances); when
#'   \code{"maxgenes"}, it is the group with the largest number of genes.
#' @param centroid_genes Optional integer. Number of top genes to use per
#'   group when computing centroids (default = use all genes).
#' @param max_n Optional integer. Maximum number of cell types to rank.
#'   When provided, the algorithm stops once this number of sources has
#'   been selected.
#'
#' @return
#' A \code{data.frame} with four columns:
#' \itemize{
#'   \item \code{group} – Name of the gene group.
#'   \item \code{rank} – Selection order (1 = most representative).
#'   \item \code{score} – Minimum distance to previously selected groups
#'         at the time of selection (larger values indicate stronger
#'         distinctness).
#'   \item \code{n_genes} – Number of genes contained in the group
#'         (for reference only, not used in scoring).
#' }
#' Groups are ordered by \code{rank} in the output.
#'
#' @examples
#' \dontrun{
#' # Rank groups by cosine distance, keeping only the top 6 cell types
#' ranks <- rank_sources(groups, data, method = "cosine", max_n = 6)
#'
#' # View the ranked order and scores
#' head(ranks)
#' }
#'
#' @export
rank_sources <- function(groups, data,
                         method = c("cosine", "spearman"),
                         first = c("maxdist", "maxgenes"),
                         centroid_genes = NULL,
                         max_n = NULL) {
  method <- match.arg(method)
  first  <- match.arg(first)
  
  # ----- Step 1: Compute centroids for all groups -----
  centroids <- t(sapply(groups, function(g) {
    if (!is.null(centroid_genes)) {
      g <- head(g, centroid_genes)
    }
    colMeans(data[g, , drop = FALSE])
  }))
  
  # Gene counts per group
  gene_counts <- sapply(groups, length)
  
  # ----- Step 2: Compute pairwise similarity and distance -----
  sim_mat <- switch(method,
                    cosine = {
                      num <- tcrossprod(centroids)
                      denom <- sqrt(rowSums(centroids^2))
                      denom_mat <- outer(denom, denom)
                      num / denom_mat
                    },
                    spearman = cor(t(centroids), method = "spearman", use = "pairwise.complete.obs")
  )
  dist_mat <- 1 - sim_mat
  
  # ----- Step 3: Select the first group -----
  remaining <- rownames(centroids)
  if (first == "maxgenes") {
    selected <- names(which.max(gene_counts))
  } else {
    selected <- remaining[which.max(rowSums(dist_mat))]
  }
  
  # Initialize output table
  ranks <- data.frame(
    group = selected,
    rank = 1,
    score = Inf,
    n_genes = gene_counts[selected],
    stringsAsFactors = FALSE
  )
  
  # ----- Step 4: Iteratively select the most distinct groups -----
  while (length(selected) < length(remaining)) {
    if (!is.null(max_n) && length(selected) >= max_n) break
    candidates <- setdiff(remaining, selected)
    scores <- sapply(candidates, function(c) {
      min(dist_mat[c, selected])
    })
    best <- names(which.max(scores))
    selected <- c(selected, best)
    ranks <- rbind(ranks,
                   data.frame(group = best,
                              rank = length(selected),
                              score = max(scores),
                              n_genes = gene_counts[best],
                              stringsAsFactors = FALSE))
  }
  
  # ----- Step 5: Format and return -----
  ranks <- ranks[order(ranks$rank), ]
  rownames(ranks) <- NULL
  return(ranks)
}



#### K estimation ####
#'
#' @title Detect the primary elbow in a score–gene trade-off curve
#'
#' @description
#' Identify the strongest fold-change drop in model score jointly with a
#' reduction in the number of genes. Uses multiplicative (ratio-based)
#' changes rather than absolute differences.
#'
#' @param df A data frame with fixed columns:
#'   \itemize{
#'     \item \code{score} – numeric model score (to be minimized)
#'     \item \code{n_genes} – numeric gene count (complexity measure)
#'     \item \code{rank} – rank or order index
#'   }
#' @param gene_weight Numeric scalar controlling the influence of gene-count
#'   reduction in the composite metric. Default is 0.5.
#' @param min_rank Integer. The first rank to consider (skip early unstable
#'   values such as Inf baseline). Default is 2.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{primary_rank} – rank value of the largest composite drop
#'     \item \code{primary_index} – index of that rank in the data frame
#'     \item \code{composite} – full composite fold-change vector
#'     \item \code{fold_score} – raw fold change in score
#'     \item \code{fold_genes} – raw fold change in gene count
#'   }
#'
#' @examples
#' \dontrun{
#' result <- detect_primary_elbow(ranks, gene_weight = 0.5)
#' result$primary_rank}
#'
#' @export
detect_primary_elbow <- function(df,
                                 gene_weight = 0.5,
                                 min_rank = 2) {
  
  # === Extract fixed columns ===
  s <- df$score
  g <- df$n_genes
  r <- df$rank
  n <- length(s)
  
  # === Replace NA with finite values, but keep Inf for skipping ===
  s[is.na(s)] <- max(s[is.finite(s)], na.rm = TRUE)
  g[is.na(g)] <- max(g[is.finite(g)], na.rm = TRUE)
  
  # === Compute fold changes ===
  fold_score <- s[-n] / s[-1]       # ratio >1 means score decreased
  fold_gene  <- g[-n] / g[-1]       # ratio >1 means gene count decreased
  
  # === Replace invalid or infinite fold changes with NA ===
  invalid_idx <- !is.finite(fold_score) | !is.finite(fold_gene) |
    fold_score <= 0 | fold_gene <= 0
  fold_score[invalid_idx] <- NA
  fold_gene[invalid_idx]  <- NA
  
  # === Compute composite metric (log ratio combination) ===
  composite <- log(fold_score) + gene_weight * log(fold_gene)
  
  # === Pad NA for alignment (since fold change has n-1 elements) ===
  composite  <- c(composite, NA)
  fold_score <- c(fold_score, NA)
  fold_gene  <- c(fold_gene, NA)
  
  # === Ignore early unstable ranks (e.g., baseline Inf) ===
  if (min_rank > 1) {
    composite[1:(min_rank - 1)] <- NA
  }
  
  # === Find index of largest composite (strongest multiplicative drop) ===
  primary_idx <- which.max(composite)
  primary_rank <- r[primary_idx]
  
  # === Return ===
  list(
    primary_rank = primary_rank,
    primary_index = primary_idx,
    composite = composite,
    fold_score = fold_score,
    fold_genes = fold_gene
  )
}

#### Rescale A into proportions ####
#'
#' @title Rescale csCAM marker-gene expression into cell-type proportions
#'
#' @description
#' Construct an estimated proportion matrix from csCAM-selected marker genes.
#' For each source, one representative marker gene is used to form a sample-by-source
#' expression matrix, which is then column-rescaled using non-negative least squares
#' to best approximate row sums of 1 (valid proportion vectors).
#'
#' This produces:
#' \itemize{
#'   \item a full estimated proportion matrix (\code{A_csCAM})
#'   \item source-specific scaling coefficients (\code{source_scaler})
#'   \item marker gene identifiers actually used (\code{markers})
#' }
#'
#' @param X A numeric matrix of gene expression with genes in rows and samples in columns.
#' @param sel A list or character vector representing the csCAM-selected marker genes.
#'   If \code{sel} is a list (e.g., output of \code{select_sources}), the first gene in
#'   each element is used as the representative marker for that source.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{A_csCAM} – estimated sample-by-source proportion matrix
#'   \item \code{source_scaler} – numeric vector of scaling coefficients, one per source
#'   \item \code{markers} – character vector of marker genes actually used
#' }
#'
#' @examples
#' # Suppose sel is output from select_sources()
#' # cs_out <- proportion_rescaling(X, sel)
#' # A_hat <- cs_out$A_csCAM
#' # scalers <- cs_out$source_scaler
#'
#' @details
#' **Important note for proportion estimation:**  
#' The downstream function \code{proportion_rescaling()} requires each selected
#' source to contain exactly one representative marker gene.  
#' Therefore, when csCAM is used for proportion estimation, \strong{per must be set to 1}
#' so that each source contributes only its top-ranked gene.  
#' Using \code{per > 1} will produce groups containing multiple genes per source,
#' which makes proportion estimation invalid.
#'
#' @export
proportion_rescaling <- function(X, sel) {
  # 1. parse marker genes
  if (is.list(sel)) {
    marker_genes <- vapply(sel, function(g) g[1], character(1))
  } else {
    marker_genes <- as.character(sel)
  }
  
  # ensure markers exist in X
  marker_genes <- intersect(marker_genes, rownames(X))
  if (length(marker_genes) == 0L) {
    stop("No marker genes from 'sel' found in rownames(X).")
  }
  
  K <- length(marker_genes)
  
  # 2. expression of marker genes: sample × K
  V <- t(X[marker_genes, , drop = FALSE])
  V[V < 0] <- 0
  
  # 3. non-negative least squares to solve V %*% s ≈ 1
  b <- rep(1, nrow(V))
  nnls_fit <- nnls::nnls(as.matrix(V), b)
  s_hat <- nnls_fit$x
  
  # 4. construct A_csCAM = V * diag(s_hat)
  A_csCAM <- sweep(V, 2, s_hat, `*`)
  
  # (optional) ensure valid row sums
  # A_csCAM <- A_csCAM / pmax(rowSums(A_csCAM), 1e-12)
  
  # add dimnames
  rownames(A_csCAM) <- colnames(X)
  colnames(A_csCAM) <- paste0("csCAM_", seq_len(K))
  
  list(
    A_csCAM = A_csCAM,
    source_scaler = s_hat,
    markers = marker_genes
  )
}
