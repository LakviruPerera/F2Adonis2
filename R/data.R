#' Polycystic Kidney Disease (PKD) Metadata
#'
#' A dataset containing group categories for the Polycystic Kidney Disease (PKD)
#' cohort, used to demonstrate multivariate testing under biological heterogeneity.
#'
#' @format A data frame with 95 rows (samples) and 1 column:
#' \describe{
#'   \item{Group}{A factor indicating disease classification status (Ctrl vs PKD).}
#' }
#' @source Derived from the urinary metabolomics cohort in Houske et al. (2023).
#' @keywords datasets
"pkd_metadata"

#' Polycystic Kidney Disease (PKD) Metabolites
#'
#' An auto-scaled, log-transformed intensity matrix of urinary profiles (metabolites) from
#' early-stage Polycystic Kidney Disease patients and healthy controls.
#'
#' @format A matrix with 95 rows (samples) and 1554 columns (metabolites).
#' @source Derived from the urinary metabolomics cohort in Houske et al. (2023).
#' @keywords datasets
"pkd_metabolites"

#' Norway Benthic Macrofauna Presence-Absence Data
#'
#' A presence-absence matrix for 809 marine macrofauna species across 101 sites,
#' utilized as a baseline biological validation dataset.
#'
#' @format A matrix with 101 rows (sites) and 809 columns (species).
#' @source \doi{10.1111/anzs.12176}
#' @keywords datasets
"norway_pa"

#' Norway Benthic Macrofauna Metadata
#'
#' Metadata containing geographical location assignments for the 101 sampled sites
#' in the Norway benthic community dataset.
#'
#' @format A data frame with 101 rows and 1 column:
#' \describe{
#'   \item{Area}{A grouping factor with 5 levels indicating the specific geographical region.}
#' }
#' @source \doi{10.1111/anzs.12176}
#' @keywords datasets
"norway_metadata"
