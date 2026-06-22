#' Generate Catch Advice for Management Strategy Evaluation
#'
#' Generate catch advice based on harvest control rules (HCRs) using a fitted WHAM
#' estimation model. This function extends \code{\link{project_wham}} by embedding
#' HCR-specific logic (F\%SPR, constant catch, hockey-stick), while also allowing
#' users to supply custom projection options via \code{proj.opts} and optional
#' environmental covariate projection controls via \code{ecov_em_opts}.
#'
#' @param em A fitted WHAM estimation model object.
#' @param pro.yr Integer. Number of years for projection. Default is taken from
#'   \code{assess.interval}.
#' @param hcr A list specifying the harvest control rule:
#'   \describe{
#'     \item{\code{hcr.type}}{Integer specifying the HCR type:
#'       \enumerate{
#'         \item F\_XSPR (or F\_MSY if \code{use_FMSY = TRUE}).
#'         \item Constant catch.
#'         \item Hockey-stick scaling based on biomass thresholds.
#'         \item MAFMC P* rule
#'       }}
#'     \item{\code{hcr.opts}}{List of options controlling HCR behavior:
#'       \itemize{
#'         \item \code{use_FXSPR}, \code{percentFXSPR}
#'         \item \code{use_FMSY}, \code{percentFMSY}
#'         \item \code{avg_yrs}, \code{cont.M.re}, \code{cont.move.re}
#'         \item \code{max_percent}, \code{min_percent},
#'           \code{BThresh_up}, \code{BThresh_low} (for HCR 3)
#'         \item \code{max_pstar}, \code{mid_pstar}, \code{min_pstar},
#'           \code{BThresh_high}, \code{OFL_CV} (for HCR 4)
#'       }}
#'   }
#' @param proj.opts A named list of projection options passed directly to
#'   \code{\link{project_wham}}. Any field here overrides values set internally
#'   from defaults, HCR settings, and \code{ecov_em_opts}.
#' @param ecov_em_opts List (optional). Options for projecting environmental
#'   covariates in the estimation model during catch advice. Expected components:
#'   \itemize{
#'     \item \code{use_ecov_em}: Logical. If \code{TRUE}, use EM-based/projected Ecov values.
#'     \item \code{lag}: Integer lag used to align Ecov values.
#'     \item \code{period}: Integer vector (optional). Projection indices to override.
#'   }
#'   If \code{NULL}, environmental covariates are handled through the default
#'   projection settings or any user-supplied \code{proj.opts}.
#'
#' @details
#' Projection options are constructed in four steps:
#' \enumerate{
#'   \item Defaults are set internally.
#'   \item HCR logic modifies relevant fields.
#'   \item Optional \code{ecov_em_opts} modifies \code{proj.ecov}.
#'   \item User-supplied \code{proj.opts} overrides all previous settings.
#' }
#'
#' \strong{HCR 1}: Uses F\_XSPR (or F\_MSY if requested) and returns projected catch.
#'
#' \strong{HCR 2}: Projects first, then repeats the mean projected catch across
#' all projection years.
#'
#' \strong{HCR 3}: Hockey-stick rule. The terminal biomass ratio
#' \eqn{SSB_t / SSB_x} is mapped linearly to \code{percentFXSPR} or
#' \code{percentFMSY} between \code{BThresh_low} and \code{BThresh_up}, bounded by
#' \code{min_percent} and \code{max_percent}, and catch advice is then projected
#' under that scaled reference point.
#' 
#' \strong{HCR 4}: MAFMC rule. The terminal biomass ratio
#' \eqn{SSB_t / SSB_x} is mapped to the risk policy specified probability of 
#' overfishing (pstar) based on three B thresholds and three pstar levels:
#' \code{min_pstar} is specified below \code{BThresh_low}, a linear ramp between
#' \code{min_pstar} and \code{mid_pstar} between \code{BThresh_low} and 
#' \code{BThresh_up}, a linear ramp between \code{mid_pstar} and \code{max_pstar}
#' between \code{BThresh_up} and \code{BThresh_high}, and \code{max_pstar} above
#' \code{BThresh_high}. The OFL is calculated as in HCR 1, then the ABC (catch
#' advice) is calculated using the resulting pstar as the percentile of OFL 
#' assuming that the OFL is lognormally distributed with CV \code{OFL_CV}. Default
#' values are 0.0, 0.45, and 0.49 for min, mid, and max pstar, 0.1, 1.0, and 1.5
#' time B_msy for min, up, and high B thresholds, and 1.0 for OFL CV.
#' 
#'
#' @return A matrix of projected catch advice for \code{pro.yr} years.
#'
#' @examples
#' # Example 1: F40%SPR advice
#' hcr <- list(
#'   hcr.type = 1,
#'   hcr.opts = list(
#'     use_FXSPR = TRUE,
#'     percentFXSPR = 100,
#'     avg_yrs = 5,
#'     cont.M.re = FALSE,
#'     cont.move.re = FALSE
#'   )
#' )
#' # advice <- advice_fn(em, pro.yr = 3, hcr = hcr)
#'
#' # Example 2: user override with proj.opts
#' # advice <- advice_fn(
#' #   em, pro.yr = 3, hcr = hcr,
#' #   proj.opts = list(proj_F_opt = c(4, 4, 4), proj_Fcatch = c(0.1, 0.2, 0.3))
#' # )
#'
#' # Example 3: EM Ecov projection override
#' # advice <- advice_fn(
#' #   em, pro.yr = 5, hcr = hcr,
#' #   ecov_em_opts = list(use_ecov_em = TRUE, lag = 1)
#' # )
#'
#' @seealso \code{\link{project_wham}}, \code{\link{fit_wham}}
#' @export
advice_fn <- function(em,
                      pro.yr = 1,
                      hcr = NULL,
                      proj.opts = list(),
                      ecov_em_opts = NULL) {
  
  ## ------------------------------
  ## internal helpers
  ## ------------------------------
  get_last_proj_catch <- function(em_proj, pro.yr) {
    pc <- em_proj$rep$pred_catch
    
    if (is.null(dim(pc))) {
      pc <- matrix(pc, ncol = 1)
    } else {
      pc <- as.matrix(pc)
    }
    
    nr <- nrow(pc)
    if (pro.yr > nr) {
      stop("pro.yr exceeds the number of available rows in em_proj$rep$pred_catch.",
           call. = FALSE)
    }
    
    out <- pc[(nr - pro.yr + 1):nr, , drop = FALSE]
    return(out)
  }
  
  validate_single_fspec <- function(proj_opts, hcr.type) {
    flags <- c(
      use.last.F = isTRUE(proj_opts$use.last.F),
      use.avg.F  = isTRUE(proj_opts$use.avg.F),
      use.FXSPR  = isTRUE(proj_opts$use.FXSPR),
      use.FMSY   = isTRUE(proj_opts$use.FMSY),
      proj.F     = !is.null(proj_opts$proj.F),
      proj.catch = !is.null(proj_opts$proj.catch)
    )
    
    if (sum(flags) != 1) {
      stop(
        sprintf(
          "Exactly one F/catch specification must be active for HCR type %s. Active settings: %s",
          hcr.type,
          paste(names(flags)[flags], collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  
  print_proj_opts <- function(proj_opts, title = "Final Projection Options:") {
    cat(title, "\n", sep = "")
    for (nm in names(proj_opts)) {
      if (!is.null(proj_opts[[nm]])) {
        cat(sprintf(" %s: %s\n", nm, toString(proj_opts[[nm]])))
      }
    }
  }
  
  ## ------------------------------
  ## Step 1. normalize inputs
  ## ------------------------------
  if (is.null(hcr)) {
    hcr <- list(hcr.type = 1, hcr.opts = list())
  }
  
  hcr.type <- if (is.null(hcr$hcr.type)) 1 else hcr$hcr.type
  hcr.opts <- if (is.null(hcr$hcr.opts)) list() else hcr$hcr.opts
  
  if (!hcr.type %in% c(1, 2, 3, 4)) {
    stop("hcr$hcr.type must be one of 1, 2, or 3.", call. = FALSE)
  }
  
  cat(paste0("\nHarvest Control Rule type ", hcr.type, "\n"))
  
  ## ------------------------------
  ## Step 2. HCR-level defaults
  ## ------------------------------
  use_FXSPR    <- if (is.null(hcr.opts$use_FXSPR)) TRUE else hcr.opts$use_FXSPR
  percentFXSPR <- if (is.null(hcr.opts$percentFXSPR)) 75 else hcr.opts$percentFXSPR
  use_FMSY     <- if (is.null(hcr.opts$use_FMSY)) FALSE else hcr.opts$use_FMSY
  percentFMSY  <- if (is.null(hcr.opts$percentFMSY)) 75 else hcr.opts$percentFMSY
  avg_yrs_n    <- if (is.null(hcr.opts$avg_yrs)) 5 else hcr.opts$avg_yrs
  
  if (length(avg_yrs_n) != 1 || !is.numeric(avg_yrs_n) || is.na(avg_yrs_n) || avg_yrs_n < 1) {
    stop("hcr.opts$avg_yrs must be a positive integer.", call. = FALSE)
  }
  avg_yrs_n <- min(as.integer(avg_yrs_n), length(em$years))
  
  cont.M.re <- if (is.null(hcr.opts$cont.M.re)) FALSE else hcr.opts$cont.M.re
  cont.move.re <- if (is.null(hcr.opts$cont.move.re) || em$input$data$n_regions == 1) {
    NULL
  } else {
    hcr.opts$cont.move.re
  }
  
  ## ------------------------------
  ## Step 3. default projection options
  ## ------------------------------
  defaults <- list(
    n.yrs          = pro.yr,
    use.last.F     = FALSE,
    use.avg.F      = FALSE,
    use.FXSPR      = use_FXSPR,
    percentFXSPR   = percentFXSPR,
    use.FMSY       = use_FMSY,
    percentFMSY    = percentFMSY,
    proj.F         = NULL,
    proj.catch     = NULL,
    avg.yrs        = tail(em$years, avg_yrs_n),
    cont.ecov      = TRUE,
    use.last.ecov  = FALSE,
    avg.ecov.yrs   = NULL,
    proj.ecov      = NULL,
    cont.M.re      = cont.M.re,
    cont.move.re   = cont.move.re,
    cont.L.re      = FALSE,
    avg.rec.yrs    = NULL,
    proj_F_opt     = NULL,
    proj_Fcatch    = NULL,
    proj_mature    = NULL,
    proj_waa       = NULL,
    proj_R_opt     = NULL,
    proj_NAA_opt   = NULL,
    proj_NAA_init  = NULL,
    proj_F_init    = NULL,
    avg.yrs.sel      = NULL,
    avg.yrs.waacatch = NULL,
    avg.yrs.waassb   = NULL,
    avg.yrs.mature   = NULL,
    avg.yrs.L        = NULL,
    avg.yrs.M        = NULL,
    avg.yrs.move     = NULL,
    avg.yrs.R        = NULL,
    avg.yrs.NAA      = NULL
  )
  
  ## ------------------------------
  ## Step 4. optional EM Ecov projection logic
  ## ------------------------------
  if (!is.null(ecov_em_opts) && isTRUE(ecov_em_opts$use_ecov_em)) {
    
    if (is.null(ecov_em_opts$lag)) {
      stop("Must specify ecov_em_opts$lag when ecov_em_opts$use_ecov_em = TRUE.",
           call. = FALSE)
    }
    
    if (is.null(em$rep$Ecov_x)) {
      stop("ecov_em_opts$use_ecov_em = TRUE but em$rep$Ecov_x is not available.",
           call. = FALSE)
    }
    
    lag_val <- as.integer(ecov_em_opts$lag)
    ecov_x <- em$rep$Ecov_x
    
    if (is.null(dim(ecov_x))) {
      ecov_x <- matrix(ecov_x, ncol = 1)
    } else {
      ecov_x <- as.matrix(ecov_x)
    }
    
    start_id <- length(em$input$years) - lag_val
    end_id   <- start_id + pro.yr - 1
    
    if (start_id < 1 || end_id > nrow(ecov_x)) {
      stop("Requested Ecov projection period is outside the available range of em$rep$Ecov_x.",
           call. = FALSE)
    }
    
    proj_ecov <- ecov_x[start_id:end_id, , drop = FALSE]
    
    if (!is.null(ecov_em_opts$period)) {
      period_id <- ecov_em_opts$period
      if (any(period_id < 1 | period_id > pro.yr)) {
        stop("ecov_em_opts$period must contain indices between 1 and pro.yr.",
             call. = FALSE)
      }
      tmp_ecov <- defaults$proj.ecov
      if (is.null(tmp_ecov)) {
        tmp_ecov <- matrix(NA_real_, nrow = pro.yr, ncol = ncol(proj_ecov))
      }
      tmp_ecov[period_id, ] <- proj_ecov[period_id, , drop = FALSE]
      defaults$proj.ecov <- tmp_ecov
    } else {
      defaults$proj.ecov <- proj_ecov
    }
  }
  
  ## ------------------------------
  ## Step 5. merge user overrides
  ## user proj.opts has highest priority
  ## ------------------------------
  proj_opts <- modifyList(defaults, proj.opts)
  
  if (!is.null(proj.opts$avg.yrs) &&
      is.numeric(proj.opts$avg.yrs) &&
      length(proj.opts$avg.yrs) == 1) {
    n_tail <- min(as.integer(proj.opts$avg.yrs), length(em$years))
    proj_opts$avg.yrs <- tail(em$years, n_tail)
  }
  
  proj_opts <- Filter(Negate(is.null), proj_opts)
  
  ## ------------------------------
  ## Step 6. HCR logic
  ## ------------------------------
  if (hcr.type %in% c(1, 2)) {
    
    validate_single_fspec(proj_opts, hcr.type = hcr.type)
    
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
    advice_mat <- get_last_proj_catch(em_proj, pro.yr = pro.yr)
    
    if (hcr.type == 1) {
      advice <- advice_mat
    }
    
    if (hcr.type == 2) {
      catch_mean <- if (nrow(advice_mat) == 1) {
        as.numeric(advice_mat[1, ])
      } else {
        colMeans(advice_mat)
      }
      
      advice <- matrix(
        rep(catch_mean, each = pro.yr),
        nrow = pro.yr,
        byrow = TRUE
      )
      colnames(advice) <- colnames(advice_mat)
      rownames(advice) <- NULL
    }
  }
  
  if (hcr.type == 3) {
    
    max_percent <- if (is.null(hcr.opts$max_percent)) 75 else hcr.opts$max_percent
    min_percent <- if (is.null(hcr.opts$min_percent)) 0.01 else hcr.opts$min_percent
    BThresh_up  <- if (is.null(hcr.opts$BThresh_up)) 0.5 else hcr.opts$BThresh_up
    BThresh_low <- if (is.null(hcr.opts$BThresh_low)) 0.1 else hcr.opts$BThresh_low
    
    if (BThresh_low >= BThresh_up) {
      stop("BThresh_low must be smaller than BThresh_up.", call. = FALSE)
    }
    
    if (isTRUE(use_FXSPR) && isTRUE(use_FMSY)) {
      stop("For HCR type 3, choose only one reference path: use_FXSPR = TRUE or use_FMSY = TRUE.",
           call. = FALSE)
    }
    
    if (!isTRUE(use_FXSPR) && !isTRUE(use_FMSY)) {
      stop("For HCR type 3, one of use_FXSPR or use_FMSY must be TRUE.",
           call. = FALSE)
    }
    
    if (isTRUE(use_FXSPR)) {
      
      if (is.null(em$rep$log_SSB_FXSPR)) {
        stop("HCR type 3 with use_FXSPR = TRUE requires em$rep$log_SSB_FXSPR.",
             call. = FALSE)
      }
      
      if (is.null(em$rep$SSB)) {
        stop("HCR type 3 requires em$rep$SSB.", call. = FALSE)
      }
      
      SSB_x <- exp(em$rep$log_SSB_FXSPR[nrow(em$rep$log_SSB_FXSPR),
                                        ncol(em$rep$log_SSB_FXSPR)])
      SSB_t <- sum(em$rep$SSB[nrow(em$rep$SSB), ])
      ratio <- SSB_t / SSB_x
      
      if (ratio >= BThresh_up) {
        percent <- max_percent
      } else if (ratio > BThresh_low) {
        slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
        percent <- slope * (ratio - BThresh_low) + min_percent
      } else {
        percent <- min_percent
      }
      
      cat(sprintf("SSB_t / SSB_XSPR = %.3f -> percentFXSPR = %.2f\n", ratio, percent))
      
      proj_opts$use.last.F   <- FALSE
      proj_opts$use.avg.F    <- FALSE
      proj_opts$use.FXSPR    <- TRUE
      proj_opts$use.FMSY     <- FALSE
      proj_opts$proj.F       <- NULL
      proj_opts$proj.catch   <- NULL
      proj_opts$percentFXSPR <- as.numeric(percent)
      proj_opts$percentFMSY  <- NULL
    }
    
    if (isTRUE(use_FMSY)) {
      
      if (is.null(em$rep$log_SSB_MSY)) {
        stop("HCR type 3 with use_FMSY = TRUE requires em$rep$log_SSB_MSY.",
             call. = FALSE)
      }
      
      if (is.null(em$rep$SSB)) {
        stop("HCR type 3 requires em$rep$SSB.", call. = FALSE)
      }
      
      SSB_x <- exp(em$rep$log_SSB_MSY[nrow(em$rep$log_SSB_MSY),
                                      ncol(em$rep$log_SSB_MSY)])
      SSB_t <- sum(em$rep$SSB[nrow(em$rep$SSB), ])
      ratio <- SSB_t / SSB_x
      
      if (ratio >= BThresh_up) {
        percent <- max_percent
      } else if (ratio > BThresh_low) {
        slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
        percent <- slope * (ratio - BThresh_low) + min_percent
      } else {
        percent <- min_percent
      }
      
      cat(sprintf("SSB_t / SSB_MSY = %.3f -> percentFMSY = %.2f\n", ratio, percent))
      
      proj_opts$use.last.F   <- FALSE
      proj_opts$use.avg.F    <- FALSE
      proj_opts$use.FXSPR    <- FALSE
      proj_opts$use.FMSY     <- TRUE
      proj_opts$proj.F       <- NULL
      proj_opts$proj.catch   <- NULL
      proj_opts$percentFMSY  <- as.numeric(percent)
      proj_opts$percentFXSPR <- NULL
    }
    
    validate_single_fspec(proj_opts, hcr.type = hcr.type)
    
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
    advice <- get_last_proj_catch(em_proj, pro.yr = pro.yr)
  }
  
  if (hcr.type == 4) {
    
    max_pstar <- if (is.null(hcr.opts$max_pstar)) 0.49 else hcr.opts$max_pstar
    mid_pstar <- if (is.null(hcr.opts$mid_pstar)) 0.45 else hcr.opts$mid_pstar
    min_pstar <- if (is.null(hcr.opts$min_pstar)) 0.0 else hcr.opts$min_pstar
    
    BThresh_high <- if (is.null(hcr.opts$BThresh_high)) 1.5 else hcr.opts$BThresh_high
    BThresh_up  <- if (is.null(hcr.opts$BThresh_up)) 1.0 else hcr.opts$BThresh_up
    BThresh_low <- if (is.null(hcr.opts$BThresh_low)) 0.1 else hcr.opts$BThresh_low
    
    if (BThresh_low >= BThresh_up) {
      stop("BThresh_low must be smaller than BThresh_up.", call. = FALSE)
    }

    if (BThresh_up >= BThresh_high) {
      stop("BThresh_up must be smaller than BThresh_high.", call. = FALSE)
    }

    
    ############### UNMODIFIED BELOW HERE #######################
    
    if (isTRUE(use_FXSPR) && isTRUE(use_FMSY)) {
      stop("For HCR type 3, choose only one reference path: use_FXSPR = TRUE or use_FMSY = TRUE.",
           call. = FALSE)
    }
    
    if (!isTRUE(use_FXSPR) && !isTRUE(use_FMSY)) {
      stop("For HCR type 3, one of use_FXSPR or use_FMSY must be TRUE.",
           call. = FALSE)
    }
    
    if (isTRUE(use_FXSPR)) {
      
      if (is.null(em$rep$log_SSB_FXSPR)) {
        stop("HCR type 3 with use_FXSPR = TRUE requires em$rep$log_SSB_FXSPR.",
             call. = FALSE)
      }
      
      if (is.null(em$rep$SSB)) {
        stop("HCR type 3 requires em$rep$SSB.", call. = FALSE)
      }
      
      SSB_x <- exp(em$rep$log_SSB_FXSPR[nrow(em$rep$log_SSB_FXSPR),
                                        ncol(em$rep$log_SSB_FXSPR)])
      SSB_t <- sum(em$rep$SSB[nrow(em$rep$SSB), ])
      ratio <- SSB_t / SSB_x
      
      if (ratio >= BThresh_up) {
        percent <- max_percent
      } else if (ratio > BThresh_low) {
        slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
        percent <- slope * (ratio - BThresh_low) + min_percent
      } else {
        percent <- min_percent
      }
      
      cat(sprintf("SSB_t / SSB_XSPR = %.3f -> percentFXSPR = %.2f\n", ratio, percent))
      
      proj_opts$use.last.F   <- FALSE
      proj_opts$use.avg.F    <- FALSE
      proj_opts$use.FXSPR    <- TRUE
      proj_opts$use.FMSY     <- FALSE
      proj_opts$proj.F       <- NULL
      proj_opts$proj.catch   <- NULL
      proj_opts$percentFXSPR <- as.numeric(percent)
      proj_opts$percentFMSY  <- NULL
    }
    
    if (isTRUE(use_FMSY)) {
      
      if (is.null(em$rep$log_SSB_MSY)) {
        stop("HCR type 3 with use_FMSY = TRUE requires em$rep$log_SSB_MSY.",
             call. = FALSE)
      }
      
      if (is.null(em$rep$SSB)) {
        stop("HCR type 3 requires em$rep$SSB.", call. = FALSE)
      }
      
      SSB_x <- exp(em$rep$log_SSB_MSY[nrow(em$rep$log_SSB_MSY),
                                      ncol(em$rep$log_SSB_MSY)])
      SSB_t <- sum(em$rep$SSB[nrow(em$rep$SSB), ])
      ratio <- SSB_t / SSB_x
      
      if (ratio >= BThresh_up) {
        percent <- max_percent
      } else if (ratio > BThresh_low) {
        slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
        percent <- slope * (ratio - BThresh_low) + min_percent
      } else {
        percent <- min_percent
      }
      
      cat(sprintf("SSB_t / SSB_MSY = %.3f -> percentFMSY = %.2f\n", ratio, percent))
      
      proj_opts$use.last.F   <- FALSE
      proj_opts$use.avg.F    <- FALSE
      proj_opts$use.FXSPR    <- FALSE
      proj_opts$use.FMSY     <- TRUE
      proj_opts$proj.F       <- NULL
      proj_opts$proj.catch   <- NULL
      proj_opts$percentFMSY  <- as.numeric(percent)
      proj_opts$percentFXSPR <- NULL
    }
    
    validate_single_fspec(proj_opts, hcr.type = hcr.type)
    
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
    advice <- get_last_proj_catch(em_proj, pro.yr = pro.yr)
    
    
  }
    
  
  ## ------------------------------
  ## Step 7. print and return
  ## ------------------------------
  print_proj_opts(proj_opts)
  
  return(advice)
}