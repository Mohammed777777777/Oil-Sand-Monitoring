#PAPER 2 — TABLE T3 — DETECTION GRADIENT

hfi_year <- rasterize(
    g,
    fb,
    field = "hfi_yr",
    fun = "min",
    background = NA
  )

#Then it calculates:
  
  year_difference <- (
    fb -
      hfi_year
  )

temporal_match <- (
  reliable_loss &
    !is.na(hfi_year) &
    abs(year_difference) <= 2
)


library(terra)

# ============================================================
# PAPER 2 — TABLE T3
# DETECTION GRADIENT
#
# Final locked definitions:
#
# Any break:
#   first_break_year > 0
#
# Loss-type break:
#   any break AND change_type_loss > 0
#
# Reliable loss:
#   loss-type break AND last_segment_reliable > 0
#
# Temporal consistency:
#   |CCDC first break year - HFI year| <= 2 years
#
# IMPORTANT:
# Severity >= 0.249 is NOT used here.
# That threshold belongs to the independent
# high-confidence landscape characterization.
# ============================================================


# ============================================================
# 1. OUTPUT FOLDERS
# ============================================================

base <- "E:/NEW_CCDC_Augs/paper2_analysis"

dir.create(
  file.path(base, "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(base, "data"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. LOAD FINAL 39-BAND CCDC MOSAIC
# ============================================================

m <- rast(
  "E:/NEW_CCDC_Augs/ccdc_dense_summary_mosaic_04Aug26_MERGED.tif"
)

if (nlyr(m) != 39) {
  stop("ERROR: Expected the verified 39-band CCDC mosaic.")
}

cat("============================================\n")
cat(" PAPER 2 — DETECTION GRADIENT\n")
cat("============================================\n")

cat("Bands:", nlyr(m), "\n")
cat("Resolution:", res(m)[1], "x", res(m)[2], "m\n\n")


# ============================================================
# 3. LOCKED CCDC BAND INDICES
# ============================================================

# Band 1  = first_break_year
# Band 3  = change_count
# Band 31 = change_type_loss
# Band 38 = last_segment_reliable

fb   <- m[[1]]
cc   <- m[[3]]
loss <- m[[31]]
rel  <- m[[38]]

NAflag(fb)   <- -9999
NAflag(cc)   <- -9999
NAflag(loss) <- -9999
NAflag(rel)  <- -9999


# ============================================================
# 4. HFI GROUPS
# ============================================================

group_paths <- list(
  
  clearing =
    "E:/hfi/analysis_groups/clearing.shp",
  
  oilsands =
    "E:/hfi/analysis_groups/oilsands.shp",
  
  wellpad =
    "E:/hfi/analysis_groups/wellpads.shp"
)


# Known results from the validated earlier analysis
known_results <- data.frame(
  
  group = c(
    "clearing",
    "oilsands",
    "wellpad"
  ),
  
  covered_known = c(
    9621048,
    719513,
    535817
  ),
  
  any_known = c(
    92.5,
    98.0,
    80.8
  ),
  
  loss_known = c(
    89.1,
    94.7,
    70.3
  ),
  
  reliable_known = c(
    78.0,
    60.8,
    62.4
  ),
  
  temporal_all_known = c(
    67.1,
    41.4,
    41.8
  ),
  
  temporal_reliable_known = c(
    86.0,
    68.2,
    66.9
  )
)


# ============================================================
# 5. DETECTION FUNCTION
# ============================================================

analyze_detection <- function(
    shp,
    label
) {
  
  cat("\n--------------------------------------------\n")
  cat("Processing:", label, "\n")
  cat("--------------------------------------------\n")
  
  g <- vect(shp)
  
  # ----------------------------------------------------------
  # CRS alignment
  # ----------------------------------------------------------
  
  if (!same.crs(g, m)) {
    
    g <- project(
      g,
      crs(m)
    )
  }
  
  
  # ----------------------------------------------------------
  # Check HFI YEAR field
  # ----------------------------------------------------------
  
  if (!"YEAR" %in% names(g)) {
    
    stop(
      paste(
        "ERROR: YEAR field not found in",
        label
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Convert HFI year to numeric
  # ----------------------------------------------------------
  
  g$hfi_yr <- suppressWarnings(
    as.integer(
      as.character(
        g$YEAR
      )
    )
  )
  
  
  # ----------------------------------------------------------
  # Enforce HFI comparison period
  # 1985–2022
  # ----------------------------------------------------------
  
  keep <- (
    !is.na(g$hfi_yr) &
      g$hfi_yr >= 1985 &
      g$hfi_yr <= 2022
  )
  
  g <- g[keep, ]
  
  cat(
    "HFI features retained (1985-2022):",
    nrow(g),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Rasterize HFI footprint
  # ----------------------------------------------------------
  
  footprint <- rasterize(
    g,
    fb,
    field = 1,
    background = NA
  )
  
  inside <- !is.na(
    footprint
  )
  
  
  # ----------------------------------------------------------
  # Rasterize HFI year
  #
  # If polygons overlap, use earliest mapped HFI year.
  # This matches the previous analysis.
  # ----------------------------------------------------------
  
  hfi_year <- rasterize(
    g,
    fb,
    field = "hfi_yr",
    fun = "min",
    background = NA
  )
  
  
  # ==========================================================
  # DETECTION DEFINITIONS
  # ==========================================================
  
  # Any CCDC break
  any_break <- (
    inside &
      !is.na(fb) &
      fb > 0
  )
  
  
  # Loss-type break
  loss_break <- (
    any_break &
      !is.na(loss) &
      loss > 0
  )
  
  
  # Reliable loss
  reliable_loss <- (
    loss_break &
      !is.na(rel) &
      rel > 0
  )
  
  
  # ----------------------------------------------------------
  # ±2-year timing consistency
  # ----------------------------------------------------------
  
  year_difference <- (
    fb -
      hfi_year
  )
  
  temporal_match <- (
    reliable_loss &
      !is.na(hfi_year) &
      abs(year_difference) <= 2
  )
  
  
  # ==========================================================
  # PIXEL COUNTS
  # ==========================================================
  
  n_covered <- global(
    ifel(inside, 1, 0),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  
  n_any <- global(
    ifel(any_break, 1, 0),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  
  n_loss <- global(
    ifel(loss_break, 1, 0),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  
  n_reliable <- global(
    ifel(reliable_loss, 1, 0),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  
  n_temporal <- global(
    ifel(temporal_match, 1, 0),
    "sum",
    na.rm = TRUE
  )[1,1]
  
  
  # ==========================================================
  # PERCENTAGES
  # ==========================================================
  
  pct_any <- (
    100 *
      n_any /
      n_covered
  )
  
  
  pct_loss <- (
    100 *
      n_loss /
      n_covered
  )
  
  
  pct_reliable <- (
    100 *
      n_reliable /
      n_covered
  )
  
  
  # Temporal match relative to ALL HFI pixels
  pct_temporal_all <- (
    100 *
      n_temporal /
      n_covered
  )
  
  
  # Temporal match relative to RELIABLE-LOSS pixels
  pct_temporal_reliable <- ifelse(
    n_reliable > 0,
    100 *
      n_temporal /
      n_reliable,
    NA_real_
  )
  
  
  # ==========================================================
  # PRINT
  # ==========================================================
  
  cat(
    "Covered pixels                    :",
    format(
      n_covered,
      big.mark = ","
    ),
    "\n"
  )
  
  cat(
    sprintf(
      "Any CCDC break                   : %.1f%%\n",
      pct_any
    )
  )
  
  cat(
    sprintf(
      "Loss-type break                  : %.1f%%\n",
      pct_loss
    )
  )
  
  cat(
    sprintf(
      "Reliable loss                    : %.1f%%\n",
      pct_reliable
    )
  )
  
  cat(
    sprintf(
      "Temporal match ±2 yr / all       : %.1f%%\n",
      pct_temporal_all
    )
  )
  
  cat(
    sprintf(
      "Temporal match ±2 yr / reliable  : %.1f%%\n",
      pct_temporal_reliable
    )
  )
  
  
  # ==========================================================
  # RETURN TABLE ROW
  # ==========================================================
  
  data.frame(
    
    group = label,
    
    hfi_features =
      nrow(g),
    
    covered_pixels =
      n_covered,
    
    any_break_pixels =
      n_any,
    
    any_break_percent =
      pct_any,
    
    loss_break_pixels =
      n_loss,
    
    loss_break_percent =
      pct_loss,
    
    reliable_loss_pixels =
      n_reliable,
    
    reliable_loss_percent =
      pct_reliable,
    
    temporal_match_pixels =
      n_temporal,
    
    temporal_match_all_percent =
      pct_temporal_all,
    
    temporal_match_reliable_percent =
      pct_temporal_reliable
  )
}


# ============================================================
# 6. RUN THE THREE DISTURBANCE GROUPS
# ============================================================

r1 <- analyze_detection(
  group_paths$clearing,
  "clearing"
)

r2 <- analyze_detection(
  group_paths$oilsands,
  "oilsands"
)

r3 <- analyze_detection(
  group_paths$wellpad,
  "wellpad"
)

T3 <- rbind(
  r1,
  r2,
  r3
)

rownames(T3) <- NULL


# ============================================================
# 7. COMPARE AGAINST LOCKED RESULTS
# ============================================================

check <- merge(
  T3,
  known_results,
  by = "group"
)

check$delta_any <- (
  check$any_break_percent -
    check$any_known
)

check$delta_loss <- (
  check$loss_break_percent -
    check$loss_known
)

check$delta_reliable <- (
  check$reliable_loss_percent -
    check$reliable_known
)

check$delta_temporal_all <- (
  check$temporal_match_all_percent -
    check$temporal_all_known
)

check$delta_temporal_reliable <- (
  check$temporal_match_reliable_percent -
    check$temporal_reliable_known
)


# ============================================================
# 8. PRINT FINAL TABLE
# ============================================================

cat("\n\n============================================\n")
cat(" TABLE T3 — DETECTION GRADIENT\n")
cat("============================================\n")

T3_print <- T3[
  ,
  c(
    "group",
    "hfi_features",
    "covered_pixels",
    "any_break_percent",
    "loss_break_percent",
    "reliable_loss_percent",
    "temporal_match_all_percent",
    "temporal_match_reliable_percent"
  )
]

print(
  T3_print,
  digits = 4
)


# ============================================================
# 9. SELF-CHECK
# ============================================================

cat("\n============================================\n")
cat(" SELF-CHECK AGAINST EARLIER LOCKED RESULTS\n")
cat("============================================\n")

print(
  check[
    ,
    c(
      "group",
      "delta_any",
      "delta_loss",
      "delta_reliable",
      "delta_temporal_all",
      "delta_temporal_reliable"
    )
  ],
  digits = 4
)


max_difference <- max(
  abs(
    c(
      check$delta_any,
      check$delta_loss,
      check$delta_reliable,
      check$delta_temporal_all,
      check$delta_temporal_reliable
    )
  ),
  na.rm = TRUE
)

if (max_difference <= 0.5) {
  
  cat(
    "\nSELF-CHECK PASSED.\n"
  )
  
  cat(
    "Results reproduce the earlier detection-gradient analysis within 0.5 percentage points.\n"
  )
  
} else {
  
  cat(
    "\nWARNING: SELF-CHECK DID NOT FULLY MATCH.\n"
  )
  
  cat(
    "Do NOT change the definitions automatically.\n"
  )
  
  cat(
    "Send me the printed results so we can identify exactly which quantity differs.\n"
  )
}


# ============================================================
# 10. SAVE FULL TABLE
# ============================================================

full_file <- file.path(
  base,
  "tables",
  "T3_detection_gradient_FULL.csv"
)

write.csv(
  T3,
  full_file,
  row.names = FALSE
)


# ============================================================
# 11. SAVE PAPER-READY TABLE
# ============================================================

T3_paper <- data.frame(
  
  Disturbance_group = c(
    "Vegetation clearing",
    "Oil-sands industrial",
    "Well pads"
  ),
  
  HFI_features =
    T3$hfi_features,
  
  Covered_pixels =
    T3$covered_pixels,
  
  Any_break_percent =
    round(
      T3$any_break_percent,
      1
    ),
  
  Loss_type_break_percent =
    round(
      T3$loss_break_percent,
      1
    ),
  
  Reliable_loss_percent =
    round(
      T3$reliable_loss_percent,
      1
    ),
  
  Timing_consistency_2yr_percent =
    round(
      T3$temporal_match_reliable_percent,
      1
    )
)

paper_file <- file.path(
  base,
  "tables",
  "T3_detection_gradient_PAPER.csv"
)

write.csv(
  T3_paper,
  paper_file,
  row.names = FALSE
)


# ============================================================
# 12. SAVE VALIDATION CHECK
# ============================================================

check_file <- file.path(
  base,
  "tables",
  "T3_detection_gradient_selfcheck.csv"
)

write.csv(
  check,
  check_file,
  row.names = FALSE
)


# ============================================================
# 13. FINAL MESSAGE
# ============================================================

cat("\n============================================\n")
cat(" DETECTION GRADIENT COMPLETE\n")
cat("============================================\n")

cat("\nSaved:\n")

cat(
  full_file,
  "\n"
)

cat(
  paper_file,
  "\n"
)

cat(
  check_file,
  "\n"
)

cat("\nDONE.\n")