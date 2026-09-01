# Family line plots for ASVs with LOO mean importance > 0.001 in ALL 5 of
#   clr_0.1, clr_0.5, clr_1, rclr, rclr_optspace.
# One PDF per transform (clr_1, rclr, clr_mean4).
#
# Run: Rscript R/21_family_lineplots_all5_gt0.001.R

lines <- readLines("R/16_rclr_family_lineplots.R")
cut <- which(grepl("^# ---- CLR \\(pc=1\\)", lines))[1]
if (is.na(cut) || cut < 2L) {
  stop("Could not find CLR run block in R/16_rclr_family_lineplots.R")
}
eval(parse(text = lines[seq_len(cut - 1L)]), envir = environment())

TRANSFORM_ALL5 <- c("clr_0.1", "clr_0.5", "clr_1", "rclr", "rclr_optspace")
MIN_TRANSFORMS <- length(TRANSFORM_ALL5)

message(glue(
  "Subset: mean_loo_imp > {IMP_THRESH} in {MIN_TRANSFORMS}/{length(TRANSFORM_ALL5)} ",
  "(all 5 transforms); facets = Family"
))

res_16s_clr <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_ALL5,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_1",
  pdf_name = "clr1_family_16S_all5_gt0.001.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_clr <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_ALL5,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_1",
  pdf_name = "clr1_family_18S_all5_gt0.001.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

res_16s_rclr <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_ALL5,
  min_transforms = MIN_TRANSFORMS,
  method = "rclr",
  pdf_name = "rclr_family_16S_all5_gt0.001.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_rclr <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_ALL5,
  min_transforms = MIN_TRANSFORMS,
  method = "rclr",
  pdf_name = "rclr_family_18S_all5_gt0.001.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

res_16s_sum4 <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_ALL5,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_sum4",
  pdf_name = "clr_mean4_family_16S_all5_gt0.001.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_sum4 <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_ALL5,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_sum4",
  pdf_name = "clr_mean4_family_18S_all5_gt0.001.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

message("Done (all5, family facets).")
message(glue("CLR  16S: {res_16s_clr$pdf} (n_ASV={res_16s_clr$n_asv})"))
message(glue("CLR  18S: {res_18s_clr$pdf} (n_ASV={res_18s_clr$n_asv})"))
message(glue("rCLR 16S: {res_16s_rclr$pdf}"))
message(glue("rCLR 18S: {res_18s_rclr$pdf}"))
message(glue("mean4 16S: {res_16s_sum4$pdf}"))
message(glue("mean4 18S: {res_18s_sum4$pdf}"))
