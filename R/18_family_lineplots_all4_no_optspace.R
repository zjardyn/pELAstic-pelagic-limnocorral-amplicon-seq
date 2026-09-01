# Family line plots for ASVs with LOO mean importance > 0.001 in ALL 4 of
#   clr_0.1, clr_0.5, clr_1, rclr  (excludes rclr_optspace).
# Does not change R/16 (ge3of4 / all5) outputs — writes new PDFs only.
#
# Run: Rscript R/18_family_lineplots_all4_no_optspace.R

lines <- readLines("R/16_rclr_family_lineplots.R")
cut <- which(grepl("^# ---- CLR \\(pc=1\\)", lines))[1]
if (is.na(cut) || cut < 2L) {
  stop("Could not find CLR run block in R/16_rclr_family_lineplots.R")
}
eval(parse(text = lines[seq_len(cut - 1L)]), envir = environment())

TRANSFORM_CORE4 <- c("clr_0.1", "clr_0.5", "clr_1", "rclr")
MIN_TRANSFORMS <- length(TRANSFORM_CORE4)

message(glue(
  "Subset: mean_loo_imp > {IMP_THRESH} in {MIN_TRANSFORMS}/{length(TRANSFORM_CORE4)} ",
  "(clr_0.1, clr_0.5, clr_1, rclr; no optspace)"
))

res_16s_clr <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_1",
  pdf_name = "clr1_family_16S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_clr <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_1",
  pdf_name = "clr1_family_18S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

res_16s_rclr <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "rclr",
  pdf_name = "rclr_family_16S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_rclr <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "rclr",
  pdf_name = "rclr_family_18S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

res_16s_sum4 <- run_marker(
  marker = "16S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_sum4",
  pdf_name = "clr_mean4_family_16S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 4L,
  nrow_pdf = 3L
)

res_18s_sum4 <- run_marker(
  marker = "18S",
  transforms = TRANSFORM_CORE4,
  min_transforms = MIN_TRANSFORMS,
  method = "clr_sum4",
  pdf_name = "clr_mean4_family_18S_all4_gt0.001_no_optspace.pdf",
  ncol_pdf = 3L,
  nrow_pdf = 3L
)

message("Done (all4, no optspace).")
message(glue("CLR  16S: {res_16s_clr$pdf} (n_ASV={res_16s_clr$n_asv})"))
message(glue("CLR  18S: {res_18s_clr$pdf} (n_ASV={res_18s_clr$n_asv})"))
message(glue("rCLR 16S: {res_16s_rclr$pdf}"))
message(glue("rCLR 18S: {res_18s_rclr$pdf}"))
message(glue("mean4 16S: {res_16s_sum4$pdf}"))
message(glue("mean4 18S: {res_18s_sum4$pdf}"))
