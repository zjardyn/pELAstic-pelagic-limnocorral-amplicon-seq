source("R/01_load_files2.R")
source("R/functions.R")

# theme_set(theme_linedraw(base_size = 12))
theme_set(theme_bw(base_size = 14))

# Drop samples with library size below this (total sum of taxa counts). Set to 0 to keep all.
# 5000 drops the two lowest 16S libs (~4.7–4.8k); 18S also drops WS_A_6 (failed) and WS_G_6 (~4.1k).
min_sum_taxa <- 5000L

# Sum of taxa abundances per sample (library size)
taxa_sums_per_sample <- function(phy, label) {
  tibble(
    dataset = label,
    sample = sample_names(phy),
    sum_taxa = as.numeric(sample_sums(phy))
  ) %>%
    arrange(sample)
}

prune_samples_min_depth <- function(phy, label, min_reads) {
  if (min_reads <= 0) {
    return(phy)
  }
  ss <- sample_sums(phy)
  keep <- unname(ss) >= min_reads
  dropped <- names(ss)[!keep]
  if (length(dropped)) {
    message(
      label, ": removing ", length(dropped), " sample(s) with sum_taxa < ", min_reads, ": ",
      paste(dropped, collapse = ", ")
    )
  } else {
    message(label, ": no samples below sum_taxa = ", min_reads)
  }
  phy <- prune_samples(keep, phy)
  prune_taxa(taxa_sums(phy) > 0, phy)
}

message("Taxa sums per sample (before depth filter)")
print(taxa_sums_per_sample(phy_16S, "16S"), n = Inf)
print(taxa_sums_per_sample(phy_18S, "18S"), n = Inf)

phy_16S <- prune_samples_min_depth(phy_16S, "16S", min_sum_taxa)
phy_18S <- prune_samples_min_depth(phy_18S, "18S", min_sum_taxa)

message("Taxa sums per sample (after depth filter)")
print(taxa_sums_per_sample(phy_16S, "16S"), n = Inf)
print(taxa_sums_per_sample(phy_18S, "18S"), n = Inf)

# Rarefy to even depth (reads per sample). Set equal to min_sum_taxa, or lower, so all samples qualify.
rarefy_depth <- min_sum_taxa
if (rarefy_depth > 0L) {
  rarefy_one <- function(phy, label) {
    m <- min(sample_sums(phy))
    if (m < rarefy_depth) {
      stop(label, ": min sum_taxa (", m, ") is below rarefy_depth (", rarefy_depth, ")")
    }
    rarefy_even_depth(
      phy,
      sample.size = rarefy_depth,
      rngseed = 1L,
      replace = FALSE, # phyloseq default is TRUE; FALSE = classical rarefaction (subsample w/o replacement)
      trimOTUs = TRUE,
      verbose = FALSE
    )
  }
  message("Rarefying to ", rarefy_depth, " reads per sample")
  phy_16S <- rarefy_one(phy_16S, "16S")
  phy_18S <- rarefy_one(phy_18S, "18S")
  message("Taxa sums per sample (after rarefaction)")
  print(taxa_sums_per_sample(phy_16S, "16S"), n = Inf)
  print(taxa_sums_per_sample(phy_18S, "18S"), n = Inf)
}

# 16S
meta <- sample_data(phy_16S)
r.names <- row.names(meta)
meta <- meta %>%
    as_tibble() %>%
    bind_cols(sample = r.names)

# Berger-Parker = max taxon abundance / total (same taxonomic rank as Observed/Shannon)
berger_parker <- function(x) {
  max(x) / sum(x)
}

berger_parker_from_phy <- function(phy) {
  otu <- as(otu_table(phy), "matrix")
  if (!taxa_are_rows(phy)) {
    otu <- t(otu)
  }
  apply(otu, 2, berger_parker)
}

# Genus table shared by Observed, Shannon, and Berger-Parker
phy_16S_genus <- phy_16S %>%
    phyloseq_validate(remove_undetected = TRUE) %>%
    tax_fix() %>%
    tax_transform("identity", "Genus") %>%
    ps_arrange(Date)

bp_16S <- data.frame(
    sample = sample_names(phy_16S_genus),
    Berger_Parker = unname(berger_parker_from_phy(phy_16S_genus)),
    stringsAsFactors = FALSE
)

alpha <- phy_16S_genus %>%
    phyloseq::estimate_richness(measures = c("Observed", "Shannon")) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    mutate(sample = str_replace_all(sample, "\\.", "-")) %>%
    left_join(meta, by = "sample") %>%
    left_join(bp_16S, by = "sample") %>%
    as_tibble() %>%
    # Arrange sample by plastic_concentration (same as stacked bar charts)
    mutate(sample = factor(sample)) %>%
    mutate(Location = factor(Location, levels = c("WS", "MS"))) %>%
        #    Date = factor(Date, levels = c(3, 6, 9))) %>%
    # Reorder samples by plastic concentration within each Location
    group_by(Location) %>%
    mutate(sample = fct_reorder(sample, plastic_concentration, .desc = F)) %>%
    ungroup() %>%
    # Create separate factor levels for each location to ensure proper ordering
    group_by(Location) %>%
    mutate(sample = factor(sample, levels = unique(sample[order(plastic_concentration)]))) %>%
    ungroup()

# Calculate consistent y-axis limits
y_range <- range(alpha$Observed, na.rm = TRUE)

# Create separate plots for WS and MS
p1_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = Observed)) + 
    geom_point() + 
    facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
    )) +
                theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          plot.tag = element_text(face = "bold")) +
    ylim(y_range) +
    labs(tag = "A", x = NULL)

p1_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = Observed)) + 
    geom_point() + 
    facet_grid(Location ~ Date, labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
    )) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range) +
    labs(x = NULL)

# Combine the plots with AAAB layout
p1 <- p1_ws + p1_ms + p1_ms + p1_ms + 
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = "16S: Observed species counts")

# Calculate consistent y-axis limits for p2 and p3
y_range_p2 <- range(alpha$Berger_Parker, na.rm = TRUE)
y_range_p3 <- range(alpha$Shannon, na.rm = TRUE)

# Create separate plots for p2 (Berger-Parker)
p2_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = Berger_Parker)) + 
    geom_point() + 
    facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
    )) +
                theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          plot.tag = element_text(face = "bold")) +
    ylim(y_range_p2) +
    labs(tag = "B", x = NULL)

p2_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = Berger_Parker)) + 
    geom_point() + 
    facet_grid(Location ~ Date, labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
    )) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range_p2) +
    labs(x = NULL)

# Create separate plots for p3 (Shannon)
p3_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = Shannon)) + 
    geom_point() + 
    facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
    )) +
                theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          plot.tag = element_text(face = "bold")) +
    ylim(y_range_p3) +
    labs(tag = "C", x = "Sample name")

p3_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = Shannon)) + 
    geom_point() + 
    facet_grid(Location ~ Date, labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
    )) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range_p3) +
    labs(x = "Sample name")

# Combine the plots with AAAB layout
p2 <- p2_ws + p2_ms + p2_ms + p2_ms + 
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = "16S: Berger-Parker index")

p3 <- p3_ws + p3_ms + p3_ms + p3_ms + 
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = "16S: Shannon diversity index")

p_comb <- wrap_plots(p1, p2, p3, ncol = 1, axes = "collect")
ggsave("figures/fig_S3_alpha_diversity_16S.png", p_comb, width = 10, height = 12, dpi = 300, scale = 0.9)

# 18S
meta <- sample_data(phy_18S)
r.names <- row.names(meta)
meta <- meta %>%
    as_tibble() %>%
    bind_cols(sample = r.names)

# Genus table shared by Observed, Shannon, and Berger-Parker
phy_18S_genus <- phy_18S %>%
    phyloseq_validate(remove_undetected = TRUE) %>%
    tax_fix() %>%
    tax_transform("identity", "Genus") %>%
    ps_arrange(Date)

bp_18S <- data.frame(
    sample = sample_names(phy_18S_genus),
    Berger_Parker = unname(berger_parker_from_phy(phy_18S_genus)),
    stringsAsFactors = FALSE
)

alpha <- phy_18S_genus %>%
    phyloseq::estimate_richness(measures = c("Observed", "Shannon")) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    mutate(sample = str_replace_all(sample, "\\.", "-")) %>%
    left_join(meta, by = "sample") %>%
    left_join(bp_18S, by = "sample") %>%
    as_tibble() %>%
    # Arrange sample by plastic_concentration (same as stacked bar charts)
    mutate(sample = factor(sample)) %>%
    mutate(Location = factor(Location, levels = c("WS", "MS"))) %>%
        #    Date = factor(Date, levels = c(3, 6, 9))) %>%
    # Reorder samples by plastic concentration within each Location
    group_by(Location) %>%
    mutate(sample = fct_reorder(sample, plastic_concentration, .desc = F)) %>%
    ungroup() %>%
    # Create separate factor levels for each location to ensure proper ordering
    group_by(Location) %>%
    mutate(sample = factor(sample, levels = unique(sample[order(plastic_concentration)]))) %>%
    ungroup()

# Calculate consistent y-axis limits
y_range <- range(alpha$Observed, na.rm = TRUE)

# Create separate plots for WS and MS
p1_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = Observed)) + 
    geom_point() + 
    facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
    )) +
                theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          plot.tag = element_text(face = "bold")) +
    ylim(y_range) +
    labs(tag = "A", x = NULL )

p1_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = Observed)) + 
    geom_point() + 
    facet_grid(Location ~ Date, labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
    )) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range) +
    labs(x = NULL)

# Combine the plots with AAAB layout
p1 <- p1_ws + p1_ms + p1_ms + p1_ms + 
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = "18S: Observed species counts")

# Calculate consistent y-axis limits for p2 and p3
y_range_p2 <- range(alpha$Berger_Parker, na.rm = TRUE)
y_range_p3 <- range(alpha$Shannon, na.rm = TRUE)

# Create separate plots for p2 (Berger-Parker)
p2_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = Berger_Parker)) + 
    geom_point() + 
    facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
    )) +
                theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          plot.tag = element_text(face = "bold")) +
    ylim(y_range_p2) +
    labs(tag = "B", x = NULL)

p2_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = Berger_Parker)) + 
    geom_point() + 
    facet_grid(Location ~ Date, labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
    )) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range_p2) +
    labs(x = NULL)

# Create separate plots for p3 (Shannon)
p3_ws <- alpha %>%
    filter(Location == "WS") %>%
    ggplot(aes(x = sample, y = Shannon)) + 
    geom_point() + 
    facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
        Location = function(x) "Wall strip",
        Date = function(x) paste("Week", x)
    )) +
                theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
          plot.tag = element_text(face = "bold")) +
    ylim(y_range_p3) +
    labs(tag = "C", x = "Sample name")

p3_ms <- alpha %>%
    filter(Location == "MS") %>%
    ggplot(aes(x = sample, y = Shannon)) + 
    geom_point() + 
    facet_grid(Location ~ Date, labeller = labeller(
        Location = function(x) "Microscope slide",
        Date = function(x) paste("Week", x)
    )) +
    theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
    ylim(y_range_p3) +
    labs(x = "Sample name")

# Combine the plots with AAAB layout
p2 <- p2_ws + p2_ms + p2_ms + p2_ms + 
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = "18S: Berger-Parker index")

p3 <- p3_ws + p3_ms + p3_ms + p3_ms + 
    plot_layout(design = "AAAB", guides = "collect", axes = "collect") +
    plot_annotation(title = "18S: Shannon diversity index")

p_comb <- wrap_plots(p1, p2, p3, ncol = 1, axes = "collect")
ggsave("figures/fig_S4_alpha_diversity_18S.png", p_comb, width = 10, height = 12, dpi = 300, scale = 0.9)

