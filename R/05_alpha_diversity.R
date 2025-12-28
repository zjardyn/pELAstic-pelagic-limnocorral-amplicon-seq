source("R/01_load_files.R")
source("R/functions.R")

# theme_set(theme_linedraw(base_size = 12))
theme_set(theme_bw(base_size = 14))

# 16S
meta <- sample_data(phy_16S)
r.names <- row.names(meta)
meta <- meta %>%
    as_tibble() %>%
    bind_cols(sample = r.names)

# Calculate Berger-Parker index
berger_parker <- function(x) {
  max(x) / sum(x)
}

# Calculate Berger-Parker for 16S
bp_16S <- data.frame(
    sample = sample_names(phy_16S),
    Berger_Parker = apply(otu_table(phy_16S), 2, berger_parker)
)

# alpha <- phy_16S %>%
alpha <- phy_16S %>%
    phyloseq_validate(remove_undetected = TRUE) %>%
    tax_fix() %>%
    tax_transform("identity", "Genus") %>%
    ps_arrange(Date) %>%
    phyloseq::estimate_richness(measures = c("Observed", "Shannon")) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    mutate(sample = str_replace_all(sample, "\\.", "-")) %>%
    left_join(meta, by = "sample") %>%
    left_join(bp_16S, by = "sample") %>%
    as_tibble() %>%
    # Arrange sample by plastic_concentration (same as stacked bar charts)
    mutate(sample = factor(sample)) %>%
    mutate(Location = factor(Location, levels = c("WS", "MS")),
           Date = factor(Date, levels = c(3, 6, 9))) %>%
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

# Calculate Berger-Parker for 18S
bp_18S <- data.frame(
    sample = sample_names(phy_18S),
    Berger_Parker = apply(otu_table(phy_18S), 2, berger_parker)
)

alpha <- phy_18S %>%
    phyloseq_validate(remove_undetected = TRUE) %>%
    tax_fix() %>%
    tax_transform("identity", "Genus") %>%
    ps_arrange(Date) %>%
    phyloseq::estimate_richness(measures = c("Observed", "Shannon")) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    mutate(sample = str_replace_all(sample, "\\.", "-")) %>%
    left_join(meta, by = "sample") %>%
    left_join(bp_18S, by = "sample") %>%
    as_tibble() %>%
    # Arrange sample by plastic_concentration (same as stacked bar charts)
    mutate(sample = factor(sample)) %>%
    mutate(Location = factor(Location, levels = c("WS", "MS")),
           Date = factor(Date, levels = c(3, 6, 9))) %>%
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

