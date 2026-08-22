# Functions used by other scripts 
# Libraries required, either download or use the docker container
library(tidyverse)
library(microViz)
library(bubbler)
library(vegan)
library(phyloseq)
library(ggrepel)
library(ggplot2)
library(viridis)
library(glue)
library(ggtext)
library(patchwork)
library(indicspecies)
library(janitor)
library(ANCOMBC)
library(parallel)
library(ggside)
library(magick)
library(ggpp)
library(tidyr)
library(scales)

# CLR transform for taxa-as-rows count/abundance matrix (samples as columns).
# Pseudocount of 1, then compositional (proportions), then CLR:
#   clr(x)_j = log(x_j) - mean(log(x))  for each sample (column).
# Equivalent to microViz::tax_transform("clr", add = 1) /
# microbiome::transform after otu + 1 (verified: identical matrix, PC scores cor = 1).
clr_transform <- function(otu) {
    otu <- as.matrix(otu) + 1
    # relative abundance per sample (column)
    otu <- sweep(otu, 2, colSums(otu), "/")
    # CLR per sample
    log_otu <- log(otu)
    sweep(log_otu, 2, colMeans(log_otu), "-")
}

# Unconstrained RDA (= PCA) on a taxa-as-rows matrix; returns vegan rda object.
# Same as microViz::ord_calc(method = "PCA") which wraps phyloseq::ordinate(..., "RDA").
pca_from_otu <- function(otu) {
    # vegan expects samples as rows
    vegan::rda(t(as.matrix(otu)))
}

# Aggregate to a taxonomic rank; returns taxa-as-rows OTU matrix + phyloseq.
genus_otu <- function(phy, tax_level = "Genus") {
    phy_agg <- phy %>%
        tax_fix() %>%
        tax_agg(rank = tax_level) %>%
        ps_get()
    otu <- as(otu_table(phy_agg), "matrix")
    if (!taxa_are_rows(phy_agg)) {
        otu <- t(otu)
    }
    list(phy = phy_agg, otu = otu)
}

# Aitchison distance = Euclidean distance in CLR space (pseudocount 1).
aitchison_dist <- function(otu) {
    stats::dist(t(clr_transform(otu)))
}

# 2D NMDS from a distance matrix (or sample-rows abundance matrix + method).
nmds_from_dist <- function(d, seed = 123L, trymax = 100L) {
    set.seed(seed)
    vegan::metaMDS(
        d,
        k = 2,
        trymax = trymax,
        autotransform = FALSE,
        trace = FALSE
    )
}

# Classical metric PCoA (cmdscale) from a distance matrix.
pcoa_from_dist <- function(d, k = 2L) {
    ord <- stats::cmdscale(d, k = k, eig = TRUE)
    pos_eig <- ord$eig[ord$eig > 0]
    var_explained <- 100 * ord$eig[seq_len(k)] / sum(pos_eig)
    list(
        points = ord$points,
        eig = ord$eig,
        var_explained = var_explained
    )
}

# Attach sample metadata to site scores (row names = sample ids).
ord_scores_df <- function(points, phy, axis_names = colnames(points)) {
    scores_df <- as.data.frame(points)
    colnames(scores_df) <- axis_names[seq_len(ncol(scores_df))]
    scores_df <- scores_df %>%
        tibble::rownames_to_column("Sample")
    meta <- sample_data(phy) %>%
        as_tibble(rownames = "Sample")
    scores_df %>% left_join(meta, by = "Sample")
}

# Point-only ordination plot styled like full-dataset PCA (Date / Location).
# Optional taxa scores: arrows + labels (from CA / CCA species scores).
plot_ord_points <- function(
    scores,
    x,
    y,
    colour,
    title,
    tag,
    fill_name,
    fill_labels = waiver(),
    xlab = NULL,
    ylab = NULL,
    caption = NULL,
    point_size = 5,
    taxa = NULL,
    taxa_size = 2.6,
    taxa_colour = "grey25"
) {
    colour <- enquo(colour)
    x <- enquo(x)
    y <- enquo(y)

    p <- ggplot(scores, aes(x = !!x, y = !!y, fill = !!colour, shape = Location)) +
        geom_hline(yintercept = 0, colour = "grey35", alpha = 0.65, linetype = "dotted", linewidth = 0.8) +
        geom_vline(xintercept = 0, colour = "grey35", alpha = 0.65, linetype = "dotted", linewidth = 0.8) +
        geom_point(size = point_size, color = "black", stroke = 0.5, alpha = 0.7) +
        scale_shape_manual(
            name = "Location",
            values = c("MS" = 21, "WS" = 24),
            labels = c("MS" = "Microscope slide", "WS" = "Wall strip")
        ) +
        scale_fill_viridis_d(name = fill_name, labels = fill_labels) +
        guides(
            shape = guide_legend(override.aes = list(color = "black", stroke = 0.5)),
            fill = guide_legend(override.aes = list(color = "black", stroke = 0.5, shape = 21))
        ) +
        labs(tag = tag, x = xlab, y = ylab, caption = caption) +
        ggtitle(title) +
        plastic_theme

    if (!is.null(taxa) && nrow(taxa) > 0) {
        p <- p +
            geom_segment(
                data = taxa,
                aes(x = 0, y = 0, xend = !!x, yend = !!y),
                colour = "grey50",
                linewidth = 0.35,
                alpha = 0.55,
                arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"),
                inherit.aes = FALSE
            ) +
            ggrepel::geom_text_repel(
                data = taxa,
                aes(x = !!x, y = !!y, label = Taxon),
                size = taxa_size,
                colour = taxa_colour,
                fontface = "italic",
                max.overlaps = 40,
                min.segment.length = 0,
                box.padding = 0.25,
                point.padding = 0.15,
                segment.colour = "grey60",
                inherit.aes = FALSE
            )
    }
    p
}

save_ord_pdf <- function(plot, path, width = 12 * 0.85, height = 20 * 0.85) {
    pdf(path, width = width, height = height, family = "Helvetica", useDingbats = FALSE)
    print(plot)
    dev.off()
    invisible(path)
}

# 2x2: left = colour by time (Date), right = colour by plastic_level; rows = 16S then 18S.
# Under left: time+site PERMANOVA. Under right: plastic PERMANOVA; stress (if any) on far right.
# Shared legend at the bottom.
combine_time_plastic_4panel <- function(
    p_16s_time,
    p_16s_plastic,
    p_18s_time,
    p_18s_plastic,
    caption_16s_left = NULL,
    caption_16s_right = NULL,
    caption_18s_left = NULL,
    caption_18s_right = NULL,
    stress_16s = NULL,
    stress_18s = NULL
) {
    no_cap <- function(p) {
        p + labs(caption = NULL)
    }

    # Centered text caption under a left panel
    cap_cell <- function(txt) {
        if (is.null(txt) || !nzchar(as.character(txt))) {
            return(patchwork::plot_spacer())
        }
        label <- stringr::str_wrap(as.character(txt), width = 48)
        ggplot() +
            annotate(
                "text",
                x = 0.5,
                y = 0.5,
                label = label,
                size = 3.5,
                hjust = 0.5,
                vjust = 0.5,
                lineheight = 1.1
            ) +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            theme_void() +
            theme(plot.margin = margin(2, 6, 6, 6))
    }

    # Under right panel: plastic PERMANOVA (center-left) + stress on the right edge
    right_cap_cell <- function(plastic_txt, stress = NULL) {
        has_pl <- !is.null(plastic_txt) && nzchar(as.character(plastic_txt))
        has_st <- !is.null(stress) && is.finite(stress)
        if (!has_pl && !has_st) {
            return(patchwork::plot_spacer())
        }

        p <- ggplot() +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            theme_void() +
            theme(plot.margin = margin(2, 10, 6, 6))

        if (has_pl) {
            p <- p +
                annotate(
                    "text",
                    x = 0.35,
                    y = 0.5,
                    label = stringr::str_wrap(as.character(plastic_txt), width = 40),
                    size = 3.5,
                    hjust = 0.5,
                    vjust = 0.5,
                    lineheight = 1.1
                )
        }
        if (has_st) {
            p <- p +
                annotate(
                    "text",
                    x = 1,
                    y = 0.5,
                    label = glue("Stress = {round(stress, 3)}"),
                    size = 3.5,
                    hjust = 1,
                    vjust = 0.5,
                    fontface = "plain"
                )
        }
        p
    }

    pA <- no_cap(p_16s_time)
    pB <- no_cap(p_16s_plastic)
    pC <- no_cap(p_18s_time)
    pD <- no_cap(p_18s_plastic)
    cA <- cap_cell(caption_16s_left)
    cB <- right_cap_cell(caption_16s_right, stress = stress_16s)
    cC <- cap_cell(caption_18s_left)
    cD <- right_cap_cell(caption_18s_right, stress = stress_18s)

    # AB plots / captions under A,B / CD plots / captions under C,D
    (pA + pB + cA + cB + pC + pD + cC + cD) +
        patchwork::plot_layout(
            design = "AB\nCD\nEF\nGH",
            heights = c(12, 1.5, 12, 1.5),
            guides = "collect"
        ) &
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.box.just = "center",
            legend.direction = "horizontal",
            legend.title.align = 0.5
        )
}

# 2x3: time | plastic (full) | plastic (week-9 wall strips). Rows: 16S then 18S.
# Captions under each column; full-dataset stress under middle panel (far right of that cell);
# week-9 subset stress under the third panel.
combine_time_plastic_w9_6panel <- function(
    p_16s_time,
    p_16s_plastic,
    p_16s_w9,
    p_18s_time,
    p_18s_plastic,
    p_18s_w9,
    caption_16s_time = NULL,
    caption_16s_plastic = NULL,
    caption_16s_w9 = NULL,
    caption_18s_time = NULL,
    caption_18s_plastic = NULL,
    caption_18s_w9 = NULL,
    stress_16s = NULL,
    stress_16s_w9 = NULL,
    stress_18s = NULL,
    stress_18s_w9 = NULL
) {
    no_cap <- function(p) p + labs(caption = NULL)

    cap_cell <- function(txt) {
        if (is.null(txt) || !nzchar(as.character(txt))) {
            return(patchwork::plot_spacer())
        }
        label <- stringr::str_wrap(as.character(txt), width = 42)
        ggplot() +
            annotate(
                "text",
                x = 0.5,
                y = 0.5,
                label = label,
                size = 3.2,
                hjust = 0.5,
                vjust = 0.5,
                lineheight = 1.1
            ) +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            theme_void() +
            theme(plot.margin = margin(2, 4, 6, 4))
    }

    right_cap_cell <- function(plastic_txt, stress = NULL) {
        has_pl <- !is.null(plastic_txt) && nzchar(as.character(plastic_txt))
        has_st <- !is.null(stress) && is.finite(stress)
        if (!has_pl && !has_st) {
            return(patchwork::plot_spacer())
        }
        p <- ggplot() +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
            theme_void() +
            theme(plot.margin = margin(2, 8, 6, 4))
        if (has_pl) {
            p <- p +
                annotate(
                    "text",
                    x = 0.3,
                    y = 0.5,
                    label = stringr::str_wrap(as.character(plastic_txt), width = 32),
                    size = 3.2,
                    hjust = 0.5,
                    vjust = 0.5,
                    lineheight = 1.1
                )
        }
        if (has_st) {
            p <- p +
                annotate(
                    "text",
                    x = 1,
                    y = 0.5,
                    label = glue("Stress = {round(stress, 3)}"),
                    size = 3.2,
                    hjust = 1,
                    vjust = 0.5
                )
        }
        p
    }

    pA <- no_cap(p_16s_time)
    pB <- no_cap(p_16s_plastic)
    pC <- no_cap(p_16s_w9)
    pD <- no_cap(p_18s_time)
    pE <- no_cap(p_18s_plastic)
    pF <- no_cap(p_18s_w9)
    cA <- cap_cell(caption_16s_time)
    cB <- right_cap_cell(caption_16s_plastic, stress = stress_16s)
    cC <- right_cap_cell(caption_16s_w9, stress = stress_16s_w9)
    cD <- cap_cell(caption_18s_time)
    cE <- right_cap_cell(caption_18s_plastic, stress = stress_18s)
    cF <- right_cap_cell(caption_18s_w9, stress = stress_18s_w9)

    # ABC plots / captions / DEF plots / captions
    (pA + pB + pC + cA + cB + cC + pD + pE + pF + cD + cE + cF) +
        patchwork::plot_layout(
            design = "ABC\nDEF\nGHI\nJKL",
            heights = c(12, 1.6, 12, 1.6),
            guides = "collect"
        ) &
        theme(
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.box.just = "center",
            legend.direction = "horizontal",
            legend.title.align = 0.5
        )
}

# Align sample metadata to a distance matrix sample order.
permanova_meta_aligned <- function(d, phy) {
    meta <- sample_data(phy) %>%
        as_tibble(rownames = "Sample") %>%
        mutate(
            Date = factor(Date),
            Location = factor(Location),
            plastic_level = factor(plastic_level)
        )
    sample_ids <- attr(d, "Labels")
    if (is.null(sample_ids)) {
        sample_ids <- labels(d)
    }
    if (is.null(sample_ids)) {
        stop("Distance object has no sample labels")
    }
    meta <- meta %>%
        filter(Sample %in% sample_ids) %>%
        arrange(match(Sample, sample_ids))
    stopifnot(identical(meta$Sample, as.character(sample_ids)))
    meta
}

pca_plot0 <- function(phy, colour = NULL, shape = NULL, r2_cutoff = 0.02, tax_level = "Genus", transform = "clr", title = "PCA", point_size = 4, italics = TRUE, tax_lab_size = 3, highlight_samples = NULL, show_vectors = TRUE, time_site_ellipses = FALSE) {
    
    colour <- enquo(colour)
    shape <- enquo(shape)

    # Aggregate with microViz helpers only (tax_fix / tax_agg).
    # Transform + PCA are ours: clr_transform (pseudocount 1) + vegan::rda.
    phy_agg <- phy %>%
        tax_fix() %>%
        tax_agg(rank = tax_level) %>%
        ps_get()

    otu <- as(otu_table(phy_agg), "matrix")
    if (!taxa_are_rows(phy_agg)) {
        otu <- t(otu)
    }

    if (transform == "clr") {
        otu <- clr_transform(otu)
    } else if (transform != "identity") {
        stop("Unsupported transform: ", transform, ". Use 'clr' or 'identity'.")
    }

    mod <- pca_from_otu(otu)

    smp_scores <- scores(mod, display = "sites") %>%
        as_tibble(rownames = "Sample")

    r.names <- rownames(sample_data(phy))
    meta <- sample_data(phy) %>%
        as_tibble() %>%
        bind_cols(Sample = r.names)

    smp_scores <- smp_scores %>%
        left_join(meta, by = "Sample") %>%
        mutate(time_site = interaction(Date, Location, drop = TRUE, sep = "_"))

    spp_scores <- scores(mod, display = "species") %>%
        as_tibble(rownames = "Species")

    spp_scores <- spp_scores %>%
        mutate(
            r2_PCA1 = (PC1^2) / sum(PC1^2),
            r2_PCA2 = (PC2^2) / sum(PC2^2),
            max_r2 = pmax(r2_PCA1, r2_PCA2)
        )

    selected_spp_scores <- spp_scores %>%
        filter(max_r2 >= r2_cutoff) %>%
        mutate(fontface = if_else(italics & !Species %in% c("Leptolyngbyaceae Family", "LKM11", "Glissomonadida Family"), 
                                   "italic", "plain"))

    eig <- eigenvals(mod)
    var_explained <- (eig / sum(eig)) * 100

    hi_scores <- if (!is.null(highlight_samples)) {
        smp_scores %>% filter(Sample %in% highlight_samples)
    } else {
        smp_scores %>% filter(FALSE)
    }

    base_scores <- if (!is.null(highlight_samples)) {
        smp_scores %>% filter(!Sample %in% highlight_samples)
    } else {
        smp_scores
    }

    # Ellipses need >= 3 points per group
    ellipse_data <- smp_scores %>%
        add_count(time_site, name = "n_group") %>%
        filter(n_group >= 3)

    ellipse_labels <- NULL

    p <- ggplot() +
        geom_hline(yintercept = 0, colour = "grey35", alpha = 0.65, linetype = "dotted", linewidth = 0.8) +
        geom_vline(xintercept = 0, colour = "grey35", alpha = 0.65, linetype = "dotted", linewidth = 0.8)

    if (isTRUE(time_site_ellipses) && nrow(ellipse_data) > 0) {
        # Place labels just outside each cloud, pushed away from the plot centre
        cx <- mean(range(smp_scores$PC1))
        cy <- mean(range(smp_scores$PC2))
        ellipse_labels <- ellipse_data %>%
            group_by(time_site, Date, Location) %>%
            summarise(
                x0 = mean(PC1),
                y0 = mean(PC2),
                x_half = 0.5 * diff(range(PC1)),
                y_half = 0.5 * diff(range(PC2)),
                .groups = "drop"
            ) %>%
            mutate(
                dx = x0 - cx,
                dy = y0 - cy,
                len = sqrt(dx^2 + dy^2) + 1e-6,
                ux = dx / len,
                uy = dy / len,
                pad = 1.35,
                # Segment starts at outer edge of cloud; label sits further out
                x0 = x0 + ux * (x_half * 1.15),
                y0 = y0 + uy * (y_half * 1.15),
                PC1 = x0 + ux * pad,
                PC2 = y0 + uy * pad,
                ellipse_label = case_when(
                    Location == "WS" ~ paste0("Wall strip week ", as.character(Date)),
                    Location == "MS" ~ paste0("Microscope slide week ", as.character(Date)),
                    TRUE ~ as.character(time_site)
                )
            )

        p <- p +
            stat_ellipse(
                data = ellipse_data,
                aes(x = PC1, y = PC2, group = time_site),
                type = "t",
                level = 0.95,
                geom = "polygon",
                alpha = 0.10,
                fill = "grey35",
                colour = "grey55",
                linewidth = 0.35,
                linetype = "dotted",
                show.legend = FALSE
            ) +
            geom_segment(
                data = ellipse_labels,
                aes(x = x0, y = y0, xend = PC1, yend = PC2),
                colour = "grey55",
                linewidth = 0.35,
                alpha = 0.85,
                show.legend = FALSE
            ) +
            geom_label(
                data = ellipse_labels,
                aes(x = PC1, y = PC2, label = ellipse_label),
                fill = alpha("white", 0.92),
                colour = "grey25",
                size = 3.4,
                fontface = "plain",
                linewidth = 0.2,
                label.padding = unit(0.18, "lines"),
                label.r = unit(0.1, "lines"),
                show.legend = FALSE
            )
    }

    p <- p +
        geom_point(
            data = base_scores,
            aes(x = PC1, y = PC2, fill = !!colour, shape = !!shape),
            size = point_size,
            alpha = 0.7,
            color = "black",
            stroke = 0.5
        ) +
        geom_point(
            data = hi_scores,
            aes(x = PC1, y = PC2, fill = !!colour, shape = !!shape),
            size = point_size,
            alpha = 0.9,
            color = "#D62728",
            stroke = 1.8,
            show.legend = FALSE
        )

    if (isTRUE(show_vectors) && nrow(selected_spp_scores) > 0) {
        p <- p +
            geom_segment(data = selected_spp_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                         arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7, colour = "black") +
            geom_text_repel(data = selected_spp_scores, aes(x = PC1, y = PC2, label = Species, fontface = fontface),
                size = tax_lab_size,
                alpha = 1,
                position = position_nudge_center(0.2, 0.1, 0, 0),
                box.padding = 1.5,
                point.padding = 1.5,
                min.segment.length = 0.3,
                force = 2,
                max.overlaps = Inf,
                max.iter = 15000)
    }

    # Expand limits so outward ellipse labels stay inside the panel
    if (!is.null(ellipse_labels) && nrow(ellipse_labels) > 0) {
        xr <- range(c(smp_scores$PC1, ellipse_labels$PC1), na.rm = TRUE)
        yr <- range(c(smp_scores$PC2, ellipse_labels$PC2), na.rm = TRUE)
        # Extra room for label box width/height
        x_pad <- diff(xr) * 0.12 + 0.35
        y_pad <- diff(yr) * 0.12 + 0.35
        p <- p +
            coord_cartesian(
                xlim = c(xr[1] - x_pad, xr[2] + x_pad),
                ylim = c(yr[1] - y_pad, yr[2] + y_pad),
                clip = "off"
            )
    } else {
        p <- p + coord_cartesian(clip = "off")
    }

    p +
        labs(x = glue("PC1 [{round(var_explained[1], 1)}%]"), 
            y = glue("PC2 [{round(var_explained[2], 1)}%]")) +
        ggtitle(glue("{title}"))

}

vpa_plot <- function(phy, v1, v2, v3 = NULL, label1 = NULL, label2 = NULL, label3 = NULL, title = NULL, output = "plot") {

    v1 <- enquo(v1)
    v2 <- enquo(v2)
    v3 <- if(!is.null(v3)) quo_name(enquo(v3)) else NULL

    # data-scecific step
    phy2 <- phy %>%
        ps_mutate(week_3 = if_else(Date == "1", true = 1, false = 0),
                  week_6 = if_else(Date == "2", true = 1, false = 0), 
                  week_9 = if_else(Date == "3", true = 1, false = 0),
                  wall_strip = if_else(Location == "WS", true = 1, false = 0)) %>%
        tax_fix() %>%
        tax_transform("clr", rank = "Genus")
        
    otu_table <- as(otu_table(phy2), "matrix")
    otu_table <- t(otu_table)
    meta_data <- as(sam_data(phy2), "data.frame")
    env1 <- meta_data %>% dplyr::select(!!v1) %>% dplyr::pull()
    env2 <- meta_data %>% dplyr::select(!!v2) %>% dplyr::pull()

    if(!is.null(v3)) {
        env3 <- meta_data[[v3]]
        vp <- varpart(otu_table, env1, env2, env3)
        # Round all values to 2 decimal places
        vp$part$indfract <- round(vp$part$indfract, 2)
        vp$part$fract <- round(vp$part$fract, 2)
        
        if(output == "plot") {
            # Create a custom plot
            plot(vp, Xnames = c(
                if(!is.null(label1)) label1 else "X1",
                if(!is.null(label2)) label2 else "X2",
                if(!is.null(label3)) label3 else "X3"
            ), main = title)
        } else {
            return(vp$part)
        }
    } else {
        vp <- varpart(otu_table, env1, env2)
        # Round all values to 2 decimal places
        vp$part$indfract <- round(vp$part$indfract, 2)
        vp$part$fract <- round(vp$part$fract, 2)
        
        if(output == "plot") {
            # Create a custom plot
            plot(vp, Xnames = c(
                if(!is.null(label1)) label1 else "X1",
                if(!is.null(label2)) label2 else "X2"
            ), main = title)
        } else {
            return(vp$part)
        }
    }
}

create_ancombc_plot <- function(output_data, variable_name, title = NULL) {
  # Extract variable column names
  var_cols <- grep(variable_name, names(output_data$res), value = TRUE)
  
  # Select data
  df_var <- output_data$res %>%
    dplyr::select(taxon, all_of(var_cols))
  
  # Get column name prefixes
  lfc_col <- paste0("lfc_", variable_name)
  se_col <- paste0("se_", variable_name)
  diff_col <- paste0("diff_", variable_name) 
  passed_ss_col <- paste0("passed_ss_", variable_name)
  
  # Create figure data
  df_fig_var <- output_data$res %>%
    dplyr::filter(!!sym(diff_col) == 1) %>%
    arrange(desc(!!sym(lfc_col))) %>%
    dplyr::mutate(
      direct = ifelse(!!sym(lfc_col) > 0, "Positive LFC", "Negative LFC"),
      color = ifelse(!!sym(passed_ss_col) == 1, "#000000", "#8a8a8a")
    )
  
  # Set factors
  df_fig_var$taxon <- factor(df_fig_var$taxon, levels = df_fig_var$taxon)
  df_fig_var$direct <- factor(df_fig_var$direct, levels = c("Positive LFC", "Negative LFC"))
  
  # Set default title if not provided
  if(is.null(title)) {
    title <- paste0("LFC as one unit increase of ", variable_name)
  }
  
  # Create plot
  p <- df_fig_var %>%
    ggplot(aes(x = taxon, y = !!sym(lfc_col), fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(
      aes(ymin = !!sym(lfc_col) - !!sym(se_col), 
          ymax = !!sym(lfc_col) + !!sym(se_col)), 
      width = 0.2, position = position_dodge(0.05), color = "black"
    ) + 
    labs(x = NULL, y = "Log fold change", title = title) + 
    scale_fill_viridis_d(name = NULL, option = "D", direction = -1) +
    scale_color_viridis_d(name = NULL, option = "D", direction = -1) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1,
                                     color = df_fig_var$color))
  
  return(list(plot = p, data = df_fig_var, full_data = df_var))
}

plot_stacked_barchart <- function(gene = "16S", taxa_level = "Genus", n_taxa = 20, italics = TRUE, position = "fill", tag = "A", fill = "D", x_lab = NULL) {

    phy_obj <- get(glue("phy_{gene}"))
    rel_abund <- rel_abund_phy(phy_obj, 
                               meta_data = TRUE, 
                               taxa_level = taxa_level ) %>% 
            {if(italics) taxon_italics(.) else .} %>%
            pool_taxa(n_taxa = n_taxa, keep_metadata = TRUE) %>%
            arrange_taxa() %>%
            mutate(sample_id = factor(sample_id)) %>%
            mutate(sample_id = fct_reorder(sample_id, plastic_concentration, .desc = F))
            # mutate(Date = factor(Date, levels = c(3, 6, 9)))

    unique_taxa <- rel_abund %>%
        add_other() %>%
        all_taxa()

    rel_abund_ws <- rel_abund %>%
        dplyr::filter(str_starts(sample_id, c("WS")))

    rel_abund_ms <- rel_abund %>%
        dplyr::filter(str_starts(sample_id, c("MS"))) 

    sub_unique_taxa <- intersect(rel_abund_ws[["taxon"]], rel_abund_ms[["taxon"]])
    colorscheme <- global_colour_scheme(unique_taxa, sub_unique_taxa)

    new_layer <- sum_rel_abund(rel_abund, "sample_id") 

    p1 <-inner_join(rel_abund_ws, new_layer, by = "sample_id") %>%
        bar_plot(
            x_var = "sample_id", 
            position = position,
            global_colours = colorscheme, 
            italics = italics) + 
            scale_fill_viridis_d(option = fill) +
        ggside::geom_xsidepoint(aes(x = sample_id, y = sum, group = 1), show.legend = FALSE) +
        ggside::geom_xsideline(aes(x = sample_id, y = sum, group = 1), show.legend = FALSE) +
        ggside::scale_xsidey_continuous(limits = c(0, 0.15)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        guides(fill=guide_legend(title = glue("{taxa_level}")))  +
        scale_y_continuous(expand = expansion(mult = c(0, .05)), name = "Relative abundance")  +
        facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
            Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x)),
            Date = function(x) paste("Week", x)
        )) +
        theme(strip.text = element_text(size = 12),
              axis.text.x = element_text(size = 10)) +
        labs(x = x_lab, tag = tag) +
        theme(plot.tag = element_text(face = "bold", size = 28))

    p2 <- inner_join(rel_abund_ms, new_layer, by = "sample_id") %>%
        bar_plot(
            x_var = "sample_id", 
            position = position,
            global_colours = colorscheme, 
            italics = italics) + 
            scale_fill_viridis_d(option = fill) +
        ggside::geom_xsidepoint(aes(x = sample_id, y = sum, group = 1 ), show.legend = FALSE) +
        ggside::geom_xsideline(aes(x = sample_id, y = sum, group = 1 ), show.legend = FALSE) +
        ggside::scale_xsidey_continuous(limits = c(0, 0.15)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()) + 
        guides(fill=guide_legend(title = glue("{taxa_level}")))  +
        scale_y_continuous(expand = expansion(mult = c(0, .05)), name = NULL)  +
        facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
            Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x)),
            Date = function(x) paste("Week", x)
        )) +
        theme(strip.text = element_text(size = 12),
              axis.text.x = element_text(size = 10)) +
        labs(x = x_lab)

    # Figure 2.A 
    layout <- "AAAB"
    p1 + p2 + patchwork::plot_layout(guides = "collect", axes = "collect", design = layout)
}

prepare_ancombc_plot_heatmap <- function(output_data) {

    res_trend <- output_data$res_trend

    # Analyze significant taxa
    significant_taxa <- res_trend %>%
        dplyr::filter(diff_abn == TRUE) %>%
        dplyr::select(taxon, 
                    lfc_plastic_levellow,
                    lfc_plastic_levelmedium,
                    lfc_plastic_levelhigh,
                    W, p_val, q_val, passed_ss) %>%
        dplyr::arrange(desc(W))

    # Create a heatmap of log fold changes
    # Convert to long format for ggplot
    df_long <- significant_taxa %>%
        tidyr::pivot_longer(
            cols = c(lfc_plastic_levellow, lfc_plastic_levelmedium, lfc_plastic_levelhigh),
            names_to = "plastic_level",
            values_to = "lfc"
        ) %>%
        dplyr::mutate(
            plastic_level = factor(
                plastic_level,
                levels = c("lfc_plastic_levellow", "lfc_plastic_levelmedium", "lfc_plastic_levelhigh"),
                labels = c("Low", "Medium", "High")
            )
        )

    # Order: High ascending, then refine adjacent pairs using Low / Medium when similar
    order_taxa <- order_taxa_plastic_lfc_columns(
        df_long,
        id_col = "taxon",
        level_col = "plastic_level",
        value_col = "lfc",
        col_order = c("Low", "Medium", "High")
    )

    df_long %>%
        mutate(taxon = factor(taxon, levels = order_taxa))
}

# Order taxa for plastic-level heatmaps.
# 1) Sort by far-right column (High) ascending (purple → yellow / low → high LFC).
# 2) For each adjacent pair, walk Low → Medium: if values are similar (within tol)
#    on the current column, open the next column; if that column differs,
#    put the larger LFC lower on the heatmap. Skip High for this pair cascade
#    (it is already the primary sort).
order_taxa_plastic_lfc_columns <- function(
    long_df,
    id_col = "taxon_display",
    level_col = "plastic_level",
    value_col = "lfc",
    col_order = c("Low", "Medium", "High"),
    tol = 0.15
) {
    wide <- long_df %>%
        dplyr::select(dplyr::all_of(c(id_col, level_col, value_col))) %>%
        tidyr::pivot_wider(
            names_from = dplyr::all_of(level_col),
            values_from = dplyr::all_of(value_col)
        )

    missing_cols <- setdiff(col_order, names(wide))
    if (length(missing_cols)) {
        stop("Missing plastic-level columns for ordering: ", paste(missing_cols, collapse = ", "))
    }

    # Primary sort: High ascending (ggplot factor y: first level at bottom)
    high_col <- col_order[length(col_order)]
    o <- order(wide[[high_col]], na.last = TRUE)
    mat <- as.matrix(wide[o, col_order, drop = FALSE])
    ids <- wide[[id_col]][o]
    n <- nrow(mat)
    if (n <= 1L) {
        return(ids)
    }

    # Columns used in the similarity cascade (all but High)
    n_refine <- length(col_order) - 1L
    if (n_refine < 1L) {
        return(ids)
    }

    changed <- TRUE
    max_iter <- n * n
    iter <- 0L
    while (isTRUE(changed) && iter < max_iter) {
        changed <- FALSE
        iter <- iter + 1L
        for (i in seq_len(n - 1L)) {
            for (j in seq_len(n_refine)) {
                a <- mat[i, j]
                b <- mat[i + 1L, j]
                if (is.na(a) || is.na(b)) {
                    break
                }
                if (abs(a - b) > tol) {
                    # Not similar: keep High order
                    break
                }
                # Similar on column j → compare next column (before High)
                next_j <- j + 1L
                if (next_j > n_refine) {
                    # Next would be High: skip last column
                    break
                }
                a2 <- mat[i, next_j]
                b2 <- mat[i + 1L, next_j]
                if (is.na(a2) || is.na(b2)) {
                    break
                }
                if (abs(a2 - b2) <= tol) {
                    # Still similar: continue cascade
                    next
                }
                # Differ on next column: larger value belongs lower on plot
                # (earlier in factor levels = bottom = position i)
                if (a2 < b2) {
                    mat[c(i, i + 1L), ] <- mat[c(i + 1L, i), ]
                    ids[c(i, i + 1L)] <- ids[c(i + 1L, i)]
                    changed <- TRUE
                }
                break
            }
        }
    }
    ids
}

plastic_theme <- theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    plot.margin = margin(20, 60, 20, 20),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title.align = 0.5,
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.tag = element_text(face = "bold", size = 28)
  )

plastic_theme_2 <- theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.margin = margin(20, 60, 20, 20),
    panel.grid.minor = element_blank(),
    # legend.position = "bottom",
    legend.title.align = 0.5,
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.tag = element_text(face = "bold", size = 24)
  )
