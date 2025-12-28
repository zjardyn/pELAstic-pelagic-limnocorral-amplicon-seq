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

pca_plot0 <- function(phy, colour = NULL, shape = NULL, r2_cutoff = 0.02, tax_level = "Genus", transform = "clr", title = "PCA", point_size = 3, italics = TRUE) {
    
    colour <- enquo(colour)
    shape <- enquo(shape)

    if(transform == "hellinger") {
       tax_table(phy) <- decostand(tax_table(phy), method = "hellinger", MARGIN = 2) 
       mod <- phy %>%
            tax_transform("identity", rank = "Genus") %>%
            ord_calc(method = "PCA") %>%
            ord_get()
    } else {
        mod <- phy %>%
            tax_fix() %>%
            tax_transform(transform, rank = tax_level) %>%
            # tax_transform("identity", rank = "Genus") %>%
            ord_calc(method = "PCA") %>%
            ord_get()
    }

    # mod <- phy %>% tax_fix() %>%
    #     tax_transform(transform, rank = tax_level) %>%
    #     # tax_transform("identity", rank = "Genus") %>%
    #     ord_calc(method = "PCA") %>%
    #     ord_get()

    smp_scores <- scores(mod, display = "sites") %>%
        as_tibble(rownames = "Sample")

    r.names <- rownames(sample_data(phy))
    meta <-sample_data(phy) %>%
        as_tibble() %>%
        bind_cols(Sample = r.names)

    smp_scores <- smp_scores %>%
        left_join(meta, by = "Sample")

    spp_scores <- scores(mod, display = "species") %>%
        as_tibble(rownames = "Species")

    spp_scores <- spp_scores %>%
        mutate(
            r2_PCA1 = (PC1^2) / sum(PC1^2),
            r2_PCA2 = (PC2^2) / sum(PC2^2),
            max_r2 = pmax(r2_PCA1, r2_PCA2)  # Take the highest R² value for filtering
        )

    # Select species with R² above the cutoff
    selected_spp_scores <- spp_scores %>%
        filter(max_r2 >= r2_cutoff)

    # For all PCs at once
    eig <- eigenvals(mod)
    var_explained <- (eig / sum(eig)) * 100

    ggplot() + 
        geom_point(data = smp_scores, aes(x = PC1, y = PC2, color = !!colour, shape = !!shape), size = point_size) +
        geom_segment(data = selected_spp_scores, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = unit(0.2, "cm")), alpha = 0.9) + 
        geom_text_repel(data = selected_spp_scores, aes(x = PC1, y = PC2, label = Species),
        size = 3,
        alpha = 0.9, 
        fontface = if(italics) "italic" else "plain") + 
        labs(x = glue("PC1 [{round(var_explained[1], 1)}%]"), 
            y = glue("PC2 [{round(var_explained[2], 1)}%]")) +
        # scale_colour_viridis_d() + 
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

plot_stacked_barchart <- function(gene = "16S", taxa_level = "Genus", n_taxa = 20, italics = TRUE, position = "fill", tag = "A") {

    phy_obj <- get(glue("phy_{gene}"))
    rel_abund <- rel_abund_phy(phy_obj, 
                               meta_data = TRUE, 
                               taxa_level = taxa_level ) %>% 
            {if(italics) taxon_italics(.) else .} %>%
            pool_taxa(n_taxa = n_taxa, keep_metadata = TRUE) %>%
            arrange_taxa() %>%
            mutate(sample_id = factor(sample_id)) %>%
            mutate(sample_id = fct_reorder(sample_id, plastic_concentration, .desc = F)) %>%
            mutate(Date = factor(Date, levels = c(3, 6, 9)))

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
        ggside::geom_xsidepoint(aes(x = sample_id, y = sum, group = 1), show.legend = FALSE) +
        ggside::geom_xsideline(aes(x = sample_id, y = sum, group = 1), show.legend = FALSE) +
        ggside::scale_xsidey_continuous(limits = c(0, 0.15)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        guides(fill=guide_legend(title = glue("{taxa_level}")))  +
        scale_y_continuous(expand = expansion(mult = c(0, .05)), name = "Relative Abundance")  +
        facet_grid(Location ~ Date, scales = "free_x", labeller = labeller(
            Location = function(x) ifelse(x == "WS", "Wall strip", ifelse(x == "MS", "Microscope slide", x)),
            Date = function(x) paste("Week", x)
        )) +
        theme(strip.text = element_text(size = 10),
              axis.text.x = element_text(size = 10)) +
        labs(x = "Sample name", tag = tag) +
        theme(plot.tag = element_text(face = "bold"))

    p2 <- inner_join(rel_abund_ms, new_layer, by = "sample_id") %>%
        bar_plot(
            x_var = "sample_id", 
            position = position,
            global_colours = colorscheme, 
            italics = italics) + 
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
        theme(strip.text = element_text(size = 10),
              axis.text.x = element_text(size = 10)) +
        labs(x = "Sample name")

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

    df_long %>%
        group_by(taxon) %>%
        summarise(sum_lfc = sum(lfc), .groups = "drop") %>%
        arrange(desc(sum_lfc)) %>% 
        {order_taxa <<- .$taxon}

    df_long %>%
        mutate(taxon = factor(taxon, levels = order_taxa))
}

plastic_theme <- theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.margin = margin(20, 60, 20, 20),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title.align = 0.5,
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    plot.tag = element_text(face = "bold", size = 24)
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
