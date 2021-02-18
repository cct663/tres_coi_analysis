# This is the main analysis script for processing the COI data after running through
  # AMPtk. Just a few files are saved from AMPtk to use here. Figures and tables produced
  # in this script are saved in the repository and pasted in the Markdown document with
  # explanations.

################################################################################
# Load libraries ----

  pacman::p_load("tidyverse", "phyloseq", "plyr", "vegan", "here", "ggpubr", "igraph", "data.table", "dplyr")
    # tidyverse & plyr for data wrangling
    # phyloseq & vegan for community analyses and plotting
    # here for file reference by relative paths
    # ggpubr for plotting
    # igraph for network plots

    # Will need to install these libraries first if not yet installed
    # To install phyloseq, must install BiocManager and then use the following command: BiocManager::install(c("phyloseq"))

################################################################################
# Load and wrangle data ----

    # the prefix used for AMPtk processing
      # this is set in amptk and all files produced there have this prefix
      # when adapting this to a different set of data, change prefix
        amptk_prefix <- "trescoi"   
    
    # Load the number of reads by taxa per sample table. Format for phyloseq.
        otu_ab <- read.delim(here("1_raw_data", paste0(amptk_prefix, ".cluster.otu_table.txt")))
        rownames(otu_ab) <- otu_ab$X.OTU.ID   # give rownames from sample names
        otu_ab <- otu_ab[, 2:ncol(otu_ab)]    # remove the column of sample names
    
    # Read the mapping table
      # this 'mapping' table from AMPtk is mostly blank but has all the sample
        # names so it can be joined to actual sample metadata. It's also possible
        # to merge in the sample metadata within the AMPtk pipeline.
          map <- read.delim(here("1_raw_data", paste0(amptk_prefix, ".mapping_file.txt")))
          
            for(i in 1:nrow(map)){
              map$sampleID[i] <- strsplit(map$X.SampleID[i], "x")[[1]][2]
            }
    
          # write.table(map, "map.txt", sep = "\t") # if you want to save a copy of the mapping file
    
    # Read the otu taxonomy table
      otu_tax <- read.delim(here("1_raw_data", paste0(amptk_prefix, ".cluster.taxonomy.txt")))
      rownames(otu_tax) <- otu_tax$X.OTUID
    
    # The taxonomy result from AMPtk is in one long string of text. This is splitting up the string
      # and filling in a bunch of different columns. 
        for(i in 1:nrow(otu_tax)){
          temp <- otu_tax$taxonomy[i]
          temp2 <- strsplit(temp, "\\|")[[1]][3]
          temp3 <- strsplit(temp2, ":")[[1]][2]
          otu_tax$search_hit[i] <- strsplit(temp, "\\|")[[1]][1]
          otu_tax$hit_score[i] <- strsplit(temp, "\\|")[[1]][2]
          otu_tax$database[i] <- strsplit(temp2, ":")[[1]][1]
          otu_tax$accession[i] <- strsplit(temp3, ";")[[1]][1]
          otu_tax$kingdom[i] <- strsplit(strsplit(temp2, "k:")[[1]][2], ",")[[1]][1]
          otu_tax$phylum[i] <- strsplit(strsplit(temp2, "p:")[[1]][2], ",")[[1]][1]
          otu_tax$class[i] <- strsplit(strsplit(temp2, "c:")[[1]][2], ",")[[1]][1]
          otu_tax$order[i] <- strsplit(strsplit(temp2, "o:")[[1]][2], ",")[[1]][1]
          otu_tax$family[i] <- strsplit(strsplit(temp2, "f:")[[1]][2], ",")[[1]][1]
          otu_tax$genus[i] <- strsplit(strsplit(temp2, "g:")[[1]][2], ",")[[1]][1]
          otu_tax$species[i] <- strsplit(strsplit(temp2, "s:")[[1]][2], ",")[[1]][1]		
        }
    
        # Replace database mismatches caused by matches that aren't from BOLD records
            otu_tax$database <- gsub("None;k", "None", otu_tax$database)
            otu_tax$accession <- gsub("Animalia,p", "None", otu_tax$accession)
        
          # For phyloseq this has to be added as a matrix rather than data frame    
            otu_tax <- as.matrix(otu_tax)
          # This is saving just the taxonomic ranks rather than the database info.
            otu_tax2 <- otu_tax[, 10:16]
            
    # Add in the sample information to each sample
        s_info <- read.delim(here("1_raw_data", "tres_sample_info.txt"))    # prepared outside of AMPtk
        map2 <- join(map, s_info, "sampleID", "left", "first")
        rownames(map2) <- map2$X.SampleID
        
    # Print out some summary information about this dataset
        summary_table <- map2 %>% dplyr::count(site, age.1, cap_num)
        
        # This will identify how many samples there are in each category of site,
          # age, and capture number.

        # This will also identify samples that don't match the metadata file and write them as a separate output
          # missing <- subset(map2, is.na(map2$band) == TRUE)
          # write.table(missing, "missing_info.txt", sep = "\t")      

################################################################################
# Build initial phyloseq object ----
        
    OTU = otu_table(otu_ab, taxa_are_rows = TRUE)
    TAX = tax_table(otu_tax2)
    SAM = sample_data(map2)
    
    coi_ps <- phyloseq(OTU, TAX, SAM)

# Check sample effort ----
    
  # Subset just to arthropods
      coi_ps2 <- subset_taxa(coi_ps, phylum == "Arthropoda")

  # Histogram of sample reads. This is counting a sum of arthropod reads for each sample.
    depth <- data.frame(as(sample_data(coi_ps2), "data.frame"),
                        TotalReads = sample_sums(coi_ps2), keep.rownames = TRUE)
    p <- ggplot(depth, aes(log(TotalReads))) + geom_histogram(fill = "slateblue") + 
          ylab("Count of Samples") + xlab("log(Reads)") +
          theme_classic() + geom_vline(xintercept = log(150), linetype = "dashed", col = "coral3", size = 1) + 
          geom_text(x = log(150) - 0.2, y = 40, label = "150 Reads", angle = 90)
    p2 <- p + facet_grid(~ age)     # same but splitting out adult/nestling/negative_control
    
    # Save the histograms to file to be added to the markdown
      ggsave(here("3_r_scripts/total_reads.png"), plot = p, width = 8, height = 4.5, device = "png")
      ggsave(here("3_r_scripts/total_reads_split.png"), plot = p2, width = 8.2, height = 4, device = "png")

      # Extract the unique arthropod families found in the dataset
        unique_families <- get_taxa_unique(coi_ps2, taxonomic.rank = "family")
      # Save file with list of unique families to use to research aquatic vs.
        # terrestrial families.
        write.csv(unique_families, here("5_other_outputs/unique_families.csv"))
        # This file was then taken out of R to research aquatic vs. terrestrial
        # families, and will be re-imported later.
      
  # Rarefy to even depth of 150
    # running this would rarefy to an even depth moving forward
      # coi_ps2 <- rarefy_even_depth(coi_ps2, sample.size = 150, rngseed = 92)

################################################################################
# Agglomerate taxa ----
        
  # Depending on the analyses, we may want to agglomerate to different taxonomic ranks.
    # Many of the sequence IDs do not go all the way to species, so in those cases analyses
    # at the species level wouldn't include those sequences
        coi_genus <- tax_glom(coi_ps2, taxrank = "genus")
        coi_fam <- tax_glom(coi_ps2, taxrank = "family")
        coi_ord <- tax_glom(coi_ps2, taxrank = "order")

  # Merge family to life history
        life_history <- read.csv(here("5_other_outputs", "unique_families_aquatic_terrestrial_IDs.csv"))    # prepared outside of R
        otu_lh <- plyr::join(as.data.frame(tax_table(coi_fam)), life_history, "family", "left", "first")
        coi_fam2 <- phyloseq(
          otu_table(coi_fam),
          tax_table(as.matrix(otu_lh)),
          sample_data(coi_fam)
        )

        glom_ps <- coi_fam2    # change here which agglomeration you want to use for plots below
        
################################################################################
# Filtering criteria ----
        
  # Remove negative controls
      glom_ps <- subset_samples(glom_ps, age != "neg_control")
        
  # Remove OTUs with less than 5 reads in a sample aka 5-tons
        # (could change to singletons, 50-tons, or whatever)
      coi_ps2 <- prune_taxa(taxa_sums(glom_ps) > 5, glom_ps)

  # Transform to relative abundance
      coi_ra <- transform_sample_counts(coi_ps2, function(x) x / sum(x))
      
  ## Filter out taxa with relative abundance values below some threshold
      coi_ra2 <- filter_taxa(coi_ra, function(x) mean(x) > 1e-5, TRUE)

  # Take out the old test samples
      coi_ra2 <- subset_samples(coi_ra2, age.1 != "old")

  # Take out the testing samples from 2019 (were labelled as M19NXXX and do not have nest information)
      coi_ra2 <- subset_samples(coi_ra2, is.na(nest) == FALSE)

  # Transform to presence absence
      coi_pa <- transform_sample_counts(coi_ra2, function(x) ceiling(x))
      
  # Limit to genera in 20% of samples
      coi_20 <- prune_taxa(genefilter_sample(coi_pa, filterfun_sample(function(x) x > 0.1), A = 0.2 * nsamples(coi_pa)), coi_pa)
      
  # Limit to genera in 5% of samples (for network below)
      coi_05 <- prune_taxa(genefilter_sample(coi_pa, filterfun_sample(function(x) x > 0.1), A = 0.05 * nsamples(coi_pa)), coi_pa)
      
################################################################################
# Plot Patterns: Select objects to plot ----
      
      # Note that we will use coi_ra2 for plots for relative abundance, and 
      # coi_pa for plots for presence/absence.
      
################################################################################
# Plot Patterns: Richness ----
      
      # Richness by sample type. See help for many more options of different alpha metrics
      # Cannot use coi_ra2 here because plot_richness only accepts integers
        p <- plot_richness(coi_pa, x = "age.1", measures = c("Observed", "InvSimpson", "Shannon")) + 
              geom_boxplot(alpha = 0.3, aes(fill = age.1))
        # set order for ages for neatness in plots 
        age_order = c("6", "12", "15", "SY", "ASY", "AHY")
        p$data$age.1 <- as.character(p$data$age.1)
        p$data$age.1 <- factor(p$data$age.1, levels=age_order)
        ggsave(here("3_r_scripts/age_alpha.png"), p, width = 10, height = 4, device = "png")
        
    # Richness by age (nestling days and adult years) and site
        p <- plot_richness(coi_pa, x = "age.1", measures = c("Observed")) +
            geom_boxplot(alpha = 0.3, aes(fill = age.1)) + facet_wrap(~ site) +
            xlab("Age") +
            theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12)) + 
            theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
            theme(axis.title = element_text(size = 16)) +
            theme(strip.text = element_text(size = 16)) +
            guides(fill=guide_legend(title="Age"))
        # set order for sites for neatness in plots 
        site_order = c("Turkey_Hill", "Unit_4", "Unit_1", "Unit_2")
        p$data$site <- as.character(p$data$site)
        p$data$site <- factor(p$data$site, levels=site_order)
        # set order for ages for neatness in plots 
        age_order = c("6", "12", "15", "SY", "ASY", "AHY")
        p$data$age.1 <- as.character(p$data$age.1)
        p$data$age.1 <- factor(p$data$age.1, levels=age_order)
        ggsave(here("3_r_scripts/age_alpha_site.png"), p, width = 9, height = 6.5, device = "png")
        
      # Richness by age (just nestling vs. adult) and site
        p <- plot_richness(coi_pa, x = "age", measures = c("Observed")) +
          geom_boxplot(alpha = 0.3, aes(fill = age)) + facet_wrap(~ site) +
          xlab("Age") +
          theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12)) + 
          theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
          theme(axis.title = element_text(size = 16)) +
          theme(strip.text = element_text(size = 16)) +
          guides(fill=guide_legend(title="Age"))
        # set order for sites for neatness in plots 
        site_order = c("Turkey_Hill", "Unit_4", "Unit_1", "Unit_2")
        p$data$site <- as.character(p$data$site)
        p$data$site <- factor(p$data$site, levels=site_order)
        ggsave(here("3_r_scripts/age_nva_alpha_site.png"), p, width = 7, height = 7, device = "png")
        
    # Richness by adult capture number and site
        coi_adults <- subset_samples(coi_pa, cap_num != "")
        p <- plot_richness(coi_adults, x = "cap_num", measures = c("Observed")) + 
          geom_boxplot(alpha = 0.3, aes(fill = cap_num)) + facet_wrap(~ site)
        ggsave(here("3_r_scripts/capnum_alpha_site.png"), p, width = 9, height = 6.5, device = "png")
        
    # Richness by day of year
        # First look at richness with a series of box plots
        p <- plot_richness(coi_pa, x = "cap_doy", measures = c("Observed")) + 
          geom_boxplot(alpha = 0.3, aes(fill = age.1)) + facet_wrap(~ site)
        # Now look at richness as a loess curve
        richness <- data.table(p$data)
        richness$cap_doy <- as.numeric(richness$cap_doy)
        p <- ggplot(richness, mapping = aes(x = cap_doy, y = value)) + geom_point() + geom_smooth(method = "loess") +
              xlab("Capture day of year") + ylab("Alpha Diversity Measure") + facet_wrap(~ age, ncol = 1) +
              theme_classic() +
              theme(axis.text = element_text(size = 14)) + 
              theme(axis.title = element_text(size = 16)) +
              theme(strip.text = element_text(size = 16))
        ggsave(here("3_r_scripts/age_site_richness.png"), p, width = 8.7, height = 6.8, device = "png")

################################################################################
# Plot Patterns: Taxonomic Groups ----        
    # All genera
        p <- plot_bar(coi_pa, "family") + theme_classic() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
              geom_hline(yintercept = 393 * 0.1, linetype = "dotted", col = "coral3") + 
              geom_hline(yintercept = 393 * 0.2, linetype = "dotted", col = "coral3") + 
              geom_hline(yintercept = 393 * 0.3, linetype = "dotted", col = "coral3")
        ggsave(here("3_r_scripts/family_bar.png"), width = 10, height = 4.5, device = "png")
        
    # Genera over 20% split by age
        p <- plot_bar(coi_20, "family", fill="life_history") + theme_classic() + theme(axis.text.x = element_text(angle = 90)) + 
            facet_wrap(~ age, ncol = 1)
        ggsave(here("3_r_scripts/common_families.png"), width = 10, height = 6, device = "png")
        # Create a figure that has families color coded by life history
        # phyloseq automatically draws black borders around every sample which makes it impossible to see the colors of 
        # aquatic vs. terrestrial samples.
        # We will have to redefine the function to get it to stop doing this (I found this code on stack overflow)
        plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
        {
          mdf = psmelt(physeq)
          p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
          p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
          p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
          if (!is.null(facet_grid)) {
            p <- p + facet_grid(facet_grid)
          }
          if (!is.null(title)) {
            p <- p + ggtitle(title)
          }
          return(p)
        }
        
        p <- plot_bar_2(coi_20, x="family", fill="life_history") + 
            xlab("Family") +
            ylab("Number of samples detected in") +
            theme_classic() + theme(axis.text.x = element_text(angle = 90, size = 12)) + 
            theme(axis.title = element_text(size = 14)) +
            theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12)) +
            theme(strip.text = element_text(size = 12)) +
            facet_wrap(~ age, ncol = 1) +
            scale_fill_manual("Life history", values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4", 
                                                   "unknown" = "gray"))
        ggsave(here("3_r_scripts/common_families.png"), width = 10, height = 6, device = "png")
        
################################################################################
# Plot Patterns: Aquatic vs. Terrestrial ----
    # Plot some aquatic vs. terrestrial plots across age and site
        
        # Plot aquatic vs. terrestrial for sites, adults vs. nestlings
        plot_ra <- psmelt(coi_ra2) # psmelt makes a phyloseq object into a data frame
        p <- ggplot(plot_ra, aes(x=age, y=Abundance, fill=life_history)) +
          geom_col(position="fill") +
          facet_wrap(~ site, ncol = 2) +
          xlab("Age") +
          ylab("Relative Abundance") +
          theme_classic() +
          theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
          theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
          theme(strip.text = element_text(size = 14)) +
          scale_fill_manual("Life history", values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4", 
                                                       "unknown" = "gray"))
          # set order for sites for neatness in plots 
          site_order = c("Turkey_Hill", "Unit_4", "Unit_1", "Unit_2")
          p$data$site <- as.character(p$data$site)
          p$data$site <- factor(p$data$site, levels=site_order)
          
        ggsave(here("3_r_scripts/lh_age_site.png"), p, width = 7, height = 7, device = "png")
          
        # Plot aquatic vs. terrestrial for adults and nestlings full ages
          p1 <- ggplot(plot_ra, aes(x=age.1, y=Abundance, fill=life_history)) +
            geom_col(position="fill") +
            facet_wrap(~ site, ncol = 2) +
            xlab("Age") +
            ylab("Relative Abundance") +
            theme_classic() +
            theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
            theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
            theme(strip.text = element_text(size = 14)) +
            scale_fill_manual("Life history", values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4", 
                                                       "unknown" = "gray"))
          # set order for ages for neatness in plots
          age_order = c("6", "12", "15", "SY", "ASY", "AHY")
          p1$data$age.1 <- as.character(p$data$age.1)
          p1$data$age.1 <- factor(p$data$age.1, levels=age_order) 
        
        ggsave(here("3_r_scripts/lh_age1_site.png"), p1, width = 10, height = 8, device = "png")
        
        # Plot aquatic vs. terrestrial across days
        p2 <- ggplot(plot_ra, aes(x=cap_doy, y=Abundance, fill=life_history)) +
          geom_col(position="fill") +
          xlab("Day of Year") +
          ylab("Relative Abundance") +
          theme_classic() +
          theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
          theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
          theme(strip.text = element_text(size = 14)) +
          scale_fill_manual("Life history", values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4", 
                                                       "unknown" = "gray"))
        
        ggsave(here("3_r_scripts/lh_doy.png"), p2, width = 10, height = 8, device = "png")
        
        # Plot aquatic vs. terrestrial across days and separated by adults vs. nestlings
        p3 <- ggplot(plot_ra, aes(x=cap_doy, y=Abundance, fill=life_history)) +
          geom_col(position="fill") +
          facet_wrap(~ age, ncol = 1) +
          xlab("Day of Year") +
          ylab("Relative Abundance") +
          theme_classic() +
          theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(angle = 90, size = 14)) +
          theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
          theme(strip.text = element_text(size = 14)) +
          scale_fill_manual("Life history", values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4", 
                                                       "unknown" = "gray"))
        
        ggsave(here("3_r_scripts/lh_doy_age.png"), p3, width = 10, height = 8, device = "png")

################################################################################
# Plot Patterns: Ordinations/Beta Diversity ----       
    # Plot some ordinations using relative abundance and presence/absence
        ord_data_ra <- ordinate(coi_ra2, method = "PCoA", distance = "bray")
        ord_data_pa <- ordinate(coi_pa, method = "PCoA", distance = "bray")
        
      # Site -- relative abundance
          p1_ra <- plot_ordination(coi_ra2, ord_data_ra, color = "site", title = "Bray-Curtis PCoA") + 
                geom_point(size = 2, alpha = .4) + theme_classic() +
                stat_ellipse(aes(group = site), level = 0.9) +
                theme(axis.title = element_text(size = 16)) +
                theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14))
          ggsave(here("3_r_scripts/site_ordinate_ra.png"), p1_ra, width = 7, height = 5.4, device = "png")
        
      # Site -- presence/absence
          p1_pa <- plot_ordination(coi_pa, ord_data_pa, color = "site", title = "Bray-Curtis PCoA") + 
            geom_point(size = 2, alpha = .4) + theme_classic() +
            stat_ellipse(aes(group = site), level = 0.9) +
            theme(axis.title = element_text(size = 16)) +
            theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14))
          ggsave(here("3_r_scripts/site_ordinate_pa.png"), p1_pa, width = 7, height = 5.4, device = "png")
          
      # Age -- relative abundance
          p2_ra <- plot_ordination(coi_ra2, ord_data_ra, color = "age", title = "Bray-Curtis PCoA") + 
            geom_point(size = 2, alpha = .4) + theme_classic() +
            stat_ellipse(aes(group = age), level = 0.9)
          ggsave(here("3_r_scripts/age_ordinate_ra.png"), p2_ra, width = 7, height = 5.4, device = "png")
          
      # Age -- presence absence
          p2_pa <- plot_ordination(coi_pa, ord_data_pa, color = "age", title = "Bray-Curtis PCoA") + 
            geom_point(size = 2, alpha = .4) + theme_classic() +
            stat_ellipse(aes(group = age), level = 0.9)
          ggsave(here("3_r_scripts/age_ordinate_pa.png"), p2_pa, width = 7, height = 5.4, device = "png")
          
          
          ###############
          
      # Age and site
          p3 <- p + facet_wrap(~ age) + 
            theme(axis.title = element_text(size = 16)) +
            theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
            theme(strip.text = element_text(size = 14))
          ggsave(here("3_r_scripts/age_site_ordinate.png"), p3, width = 8.2, height = 4.8, device = "png")
          
      # Site and age
          p4 <- p2 + facet_wrap(~ site) +
            theme(axis.title = element_text(size = 16)) +
            theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
            theme(strip.text = element_text(size = 14))
          ggsave(here("3_r_scripts/site_age_ordinate.png"), p4, width = 8.2, height = 6.2, device = "png")


################################################################################
################################################################################
# Point at which Jenny stopped using code
          
################################################################################
# Plot Patterns: Adult vs. Nestling Diets ----
    # Plot adult vs. nestling diets
          adult_20 <- subset_samples(coi_20, age == "Adult")
          six_20 <- subset_samples(coi_20, age.1 == "6")
          twelve_20 <- subset_samples(coi_20, age.1 == "12")
          fifteen_20 <- subset_samples(coi_20, age.1 == "15")
          compare <- data.frame(adult = taxa_sums(adult_20) / nrow(sample_data(adult_20)),
                                six = taxa_sums(six_20) / nrow(sample_data(six_20)),
                                twel = taxa_sums(twelve_20) / nrow(sample_data(twelve_20)),
                                fift = taxa_sums(fifteen_20) / nrow(sample_data(fifteen_20)))
          p1 <- ggplot(compare, aes(x = adult, y = fift)) + geom_point(col = "slateblue") + 
                geom_smooth(col = "coral3", method = "lm") + theme_classic() + xlab("Adults") + ylab("15-Day Nestlings") +
                xlim(0.05, 0.9) + ylim(0.05, 0.9)
          p2 <- ggplot(compare, aes(x = adult, y = twel)) + geom_point(col = "slateblue") + 
                geom_smooth(col = "coral3", method = "lm") + theme_classic() + xlab("Adults") + ylab("12-Day Nestlings") +
                xlim(0.05, 0.9) + ylim(0.05, 0.9)
          p3 <- ggplot(compare, aes(x = twel, y = fift)) + geom_point(col = "slateblue") + 
                geom_smooth(col = "coral3", method = "lm") + theme_classic() + xlab("12-Day Nestlings") + ylab("15-Day Nestlings") +
                xlim(0.05, 0.9) + ylim(0.05, 0.9)
          p4 <- ggplot(compare, aes(x = adult, y = six)) + geom_point(col = "slateblue") + 
              geom_smooth(col = "coral3", method = "lm") + theme_classic() + xlab("Adults") + ylab("6-Day Nestlings") +
              xlim(0.05, 0.9) + ylim(0.05, 0.9)
          p5 <- ggplot(compare, aes(x = twel, y = six)) + geom_point(col = "slateblue") + 
              geom_smooth(col = "coral3", method = "lm") + theme_classic() + xlab("12-Day Nestlings") + ylab("6-Day Nestlings") +
              xlim(0.05, 0.9) + ylim(0.05, 0.9)
          p6 <- ggplot(compare, aes(x = fift, y = six)) + geom_point(col = "slateblue") + 
             geom_smooth(col = "coral3", method = "lm") + theme_classic() + xlab("15-Day Nestlings") + ylab("6-Day Nestlings") +
              xlim(0.05, 0.9) + ylim(0.05, 0.9)
          p7 <- ggarrange(p1, p2, p3, p4, p5, p6, nrow = 1)
          ggsave(here("3_r_scripts/compare.png"), p7, width = 11, height = 3.4, device = "png")

    # Plot network
          p <- plot_net(glom_ps, maxdist = 0.4, point_label = "nest", color = "site")
          ggsave(here("3_r_scripts/network.png"), p, width = 7.5, height = 6.5, device = "png")

################################################################################
# Plot Patterns: Order over Season ----          
    # Plot season at order level (gets it down to 17 taxa like Ryan & Lily's Ecology Letters paper)
          # Remove negative controls
            ord_ps <- subset_samples(coi_ord, age != "neg_control" & age.1 != "old")
          # Transform to relative abundance
            ord_ra <- transform_sample_counts(ord_ps, function(x) x / sum(x))
          # Filter out taxa with relative abundance values below some threshold
            ord_ra2 <- filter_taxa(ord_ra, function(x) mean(x) > 1e-5, TRUE)
          # Filter out taxa that aren't in 5% of samples
            ord_ra2 <- prune_taxa(genefilter_sample(ord_ra2, filterfun_sample(function(x) x > 0.01), A = 0.05 * nsamples(ord_ra2)), ord_ra2)
          # Transform to presence absence
            ord_pa <- transform_sample_counts(ord_ra2, function(x) ceiling(x))
          
          season <- as.data.frame(otu_table(ord_pa))
          season$otu <- rownames(season, )
          season2 <- pivot_longer(data = season, !otu, names_to = "X.SampleID", values_to = "present")
          season2 <- join(season2, map2, "X.SampleID")
          map2$cap_doy <- as.numeric(map2$cap_doy)
          tax <- as.data.frame(tax_table(ord_pa))
          tax$otu <- rownames(tax)
          tax <- tax[, c("otu", "order")]
          season2 <- join(season2, tax, "otu")
          
        # I do not understand why, but these plots will only work right now if you run this chunk of the code twice.
          # You have to run it from the line that removes negative control twice.
          p1 <- ggplot(season2, mapping = aes(x = cap_doy, y = present, col = order)) + geom_smooth(method = "loess", se = FALSE) +
            theme_classic() + xlab("Capture day of year") + ylab("Percent of samples detected") +
            facet_wrap(~ age, ncol = 1) + ylim(0, 1) +
            theme(axis.title = element_text(size = 16)) + theme(axis.text.x = element_text(size = 14)) +
            theme(legend.title = element_text(size = 16), legend.text = element_text(size = 14)) +
            theme(strip.text = element_text(size = 16))
          ggsave(here("3_r_scripts/age_season.png"), p1, width = 8.7, height = 6.8, device = "png")
          
          p2 <- ggplot(season2, mapping = aes(x = cap_doy, y = present, col = order)) + geom_smooth(method = "loess", se = FALSE) +
            theme_classic() + xlab("Capture day of year") + ylab("Percent of samples detected") +
            facet_wrap(~ site) + ylim(0, 1)
          ggsave(here("3_r_scripts/site_season.png"), p2, width = 11, height = 6.5, device = "png")

################################################################################
# Plot Patterns: Network of Food Items ----
          
# Make a network of food items ----
    # This will require some wrangling so putting it in a new section
        # The phyloseq taxa table is an incidence matrix, convert to network
          
        # This is kind of ugly and not well annotated. Can be ignored.
          # There is probably a more graceful way to do this but I'm using brute force
          # to build an adjacency matrix from the otu table in the phyloseq object and then
          # plotting a network based on which taxa are found together in the same fecal samples.
          
          use_this <- coi_05
          net <- matrix(nrow = nrow(otu_table(use_this)), ncol = nrow(otu_table(use_this)))
          colnames(net) <- rownames(otu_table(use_this))
          rownames(net) <- rownames(otu_table(use_this))
          temp <- otu_table(use_this)
          
          for(i in 1:nrow(net)){
            for(k in 1:ncol(net)){
              include <- c(rownames(net)[i], colnames(net)[k])
              temp2 <- subset(temp, rownames(temp) %in% include)
              if(nrow(temp2) > 1){
                temp3 <- as.vector(temp2[1, ] + temp2[2, ])
                net[i, k] <- length(subset(temp3, temp3 == 2))}
              if(rownames(net)[i] == colnames(net)[k]){net[i, k] <- ncol(otu_table(coi_pa))}
            }
            print(paste(i, " of ", nrow(net)))
          }
          
          net2 <- net / 394
          for(i in 1:nrow(net2)){
            for(k in 1:nrow(net2)){
              if(net2[i, k] < 0.05){net2[i, k] <- 0}
            }
          }
          
          network <- graph_from_adjacency_matrix(net2, weighted = TRUE, diag = FALSE, mode = "undirected")
          cl = cluster_fast_greedy(network)
          coords = layout_nicely(network)
          plot.igraph(network, layout = coords,
                      edge.color = "black", vertex.size = 4, vertex.label = NA)
          
          
          png(here("3_r_scripts/food_network.png"), width = 8.8, height = 8, units = "in", res = 300)
            plot(cl, network, layout = coords, vertex.label = NA)
          dev.off()