# Load all packages/libraries needed to run script

library(readr) # read CSV files into R
library(plyr) # for data manipulation based on the split-apply-combine strategy
library(dplyr) # to filter, arrange, select, and summarize data
library(tidyr) # for reshaping and cleaning data frames
library(tidyverse) # data manipulation, organization, & visualization
library(ggplot2) # data visualizations 
library(car) # provide functions and tools for regression analysis
library(vegan) # contains methods of multivariate analysis
library(emmeans) # provides methods for performing various post-hoc comparisons, such as pairwise comparisons
install.packages("remotes")
remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") # post-hoc analysis when an overall PERMANOVA test indicates significant differences among groups
library(pairwiseAdonis)
library(lme4)
library(piecewiseSEM) # structural equation modeling
library(lmerTest) # provides p-values and F/t tests for linear mixed models
library(cowplot) # for combining figures
library(gridExtra) # arrange multiple grid-based plots on a page
library(eulerr) # for Venn diagram
library(grid)


# Setting Relative path name
# Define file path 
# This ensures the script can locate the data file regardless of working directory
# by referencing the relative path within the project folder structure

soil<-file.path(".","data","Ilex_Soil_Arthropods.csv")
weights<-file.path(".","data","Berlese_Soil_Weights.csv")
traits<-file.path(".","data","Ilex_traits.csv")
Fungi<-file.path(".","data","Fungal_OTUs.csv")
fungi.cla<-file.path(".","data","Fungal_taxon.csv")
microbial<-file.path(".","data","Microbial_Biomass_Data.csv")
biogeo<-file.path(".","data","soil_biogeochem_properties.csv")
organic<-file.path(".","data","Soil_organic_matter_fractions.csv")
bact<-file.path(".","data","Bacterial_OTUs.csv")
bact_tax<-file.path(".","data","Bacteria_taxon.csv")



# Load the datasets into R
# Reads the CSV file specified by the 'XXX' file path
# and stores it as a data frame called 'XXX'

arthropods<-read_csv(soil)
soil_weight<-read_csv(weights)
Ilex_traits<-read_csv(traits)
fungi<-read_csv(Fungi)
fungal_c<-read_csv(fungi.cla)
microbial_mass<-read_csv(microbial)
biogeochemical<-read_csv(biogeo)
org_matter<-read_csv(organic)
bacteria<-read_csv(bact)
tax_bac<-read_csv(bact_tax)


## selecting best traits among our Ilex vomitoria traits to account for plant size for all analyses

traits <- Ilex_traits[, c("J_Height", "Avg_Crown", "Avg_circ", "N_stems")]

cor(traits, use = "complete.obs")

pairs(traits)

# assessing whether female and male plants differ in traits
height_model <- lm(J_Height ~ Sex * Location, data = bacteria_model_data)
car::Anova(height_model, type = 3)
summary(height_model)
plot(height_model)

crown_model <- lm(Avg_Crown ~ Sex * Location, data = bacteria_model_data)
car::Anova(crown_model, type = 3)
summary(crown_model)
plot(crown_model)

circ_model <- lm(Avg_circ ~ Sex * Location, data = bacteria_model_data)
car::Anova(circ_model, type = 3)
summary(circ_model)
plot(circ_model)

stems_model <- lm(N_stems ~ Sex * Location, data = bacteria_model_data)
car::Anova(stems_model, type = 3)
summary(stems_model)
plot(stems_model)





# Question #1: How do male and female plants differ in soil microbial and arthropod diversity and community composition?  


# BACTERIA DATA

# Assign letters to sites to sort them from driest to wettest

bacteria$Location <- dplyr::recode(bacteria$Location,
                                   "Wi" = "AWi",
                                   "H" = "BH",
                                   "TC" = "CTC")

# Select all columns except Location & Sex and create new data set
bacteria_n <-subset(bacteria, select = -c(Location,Sex))


# Calculate total sequencing depth for each sample
# This sums all bacterial OTU/ASV read counts per sample
# across all taxa columns, excluding the sample identifier (Log_ID)

seq_summary <- bacteria_n %>%
  
# Apply calculations row-by-row (per sample)  
rowwise() %>%
  
# Create a new column ('total_sequences') containing
# the sum of all sequence reads in each sample
mutate(total_sequences = sum(c_across(-Log_ID))) %>%
  
# Remove rowwise grouping to restore standard data frame behavior
ungroup()

# Summarize sequencing depth across all samples
# This calculates the minimum, maximum, and median number
# of sequence reads per sample based on total sequencing depth
sequence_stats <- seq_summary %>%
  summarize(
    min_sequences = min(total_sequences),
    max_sequences = max(total_sequences),
    median_sequences = median(total_sequences)
  )

# Print results
sequence_stats


# Bacteria richness 

# Calculate bacterial taxonomic richness for each sample
# Richness is defined as the number of taxa with non-zero read counts
# Excludes metadata columns (Log_ID, Location, Sex) from the calculation
bacteria_data_r <- bacteria %>%
  
# Create a new column ('richness') containing
# the count of taxa present (abundance > 0) in each sample
mutate(richness = rowSums(across(-c(Log_ID, Location, Sex)) > 0))


# Combine bacteria richness data & ilex v. traits to run models with some covariates

new_ilex_t<-Ilex_traits%>%
  select(c(Log_ID,J_Height,Avg_Crown,Avg_circ,N_stems))

new_bacteria_r <- merge(x=bacteria_data_r,y=new_ilex_t,
                        by=c("Log_ID"))

bacteria_model_data <- new_bacteria_r %>%
  dplyr::select(
    richness,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()


# ANOVA: Bacteria richness


# Fit linear model testing Sex, Location, and interaction
modelbrs<-lm(richness ~ Sex * Location + N_stems , data = bacteria_model_data)
# Type III ANOVA (for unbalanced designs)
Anova(modelbrs, type = 3)
# Summary of model coefficients
summary(modelbrs)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelbrs)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelbrs, pairwise ~ Sex * Location, adjust = "tukey")




# Bacteria alpha diversity- Shannon

# Calculate Shannon diversity index for bacterial communities in each sample
# Shannon diversity accounts for both taxonomic richness and evenness
# Excludes metadata and previously calculated richness from the calculation

bacteria_data_s <- bacteria_data_r %>%
  
# Create a new column ('Shannon') containing
# the Shannon diversity index for each sample
mutate(Shannon = apply(across(-c(Log_ID, Location, Sex,richness)), 1, 
                         function(x) diversity(x, index = "shannon")))


# Combine bacteria richness data & ilex v. traits to run models with some covariates

new_ilex_t<-Ilex_traits%>%
  select(c(Log_ID,J_Height,Avg_Crown,Avg_circ,N_stems))

new_bacteria_s <- merge(x=bacteria_data_s,y=new_ilex_t,
                        by=c("Log_ID"))

bacteria_shannon_data <- new_bacteria_s %>%
  dplyr::select(
    Shannon,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()

# Bacterial Shannon diversity

# Fit linear model testing Sex, Location, and interaction
modelbss<-lm(Shannon ~ Sex * Location + N_stems, data = bacteria_shannon_data)
# Type III ANOVA (for unbalanced designs)
Anova(modelbss, type = 3)
# Summary of model coefficients
summary(modelbss)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelbss)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelbss, pairwise ~ Sex * Location, adjust = "tukey")





# FUNGI DATA

# Fungal richness: rarefaction and data preparation

# Remove metadata columns (retain only OTU relative abundance data)
fungi.new <- fungi[, -c(1:3)]
# Identify the minimum sequencing depth across samples
# (used as the rarefaction target)
min_val <- min(rowSums(fungi.new))
min_val
# Rarefy all samples to the same sequencing depth
# (controls for uneven sequencing effort)
rarefied_data <- rrarefy(fungi.new, sample = min_val)
# Extract metadata columns
fungi.metadata <- fungi[, c(1:3)]
# Recombine metadata with rarefied abundance data
fungi_combined <- cbind(fungi.metadata, rarefied_data)

# Recode location names 
fungi_combined$Location <- dplyr::recode(fungi_combined$Location,
                                         "Wi" = "AWi",
                                         "H" = "BH",
                                         "TC" = "CTC")
# Calculate fungal richness
# (number of taxa with nonzero abundance per sample)
fungi_data_r <- fungi_combined %>%
  mutate(richness = rowSums(across(-c(Log_ID, Location, Sex)) > 0))

# Combine fungi richness data & Ilex v. traits to run models with some covariates


new_fungi_r <- merge(x=fungi_data_r,y=new_ilex_t,
                        by=c("Log_ID"))

fungi_rich_data <- new_fungi_r %>%
  dplyr::select(
    richness,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()


# ANOVA: Fungal richness + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelFrs<-lm(richness ~ Sex * Location + N_stems , data = fungi_rich_data)
# Type III ANOVA (for unbalanced designs)
Anova(modelFrs, type = 3)
# Summary of model coefficients
summary(modelFrs)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelFrs)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelFrs, pairwise ~ Sex * Location, adjust = "tukey")




# Fungal alpha diversity: Shannon index

# Calculate Shannon diversity for each sample
fungi_data_s <- fungi_data_r %>%
  mutate(Shannon = apply(across(-c(Log_ID, Location, Sex,richness)), 1, 
                         function(x) diversity(x, index = "shannon")))


new_fungi_s <- merge(x=fungi_data_s,y=new_ilex_t,
                        by=c("Log_ID"))

fungi_shannon_data <- new_fungi_s %>%
  dplyr::select(
    Shannon,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()


# ANOVA: Fungi Shannon diversity + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelFss<-lm(Shannon ~ Sex * Location + N_stems, data = fungi_shannon_data)
# Type III ANOVA (for unbalanced designs)
Anova(modelFss, type = 3)
# Summary of model coefficients
summary(modelFss)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelFss)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelFss, pairwise ~ Sex * Location, adjust = "tukey")






# ARTHROPOD DATA

# Remove H3-F2d, this sample was left on the berlese without heat for ~8 hours.
arthropods<-arthropods[!(arthropods$Log_ID=="H3-F2"),]

# Remove samples with no arthropods
arthropods_n<-na.omit(arthropods)

# merge soil weight and arthropod data
arthropods2 <- merge(x=arthropods_n,y=soil_weight,
                     by=c("Log_ID"))

# Calculate the count of arthropod per soil weight
arthropods2$Countperweight<-arthropods2$Count/arthropods2$Mass

# Total count of arthropods

# select the needed columns
arthropod_final<-subset(arthropods2, select = -c(Collection_Date,
                                                 Analyzer,Analysis_Date))


# Sum the total count of arthropods by log_ID (sample)
Total_count <-aggregate(cbind(Countperweight) ~ Log_ID+Mass, 
                        data = arthropod_final, FUN = sum)

# recover the male and female status from Log_ID
Total_count$Sex <- ifelse(grepl("M", Total_count$Log_ID), "M","F")

#recover the location from Log_ID
# Put an A, B, and C in front of the labels to force them to sort 
# in the order that I want. I organized them to be left to right the same way
# the precipitation gradient is west to east (Wildflower,Huntsville, Turkey Creek)
Total_count$Location <- with(Total_count, 
                             ifelse(grepl("TC", Total_count$Log_ID), "CTC",
                                    ifelse(grepl("H", Total_count$Log_ID), "BH", 'AWi')))



# Arthropod richness

# Remove not needed columns
Richness <-subset(arthropod_final, select = -c(Countperweight,Mass))

# Convert data to wide format
Richness_n<-pivot_wider(Richness,names_from = Taxa,values_from = Count, 
                      values_fill = 0)

# Remove not needed cols
Richness_new <-subset(Richness_n, select = -c(Location, Sex))
# If the cell contains a number it gets counted if it doesn't it's not counted
Richness2 <-ddply(Richness_new,~Log_ID,function(x) {
  + data.frame(richness=sum(x[-1]>0))})
# Recover Location and Sex
Richness2$Location <- with(Total_count, ifelse(grepl("TC", Richness2$Log_ID), "CTC",
                                               ifelse(grepl("H", Richness2$Log_ID), "BH", 'AWi')))
Richness2$Sex <- ifelse(grepl("M", Richness2$Log_ID), "M","F")


# Combine fungi richness data & Ilex v. traits to run models with some covariates


new_arthro_r <- merge(x=Richness2,y=new_ilex_t,
                     by=c("Log_ID"))

arthro_rich_data <- new_arthro_r %>%
  dplyr::select(
    richness,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()


# ANOVA: Arthropod richness + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelAs<-lm(richness ~ Sex * Location + N_stems, data = arthro_rich_data)
# Type III ANOVA (for unbalanced designs)
Anova(modelAs, type = 3)
# Summary of model coefficients
summary(modelAs)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelAs)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelAs, pairwise ~ Sex * Location, adjust = "tukey")




# Arthropod alpha diversity

# Remove extra columns

Art_Diversity <-subset(arthropods2, select = -c(Collection_Date,Analyzer,
                                                Analysis_Date,Count, Mass))

# Spread into wide format and fill with 0
# remove duplicated rows caused by that one sample that had a hollowed out beetle
# in it
Art_Diversity_n<-pivot_wider(Art_Diversity,names_from = Taxa,values_from = Countperweight, 
                           values_fill = 0)
# Select needed column to calculate diversity
Art_Diversity_new <-subset(Art_Diversity_n, select = -c(Location,Sex))

# Use the diversity function from Vegan to calculate the Shannon wiener index
Shannon <-diversity(Art_Diversity_new[-1], index = "shannon")
# This output should just be a list of numbers that falls in the same order as
# the original Log_IDs did
# Make it into a data frame
Shannon <- as.data.frame(Shannon)
# Recover the location and sex data from the Richness section. If you
# did not run that section this wont work!
Shannon$Location <-Richness2$Location
Shannon$Sex <-Richness2$Sex
Shannon$Log_ID <-Richness2$Log_ID

new_arthro_s <- merge(x=Shannon,y=new_ilex_t,
                     by=c("Log_ID"))

arthro_shannon_data <- new_arthro_s %>%
  dplyr::select(
    Shannon,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()

# ANOVA: Arthropod alpha diversity + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelAst<-lm(Shannon ~ Sex * Location + N_stems, data = arthro_shannon_data)
# Type III ANOVA (for unbalanced designs)
Anova(modelAst, type = 3)
# Summary of model coefficients
summary(modelAst)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelAst)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelAst, pairwise ~ Sex * Location, adjust = "tukey")


# Figure #1 Boxplot multipanel richness & diversity


base_theme <- theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "lines"),
    legend.position = "none",
    plot.margin = unit(c(5, 5, 5, 35), "pt")  # increased left margin
  )

# Significance annotations
b_r_star <- data.frame(
  Location = "CTC",
  y = max(bacteria_data_r$richness, na.rm = TRUE) * 1.08,
  label = "**"
)

b_s_star <- data.frame(
  Location = "CTC",
  y = max(bacteria_data_s$Shannon, na.rm = TRUE) * 1.08,
  label = "**"
)

f_s_star <- data.frame(
  Location = "CTC",
  y = max(fungi_data_s$Shannon, na.rm = TRUE) * 1.08,
  label = "*"
)

# Bacteria richness
plot_b_r <- ggplot(bacteria_data_r,
                   aes(x = Location, y = richness, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85),
               outlier.shape = NA) +
  geom_jitter(
    aes(fill = Sex),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.85
    ),
    alpha = 0.6,
    size = 2,
    shape = 21,
    stroke = 0.4,
    color = "black"
  ) +
  geom_text(
    data = b_r_star,
    aes(x = Location, y = y, label = label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Site", y = "Richness", fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male")
  ) +
  scale_x_discrete(
    labels = c("AWi" = "LBJWC",
               "BH" = "PERL",
               "CTC" = "BTNP")
  ) +
  base_theme


# Bacteria Shannon
plot_b_s <- ggplot(bacteria_data_s,
                   aes(x = Location, y = Shannon, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85),
               outlier.shape = NA) +
  geom_jitter(
    aes(fill = Sex),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.85
    ),
    alpha = 0.6,
    size = 2,
    shape = 21,
    stroke = 0.4,
    color = "black"
  ) +
  geom_text(
    data = b_s_star,
    aes(x = Location, y = y, label = label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Site", y = "Shannon Diversity", fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male")
  ) +
  scale_x_discrete(
    labels = c("AWi" = "LBJWC",
               "BH" = "PERL",
               "CTC" = "BTNP")
  ) +
  base_theme


# Fungi richness
plot_f_r <- ggplot(fungi_data_r, aes(x = Location, y = richness, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    aes(fill = Sex),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  labs(x = "Site", y = "Richness", fill = "Sex") +
  scale_fill_manual(values = c("F" = "tomato", "M" = "darkcyan"),
                    labels = c("F" = "Female", "M" = "Male")) +
  scale_x_discrete(labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")) +
  base_theme
plot_f_r

# Fungi Shannon
plot_f_s <- ggplot(fungi_data_s,
                   aes(x = Location, y = Shannon, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85),
               outlier.shape = NA) +
  geom_jitter(
    aes(fill = Sex),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.85
    ),
    alpha = 0.6,
    size = 2,
    shape = 21,
    stroke = 0.4,
    color = "black"
  ) +
  geom_text(
    data = f_s_star,
    aes(x = Location, y = y, label = label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Site", y = "Shannon Diversity", fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male")
  ) +
  scale_x_discrete(
    labels = c("AWi" = "LBJWC",
               "BH" = "PERL",
               "CTC" = "BTNP")
  ) +
  base_theme

# Arthropod richness
plot_a_r <- ggplot(Richness2, aes(x = Location, y = richness, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    aes(fill = Sex),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  labs(x = "Site", y = "Richness", fill = "Sex") +
  scale_fill_manual(values = c("F" = "tomato", "M" = "darkcyan"),
                    labels = c("F" = "Female", "M" = "Male")) +
  scale_x_discrete(labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")) +
  base_theme

# Arthropod Shannon
plot_a_s <- ggplot(Shannon, aes(x = Location, y = Shannon, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    aes(fill = Sex),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  labs(x = "Site", y = "Shannon Diversity", fill = "Sex") +
  scale_fill_manual(values = c("F" = "tomato", "M" = "darkcyan"),
                    labels = c("F" = "Female", "M" = "Male")) +
  scale_x_discrete(labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")) +
  base_theme

# Function to format panels with individual y axis labels
format_panel <- function(p, ylab, show_x = FALSE) {
  p +
    labs(y = ylab) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text.x = if (show_x) element_text(size = 12) else element_blank(),
      axis.text.y = element_text(size = 12)
    )
}

# Apply y labels per plot
plot_b_r2 <- format_panel(plot_b_r, "Bacteria Richness", show_x = FALSE)
plot_b_s2 <- format_panel(plot_b_s, "Bacteria Diversity", show_x = FALSE)
plot_f_r2 <- format_panel(plot_f_r, "Fungi Richness", show_x = FALSE)
plot_f_s2 <- format_panel(plot_f_s, "Fungi Diversity", show_x = FALSE)
plot_a_r2 <- format_panel(plot_a_r, "Arthropod Richness", show_x = TRUE)
plot_a_s2 <- format_panel(plot_a_s, "Arthropod Diversity", show_x = TRUE)

# Extract one legend
legend_box <- get_legend(
  plot_b_r + theme(legend.position = "right",
                   legend.box.margin = margin(0, 0, 0, 0))
)

# Build the grid
alpha_div <- plot_grid(
  plot_b_r2, plot_b_s2,
  plot_f_r2, plot_f_s2,
  plot_a_r2, plot_a_s2,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E", "F"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.04,
  label_y = 1.02,
  align = "hv"
)

# attach legend
rel_widths <- c(1, 0.15)
alpha_with_legend <- plot_grid(alpha_div, legend_box, ncol = 2, rel_widths = rel_widths)

# Final plot (no shared labels since each panel has own y axis)
final_plot <- ggdraw(alpha_with_legend)

final_plot

# save plots
ggsave("alpha_multipanel.tiff",final_plot, width = 8, height = 8, device='tiff', dpi=600)

ggsave("alpha_multifinal.png",final_plot, width = 9, height = 8, device='png', dpi=600)




# Beta Diversity

# Bacteria Beta diversity

# Create a dataset containing only bacterial abundance data
# Removes metadata variables and diversity metrics
# Retains only taxon-level read counts for analyses


beta_d_b <- new_bacteria_s %>%
  dplyr::select(
    Log_ID,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()

Beta_perm_b <- merge(x=bacteria,y=beta_d_b,
                      by=c("Log_ID"))


Bacterial_data <-subset(Beta_perm_b, select = -c(Log_ID,Location,Sex,J_Height, Avg_Crown,Avg_circ,N_stems))

# Calculate Bray–Curtis dissimilarity among bacterial communities
# This metric quantifies differences in community composition
# based on relative abundances of taxa between samples
B.dis_b <- vegdist(Bacterial_data, method = "bray")


# Perform Principal Coordinates Analysis (PCoA) on the Bray–Curtis distance matrix
# PCoA reduces multivariate community data into a low-dimensional space
# that summarizes patterns of similarity among samples
PCoA_B <- cmdscale(B.dis_b, k = 2, eig = TRUE)

# Extract eigenvalues from the PCoA results
# Eigenvalues represent the amount of variation explained by each axis
# and are used to calculate the percentage of variance explained
eig_values_b <- PCoA_B$eig

# Calculate percentage explained by each axis
variance_explained_B <- eig_values_b / sum(eig_values_b)

# Print percentages for the first two axes
percent_explained_B <- round(variance_explained_B[1:2] * 100, 2)
names(percent_explained_B) <- c("PCoA1", "PCoA2")
percent_explained_B

# Extract scores for the first two PCoA axes
# These coordinates represent each sample's position
# in multivariate community composition space
LocSex_PCoA1_B <- PCoA_B$points[,1]
LocSex_PCoA2_B <- PCoA_B$points[,2]

# Combine PCoA coordinates with bacterial abundance data
# This creates a single data frame 
PCoA_plot_B <- cbind(Bacterial_data, LocSex_PCoA1_B, LocSex_PCoA2_B)

# Recover Sex, Location for plotting
PCoA_plot_B$Log_ID <-paste(Beta_perm_b$Log_ID)
PCoA_plot_B$Sex <-paste(Beta_perm_b$Sex)
PCoA_plot_B$Location<-paste(Beta_perm_b$Location)
PCoA_plot_B$J_Height<-paste(Beta_perm_b$J_Height)
PCoA_plot_B$Avg_Crown<-paste(Beta_perm_b$Avg_Crown)
PCoA_plot_B$Avg_circ<-paste(Beta_perm_b$Avg_circ)
PCoA_plot_B$N_stems<-paste(Beta_perm_b$N_stems)


PCoA_plot_B <- PCoA_plot_B %>%
  mutate(
    J_Height = as.numeric(as.character(J_Height)),
    Avg_Crown = as.numeric(as.character(Avg_Crown)),
    Avg_circ = as.numeric(as.character(Avg_circ)),
    N_stems = as.numeric(as.character(N_stems))
  )


# Test for multivariate homogeneity of dispersion (PERMDISP)
# This tests whether group variances (spread) are equal across Sex × Site groups

# Create interaction group (Sex + Location)
PCoA_plot_B$Group <- interaction(PCoA_plot_B$Sex,
                                 PCoA_plot_B$Location)

# Dispersion test for Sex
disp_sex <- betadisper(B.dis_b, PCoA_plot_B$Sex)

# Permutation test
permutest(disp_sex, permutations = 999)

# pairwise comparisons
TukeyHSD(disp_sex)

# Dispersion test for Location

disp_location <- betadisper(B.dis_b, PCoA_plot_B$Location)

permutest(disp_location, permutations = 999)

TukeyHSD(disp_location)

# Dispersion test for interaction

disp_group <- betadisper(B.dis_b, PCoA_plot_B$Group)

permutest(disp_group, permutations = 999)

TukeyHSD(disp_group)




#Bacteria PCoA

# Create a Principal Coordinates Analysis (PCoA) ordination plot
# showing patterns in bacterial community composition
# Samples are colored by plant sex and shaped by sampling site
# Ellipses represent confidence regions for each Sex × Site group



# Perform PERMANOVA to test for differences in bacterial community composition
# as a function of plant sex, sampling location, and their interaction
# using Bray–Curtis dissimilarity and permutation-based significance testing



permanova_result_bst <- adonis2(B.dis_b ~ Sex * Location + N_stems, data = PCoA_plot_B,
                                permutations = 999)
permanova_result_bst



# Fungal beta diversity: data preparation

# Subset abundance-only matrix for Bray–Curtis analysis

beta_f_b <- new_fungi_s %>%
  dplyr::select(
    Log_ID,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems
  ) %>%
  na.omit()

Beta_perm_f <- merge(x=fungi,y=beta_f_b,
                     by=c("Log_ID"))

B.Fungidata <-subset(Beta_perm_f, select = -c(Log_ID,Location,Sex, 
                                              J_Height, Avg_Crown, Avg_circ, N_stems))

# Fungal beta diversity: Bray–Curtis PCoA

# Calculate Bray–Curtis dissimilarity matrix
B.dist_F <- vegdist(B.Fungidata, method = "bray")

# Perform Principal Coordinates Analysis (PCoA)
PCoA_F <- cmdscale(B.dist_F, k = 2, eig = TRUE)
# Extract eigenvalues
eig_values <- PCoA_F$eig

# Calculate percentage explained by each axis
variance_explained_F <- eig_values / sum(eig_values)

# Print percentages for the first two axes
percent_explained_F <- round(variance_explained_F[1:2] * 100, 2)
names(percent_explained_F) <- c("PCoA1", "PCoA2")
percent_explained_F

# Extract PCoA coordinates
LocSex_PCoA1_F <- PCoA_F$points[,1]
LocSex_PCoA2_F <- PCoA_F$points[,2]
# Combine ordination scores with abundance data
PCoA_plot_F <- cbind(B.Fungidata, LocSex_PCoA1_F, LocSex_PCoA2_F)

# Recover Sex,Location and log_ID for plotting
PCoA_plot_F$Log_ID <-paste(Beta_perm_f$Log_ID)
PCoA_plot_F$Sex <-paste(Beta_perm_f$Sex)
PCoA_plot_F$Location<-paste(Beta_perm_f$Location)
PCoA_plot_F$J_Height<-paste(Beta_perm_f$J_Height)
PCoA_plot_F$Avg_Crown<-paste(Beta_perm_f$Avg_Crown)
PCoA_plot_F$Avg_circ<-paste(Beta_perm_f$Avg_circ)
PCoA_plot_F$N_stems<-paste(Beta_perm_f$N_stems)


PCoA_plot_Fn <- PCoA_plot_F %>%
  mutate(
    J_Height = as.numeric(as.character(J_Height)),
    Avg_Crown = as.numeric(as.character(Avg_Crown)),
    Avg_circ = as.numeric(as.character(Avg_circ)),
    N_stems = as.numeric(as.character(N_stems))
  )



# Test for multivariate homogeneity of dispersion (PERMDISP)
# This tests whether group variances (spread) are equal across Sex × Site groups

# Create interaction group (Sex + Location)
PCoA_plot_F$Group <- interaction(PCoA_plot_F$Sex,
                                 PCoA_plot_F$Location)

# Dispersion test for Sex
disp_sexf <- betadisper(B.dist_F, PCoA_plot_F$Sex)

# Permutation test
permutest(disp_sexf, permutations = 999)

# pairwise comparisons
TukeyHSD(disp_sexf)

# Dispersion test for Location

disp_locationf <- betadisper(B.dist_F, PCoA_plot_F$Location)

permutest(disp_locationf, permutations = 999)

TukeyHSD(disp_locationf)

# Dispersion test for interaction

disp_groupf <- betadisper(B.dist_F, PCoA_plot_F$Group)

permutest(disp_groupf, permutations = 999)

TukeyHSD(disp_groupf)




# PERMANOVA: fungal community composition

# Test effects of Sex, Location, and interaction

permanova_result_fst<- adonis2(B.dist_F ~ Sex * Location + N_stems, data = PCoA_plot_Fn,
                                permutations = 999)
permanova_result_fst




# Arthropod beta diversity


# Remove extra columns
arthropod_wide <-subset(arthropods2, select = -c(
  Collection_Date,Analyzer,Analysis_Date,Count, Mass))


# Spread into wide format and fill with 0
Art_Beta_data<-pivot_wider(arthropod_wide,names_from = Taxa,values_from = Countperweight, 
                           values_fill = 0)


# Remove no needed columns

new_arthro_sn<-new_arthro_s%>%
  select(-c(Sex, Location,Shannon))


Beta_perm_art <- merge(x=Art_Beta_data,y=new_arthro_sn,
                     by=c("Log_ID"))

Art_Beta_data_n <-subset(Beta_perm_art, select = -c(Log_ID,Location,Sex,J_Height,Avg_Crown,Avg_circ,N_stems))


# Calculate Bray–Curtis dissimilarity matrix
Art_B_dist <- vegdist(Art_Beta_data_n, method = "bray")

# Perform Principal Coordinates Analysis (PCoA)
PCoA <- cmdscale(Art_B_dist, k = 2, eig = TRUE)
# Extract eigenvalues
eig_values <- PCoA$eig

# Calculate percentage explained by each axis
variance_explained <- eig_values / sum(eig_values)

# Print percentages for the first two axes
percent_explained <- round(variance_explained[1:2] * 100, 2)
names(percent_explained) <- c("PCoA1", "PCoA2")
percent_explained
LocSex_PCoA1 <- PCoA$points[,1]
LocSex_PCoA2 <- PCoA$points[,2]
PCoA_plot <- cbind(Art_Beta_data_n , LocSex_PCoA1, LocSex_PCoA2)

# Recover Sex, Location and Log_ID

PCoA_plot$Log_ID <-paste(Beta_perm_art$Log_ID)
PCoA_plot$Sex <-paste(Beta_perm_art$Sex)
PCoA_plot$Location<-paste(Beta_perm_art$Location)
PCoA_plot$J_Height<-paste(Beta_perm_art$J_Height)
PCoA_plot$Avg_Crown<-paste(Beta_perm_art$Avg_Crown)
PCoA_plot$Avg_circ<-paste(Beta_perm_art$Avg_circ)
PCoA_plot$N_stems<-paste(Beta_perm_art$N_stems)


PCoA_plot_art <- PCoA_plot %>%
  mutate(
    J_Height = as.numeric(as.character(J_Height)),
    Avg_Crown = as.numeric(as.character(Avg_Crown)),
    Avg_circ = as.numeric(as.character(Avg_circ)),
    N_stems = as.numeric(as.character(N_stems))
  )

PCoA_plot_art$Location <- with(PCoA_plot_art, 
                             ifelse(grepl("TC", PCoA_plot_art$Log_ID), "CTC",
                                    ifelse(grepl("H", PCoA_plot_art$Log_ID), "BH", 'AWi')))




# Test for multivariate homogeneity of dispersion (PERMDISP)
# This tests whether group variances (spread) are equal across Sex × Site groups

# Create interaction group (Sex + Location)
PCoA_plot$Group <- interaction(PCoA_plot$Sex,
                                 PCoA_plot$Location)

# Dispersion test for Sex
disp_sex.ar <- betadisper(Art_B_dist, PCoA_plot$Sex)

# Permutation test
permutest(disp_sex.ar, permutations = 999)

# pairwise comparisons
TukeyHSD(disp_sex.ar)

# Dispersion test for Location

disp_location.ar <- betadisper(Art_B_dist, PCoA_plot$Location)

permutest(disp_location.ar, permutations = 999)

TukeyHSD(disp_location.ar)

# Dispersion test for interaction

disp_group.ar <- betadisper(Art_B_dist, PCoA_plot$Group)

permutest(disp_group.ar, permutations = 999)

TukeyHSD(disp_group.ar)


# PERMANOVA: Arthropod community composition

# Test effects of Sex, Location, and interaction

permanova_result_artst <- adonis2(Art_B_dist ~ Sex * Location + N_stems, data = PCoA_plot_art, 
                                permutations = 999)
permanova_result_artst




# Figure #2 Beta Diversity- PCOA multipanel plot

# Bacteria PCoA
plot_beta_b <- ggplot() +
  geom_point(
    data = PCoA_plot_B,
    aes(
      x = LocSex_PCoA1_B,
      y = LocSex_PCoA2_B,
      color = Location,       
      shape = Sex   
    ),
    size = 1.6, alpha = 0.5
  ) +
  stat_ellipse(
    data = PCoA_plot_B,
    aes(
      x = LocSex_PCoA1_B,
      y = LocSex_PCoA2_B,
      color = Location,
      linetype = Sex,
      group = interaction(Sex, Location)
    ),
    type = "t", size = 0.8, alpha = 0.9
  ) +
  scale_color_manual(
    values = c("BH" = "#0072B2", "CTC" = "#E69F00", "AWi" = "#009E73"),
    labels = c("BH" = "PERL", "CTC" = "BTNP", "AWi" = "LBJWC")
  ) +
  scale_shape_manual(values = c("F" = 1, "M" = 2)) +
  scale_linetype_manual(values = c("F" = "dashed", "M" = "solid")) +
  labs(
    title = "Bacteria",
    x = paste0("PCoA1 (", percent_explained_B["PCoA1"], "%)"),
    y = paste0("PCoA2 (", percent_explained_B["PCoA2"], "%)"),
    color = "Site", shape = "Sex", linetype = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),   
    axis.text = element_text(size = 10),     
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(0.9, "lines"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2, size = 12)
  ) + guides(shape = guide_legend(override.aes = list(size = 4))) +
  annotate(
    "text",
    x = 0.2,
    y = 0.45,
    label = "*",
    size = 6,
    fontface = "bold"
  ) +
  coord_fixed()

# Fungi PCoA
plot_beta_f <- ggplot() +
  geom_point(
    data = PCoA_plot_F,
    aes(
      x = LocSex_PCoA1_F,
      y = LocSex_PCoA2_F,
      color = Location,       
      shape = Sex   
    ),
    size = 1.6, alpha = 0.5
  ) +
  stat_ellipse(
    data = PCoA_plot_F,
    aes(
      x = LocSex_PCoA1_F,
      y = LocSex_PCoA2_F,
      color = Location,
      linetype = Sex,
      group = interaction(Sex, Location)
    ),
    type = "t", size = 0.8, alpha = 0.9
  ) +
  scale_color_manual(values = c("BH" = "#0072B2", "CTC" = "#E69F00", "AWi" = "#009E73")) +
  scale_shape_manual(values = c("F" = 1, "M" = 2)) +
  scale_linetype_manual(values = c("F" = "dashed", "M" = "solid")) +
  labs(
    title = "Fungi",
    x = paste0("PCoA1 (", percent_explained_F["PCoA1"], "%)"),
    y = paste0("PCoA2 (", percent_explained_F["PCoA2"], "%)"),
    color = "Site", shape = "Sex", linetype = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),   
    axis.text = element_text(size = 10),     
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(0.9, "lines"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2, size = 12)
  ) + annotate(
    "text",
    x = mean(range(PCoA_plot_F$LocSex_PCoA1_F, na.rm = TRUE)),
    y = max(PCoA_plot_F$LocSex_PCoA2_F, na.rm = TRUE) * 1.08,
    label = "***",
    size = 6,
    fontface = "bold"
  ) +
  coord_fixed()

# Arthropod PCoA
plot_beta_a <- ggplot() +
  geom_point(
    data = PCoA_plot,
    aes(
      x = LocSex_PCoA1,
      y = LocSex_PCoA2,
      color = Location,       
      shape = Sex   
    ),
    size = 1.6, alpha = 0.5
  ) +
  stat_ellipse(
    data = PCoA_plot,
    aes(
      x = LocSex_PCoA1,
      y = LocSex_PCoA2,
      color = Location,
      linetype = Sex,
      group = interaction(Sex, Location)
    ),
    type = "t", size = 0.8, alpha = 0.9
  ) +
  scale_color_manual(values = c("BH" = "#0072B2", "CTC" = "#E69F00", "AWi" = "#009E73")) +
  scale_shape_manual(values = c("F" = 1, "M" = 2)) +
  scale_linetype_manual(values = c("F" = "dashed", "M" = "solid")) +
  labs(
    title = "Arthropod",
    x = paste0("PCoA1 (", percent_explained["PCoA1"], "%)"),
    y = paste0("PCoA2 (", percent_explained["PCoA2"], "%)"),
    color = "Site", shape = "Sex", linetype = "Sex"
  ) +
  theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),   
    axis.text = element_text(size = 10),     
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 12),   
    legend.key.size = unit(0.9, "lines"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, vjust = 2, size = 12)
  ) +
  coord_fixed()

# Combine three PCoAs vertically
beta_div <- plot_grid(
  plot_beta_b, plot_beta_f, plot_beta_a,
  nrow = 3,
  labels = c("A", "B", "C"),
  align = "hv",
  label_x = 0.4,
  label_y = 0.99,   # adjust since plots are stacked
  label_size = 12
)

legend_pcoa <- get_legend(
  plot_beta_b +
    theme(
      legend.position = "right",
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0)
    )
)

beta_with_legend <- ggdraw() +
  draw_plot(beta_div, x = 0, y = 0, width = 0.85, height = 1) +
  draw_plot(legend_pcoa, x = 0.58, y = 0.25, width = 0.15, height = 0.5)



beta_with_legend

# save plots

ggsave("Pcoanew_figure.tiff",beta_with_legend, width = 12, height = 10, device='tiff', dpi=600)

ggsave("Pcoa_figure.png",beta_with_legend, width = 12, height = 10, device='png', dpi=300)





# Unique bacteria taxa
  
# Use mutate() to modify the existing 'Location' column
# Replace shorthand codes with full descriptive names
bacteria_data_combined <- bacteria %>%
mutate(Location = recode(Location,
                           "AWi" = "Wildflower",
                           "BH" = "Huntsville",
                           "CTC" = "Turkey Creek"))


# Select bacterial taxa that are present in at least 2 samples per site
taxa_to_keep_b <- bacteria_data_combined %>%
  
# Reshape the data set from wide to long format
# Keep the columns Log_ID, Location, and Sex as identifiers
# Gather all other columns (assumed to be taxa) into two columns:
# "Taxa" for the taxon name
# "Abundance" for the counts
pivot_longer(cols = -c(Log_ID, Location, Sex), 
               names_to = "Taxa", 
               values_to = "Abundance") %>%
# Keep only rows where the taxon is present in the sample (Abundance > 0)
filter(Abundance > 0) %>%
  
# Remove duplicate rows for the same sample and taxon
# (ensures each sample x taxon combination is unique)
distinct(Log_ID, Location, Taxa) %>%
  
# Group by site and taxon to count the number of samples
group_by(Location, Taxa) %>%
# Summarize by counting how many samples each taxon appears in per site
summarize(Count = n(), .groups = "drop") %>%
# Keep only taxa that are present in at least 2 samples per location
filter(Count >= 2)


# Extract the names of the taxa we decided to keep
selected_taxa_b <-unique(taxa_to_keep_b$Taxa)

# Define metadata columns that should be preserved
metadata_cols_b <- c("Log_ID", "Location", "Sex")

# Extract only the abundance data (all taxa columns) from the original dataset
abundance_data_b <- bacteria_data_combined[, !(names(bacteria_data_combined) %in% metadata_cols_b)]

# Keep only the columns corresponding to the selected taxa
abundance_data_filtered_b <- abundance_data_b[, names(abundance_data_b) %in% selected_taxa_b]

# Combine the metadata columns back with the filtered abundance data
filtered_datab<- cbind(bacteria_data_combined[metadata_cols_b], abundance_data_filtered_b)

# Function to summarize bacterial taxa by sex and location
get_summary_bacteria <- function(data) {
  data%>%
    
# Reshape data from wide to long format
# Keep metadata columns (Log_ID, Location, Sex)
# Gather all taxa columns into 'Taxa' and 'Abundance'
pivot_longer(cols = -c(Log_ID, Location, Sex), 
                 names_to = "Taxa", 
                 values_to = "Abundance") %>%
# Count the number of samples in which each taxon is present (>0)
# Group by Location, Taxa, and Sex
group_by(Location, Taxa, Sex)%>%
summarize(Count = sum(Abundance > 0), .groups = "drop") %>%
# Reshape back to wide format: separate columns for female (F) and male (M)
pivot_wider(names_from = Sex, values_from = Count)%>%
    
# Assign presence categories based on which sex the taxon occurs in
mutate(
      presence = case_when(
        F > 0 & M > 0 ~ "shared", # Taxa present in both sexes
        F > 0 & M == 0 ~ "female-only", # Taxa present only in females
        F == 0 & M > 0 ~ "male-only"     # Taxa present only in males
      )
    ) %>%
# Count the number of taxa in each presence category for each location
group_by(Location, presence) %>%
summarize(n_taxa = n(), .groups = "drop")
}
# Observed data summary
observed_summary_b <- get_summary_bacteria(filtered_datab) %>%
# Mark observed data with perm = 0
  mutate(perm = 0)



# Unique fungal OTUs

# Filter fungal taxa present in at least 2 samples per site

fungi_data_combined <- fungi %>%
  mutate(Location = recode(Location,
                           "Wi" = "Wildflower",
                           "H" = "Huntsville",
                           "TC" = "Turkey Creek"))


taxa_to_keep_f <- fungi_data_combined %>%
# Convert wide abundance table to long format
pivot_longer(cols = -c(Log_ID, Location, Sex), 
               names_to = "Taxa", 
               values_to = "Abundance") %>%
# Keep only taxa that are present
filter(Abundance > 0) %>%
# Remove duplicate occurrences within samples
distinct(Log_ID, Location, Taxa) %>%
# Count how many samples contain each taxon per site
group_by(Location, Taxa) %>%
summarize(Count = n(), .groups = "drop") %>%
# Retain taxa present in ≥ 2 samples per site
filter(Count >= 2)


# Extract names of retained taxa
selected_taxa_f <-unique(taxa_to_keep_f$Taxa)

# Separate metadata from abundance matrix

# Metadata columns
metadata_cols_f <- c("Log_ID", "Location", "Sex")
# Extract abundance-only matrix
abundance_data_f <- fungi_data_combined[, !(names(fungi_data_combined) %in% metadata_cols_f)]
# Keep only filtered taxa
abundance_data_filtered_f <- abundance_data_f[, names(abundance_data_f) %in% selected_taxa_f]

# Recombine metadata and filtered abundance data
filtered_dataf <- cbind(fungi_data_combined[metadata_cols_f], abundance_data_filtered_f)

# Function to summarize sex-specific taxon occurrence
# Counts how many samples contain each taxon
# Determines whether taxa are shared, female-only, or male-only
# Returns the number of taxa in each category
get_summary_fungi <- function(data) {
  data%>%
# Convert to long format
pivot_longer(cols = -c(Log_ID, Location, Sex), 
                 names_to = "Taxa", 
                 values_to = "Abundance") %>%
# Count presence of each taxon by sex
group_by(Location, Taxa, Sex)%>%
summarize(Count = sum(Abundance > 0), .groups = "drop") %>%
# Convert back to wide format (F vs M)
pivot_wider(names_from = Sex, values_from = Count)%>%
# Classify taxa by sex specificity
mutate(
      presence = case_when(
        F > 0 & M > 0 ~ "shared",
        F > 0 & M == 0 ~ "female-only",
        F == 0 & M > 0 ~ "male-only"
      )
    ) %>%
# Count taxa in each category per site
group_by(Location, presence) %>%
summarize(n_taxa = n(), .groups = "drop")
}

# Observed (non-permuted) summary
observed_summary_f <- get_summary_fungi(filtered_dataf) %>%
# Label as observed
mutate(perm = 0)




# Unique arthropod taxa

# Removes non-essential metadata columns
arthropods3<-subset(arthropods, select = -c(Collection_Date,
                                            Analyzer,Analysis_Date))

#  Converts the dataset from long to wide format
arthropods3_n<-pivot_wider(arthropods3,names_from = Taxa,
                         values_from = Count, values_fill = 0)



# Filter the arthropod taxa present in at least 2 samples


taxa_to_keep_art <-  arthropods3_n %>%
# Convert wide abundance table to long format
pivot_longer(cols = -c(Log_ID, Location, Sex), 
               names_to = "Taxa", 
               values_to = "Abundance") %>%
# Keep only taxa that are present
filter(Abundance > 0) %>%
# Remove duplicate occurrences within samples
distinct(Log_ID, Location, Taxa) %>%
group_by(Location, Taxa) %>%
# Count how many samples contain each taxon per site
summarize(Count = n(), .groups = "drop") %>%
# Retain taxa present in ≥ 2 samples per site
filter(Count >= 2)


# Extract names of retained taxa
selected_taxa_ar <-(taxa_to_keep_art$Taxa)

# Separate metadata from abundance matrix

# Metadata columns
metadata_cols_art <- c("Log_ID", "Location", "Sex")

# Extract abundance-only matrix
abundance_data_art <- arthropods3_n[, !(names(arthropods3_n) %in% metadata_cols_art)]
# Keep only filtered taxa
abundance_data_filtered_art <- abundance_data_art[, names(abundance_data_art) %in% selected_taxa_ar]

# Recombine metadata and filtered abundance data
filtered_data_art <- cbind(arthropods3_n[metadata_cols_art], abundance_data_filtered_art)

# Function to summarize sex-specific taxon occurrence
# Counts how many samples contain each taxon
# Determines whether taxa are shared, female-only, or male-only
# Returns the number of taxa in each category

get_summary <- function(data) {
  data %>%
# Convert to long format
pivot_longer(cols = -c(Log_ID, Location, Sex), 
                 names_to = "Taxa", 
                 values_to = "Abundance") %>%
# Count presence of each taxon by sex
group_by(Location, Taxa, Sex) %>%
summarize(Count = sum(Abundance > 0), .groups = "drop") %>%
# Convert back to wide format (F vs M)
pivot_wider(names_from = Sex, values_from = Count, values_fill = 0) %>%
# Classify taxa by sex specificity
mutate(
      presence = case_when(
        F > 0 & M > 0 ~ "shared",
        F > 0 & M == 0 ~ "female-only",
        F == 0 & M > 0 ~ "male-only"
      )
    ) %>%
# Count taxa in each category per site
group_by(Location, presence) %>%
summarize(n_taxa = n(), .groups = "drop")
}

# Observed (non-permuted) summary
observed_summary <- get_summary(filtered_data_art) %>%
# Label as observed
mutate(perm = 0)




# Figure #3 Venn Diagram Multipanel

make_euler <- function(
    male_only,
    female_only,
    shared,
    letter,
    stars = NULL   # list of stars with x/y positions
) {
  
# Use exclusive areas
  areas <- c(
    Male = male_only,
    Female = female_only,
    "Male&Female" = shared
  )
  
  fit <- euler(areas)
  
  venn_grob <- plot(
    fit,
    fills = list(
      fill = c(
        Male = "darkcyan",
        Female = "tomato",
        "Male&Female" = "beige"
      ),
      alpha = 1
    ),
    edges = list(
      col = "black",
      lwd = 1.5
    ),
    quantities = list(
      fontsize = 10,
      fontface = "bold"
    ),
    labels = list(
      labels = c("", ""),
      fontsize = 0
    ),
    return = "grob"
  )
  
# Panel letter
  letter_grob <- textGrob(
    letter,
    x = unit(0, "npc") + unit(6, "pt"),
    y = unit(1, "npc") - unit(6, "pt"),
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  
# Store grobs
  grob_list <- gList(venn_grob, letter_grob)
  
# Add stars 
  if (!is.null(stars)) {
    
    for (i in seq_along(stars)) {
      
      star_grob <- textGrob(
        "*",
        x = stars[[i]]$x,
        y = stars[[i]]$y,
        gp = gpar(fontsize = 25, fontface = "bold")
      )
      
      grob_list <- gList(grob_list, star_grob)
    }
  }
  
  gTree(children = grob_list)
}


# BACTERIA

venn_B_W  <- make_euler(949, 1010, 5037, "A")
venn_B_H  <- make_euler(1268, 1373, 6723, "B")
venn_B_TC <- make_euler(2827, 598, 4078, "C")


# FUNGI

# Wildflower = LBJWC
# One star near female-only circle
# One star near male-only circle

venn_F_W <- make_euler(
  684, 372, 2709, "D",
  stars = list(
    list(x = 0.27, y = 0.94),  # male-only star
    list(x = 0.70, y = 0.94)   # female-only star
  )
)

# Huntsville = PERL
# One star near female-only circle

venn_F_H <- make_euler(
  845, 968, 2940, "E",
  stars = list(
    list(x = 0.70, y = 0.98)   # female-only star
  )
)

venn_F_TC <- make_euler(
  1213, 296, 1131, "F"
)


# ARTHROPODS

venn_A_W <- make_euler(
  1, 2, 24, "G"
)

venn_A_H <- make_euler(
  2, 2, 21, "H"
)

venn_A_TC <- make_euler(
  4, 4, 8, "I"
)


# ROW LABELS

row_label <- function(text) {
  textGrob(
    text,
    gp = gpar(fontsize = 13, fontface = "bold")
  )
}


# COLUMN LABELS

col_labels <- arrangeGrob(
  textGrob("LBJWC", gp = gpar(fontsize = 13)),
  textGrob("PERL",  gp = gpar(fontsize = 13)),
  textGrob("BTNP",  gp = gpar(fontsize = 13)),
  ncol = 3
)


# LEGEND

legend_grob <- gTree(children = gList(
  
  circleGrob(
    x = 0.38,
    y = 0.5,
    r = 0.15,
    gp = gpar(fill = "darkcyan", col = "black")
  ),
  
  textGrob(
    "Male-only",
    x = 0.40,
    y = 0.5,
    just = "left"
  ),
  
  circleGrob(
    x = 0.48,
    y = 0.5,
    r = 0.15,
    gp = gpar(fill = "tomato", col = "black")
  ),
  
  textGrob(
    "Female-only",
    x = 0.50,
    y = 0.5,
    just = "left"
  ),
  
  circleGrob(
    x = 0.60,
    y = 0.5,
    r = 0.15,
    gp = gpar(fill = "beige", col = "black")
  ),
  
  textGrob(
    "Shared",
    x = 0.62,
    y = 0.5,
    just = "left"
  )
))


# ARRANGE ROWS

row1 <- arrangeGrob(
  row_label("Bacteria"),
  arrangeGrob(
    venn_B_W,
    venn_B_H,
    venn_B_TC,
    ncol = 3
  ),
  ncol = 1,
  heights = c(1.2, 10)
)

row2 <- arrangeGrob(
  row_label("Fungi"),
  arrangeGrob(
    venn_F_W,
    venn_F_H,
    venn_F_TC,
    ncol = 3
  ),
  ncol = 1,
  heights = c(1.2, 10)
)

row3 <- arrangeGrob(
  row_label("Arthropods"),
  arrangeGrob(
    venn_A_W,
    venn_A_H,
    venn_A_TC,
    ncol = 3
  ),
  ncol = 1,
  heights = c(1.2, 10)
)


# MAIN GRID

main_grid <- arrangeGrob(
  row1,
  row2,
  row3,
  ncol = 1
)

# FINAL PLOT

grid.newpage()

final_plot <- grid.arrange(
  legend_grob,
  main_grid,
  col_labels,
  nrow = 3,
  heights = c(2, 30, 1.5)
)


# Save
tiff("venn_eulerr_multipanel.tiff",
     width = 11, height = 10,
     units = "in", res = 600,
     compression = "lzw")

grid.draw(final_plot)
dev.off()

png("venn_eulerr_multipanelfinal.png",
    width = 11,
    height = 10,
    units = "in",
    res = 600,
    type = "cairo")  

grid.draw(final_plot)
dev.off()


# Question #2: How do male and female plants differ in soil biogeochemical properties known to affect soil communities?


# SOIL BIOGEOCHEMICAL PROPERTIES


# Our sites represent a precipitation gradient. Driest to wettest from left to right
# so to plot the sites in that order with assigned letters to them

Ilex_traits$Location <- with(Ilex_traits, 
                             ifelse(grepl("TC", Ilex_traits$Log_ID), "CTC",
                                    ifelse(grepl("H", Ilex_traits$Log_ID), "BH", 'AWi')))
# Soil moisture

soil_moisture <- Ilex_traits %>%
  filter(!is.na(Moisture_VWC))


new_soil_m<-soil_moisture%>%
  select(c(Log_ID,Moisture_VWC,Sex,Location))


new_ilex_bio<-Ilex_traits%>%
  select(c(Log_ID,J_Height,Avg_Crown,Avg_circ,N_stems))

new_biogeo <- merge(x=new_soil_m,y=new_ilex_bio,
                        by=c("Log_ID"))

biogeo_model_data <- new_biogeo %>%
  dplyr::select(
    Log_ID,
    Sex,
    Location,
    J_Height,
    Avg_Crown,
    Avg_circ,
    N_stems,
    Moisture_VWC
  ) %>%
  na.omit()

# ANOVA: Soil moisture + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelSms<-lm(Moisture_VWC ~ Sex * Location + N_stems, data = soil_moisture)
# Type III ANOVA (for unbalanced designs)
Anova(modelSms, type = 3)
# Summary of model coefficients
summary(modelSms)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelSms)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelSms, pairwise ~ Sex * Location, adjust = "tukey")


# pH

new_ilex_biogeo<-biogeo_model_data%>%
  select(c(Log_ID,J_Height,Avg_Crown,Avg_circ,N_stems))

new_biogeon <- merge(x=new_ilex_biogeo,y=biogeochemical,
                    by=c("Log_ID"))




# ANOVA: pH + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelPs<-lm(pH ~ Sex * Location + N_stems, data = new_biogeon)
# Type III ANOVA (for unbalanced designs)
Anova(modelPs, type = 3)
# Summary of model coefficients
summary(modelPs)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelPs)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelPs, pairwise ~ Sex * Location, adjust = "tukey")




# Our sites represent a precipitation gradient. Driest to wettest from left to right
# so to plot the sites in that order with assigned letters to them
org_matter$Location <- with(org_matter, 
                            ifelse(grepl("TC", org_matter$Log_ID), "CTC",
                                   ifelse(grepl("H", org_matter$Log_ID), "BH", 'AWi')))

new_org_matter<-org_matter%>%
  select(c(Log_ID,TotalN))

biogeo_n <- merge(x=new_org_matter,y=new_biogeon,
                     by=c("Log_ID"))


# Total nitrogen

# ANOVA: Total nitrogen + Number of stems

# Fit linear model testing Sex, Location, and interaction
modelNis<-lm(TotalN ~ Sex * Location + N_stems, data = biogeo_n)
# Type III ANOVA (for unbalanced designs)
Anova(modelNis, type = 3)
# Summary of model coefficients
summary(modelNis)
# Diagnostic plots (normality, homoscedasticity, leverage)
plot(modelNis)
# Post-hoc pairwise comparisons (Tukey-adjusted)
emmeans(modelNis, pairwise ~ Sex * Location, adjust = "tukey")




# Our sites represent a precipitation gradient. Driest to wettest from left to right
# so to plot the sites in that order with assigned letters to them


# AWi=Wildflower;BH=Huntsville;CTC=Turkey Creek
microbial_mass$Location <- with(microbial_mass, 
                                ifelse(grepl("TC", microbial_mass$Log_ID), "CTC",
                                       ifelse(grepl("H", microbial_mass$Log_ID), "BH", 'AWi')))


new_microbial_m<-microbial_mass%>%
  select(c(Log_ID,MBC))

biogeo_micro <- merge(x=new_microbial_m,y=new_biogeon,
                  by=c("Log_ID"))




# ANOVA: Microbial Biomass Carbon + Number of stems

modelMbs<-lm(MBC ~ Sex * Location + N_stems, data = biogeo_micro)
Anova(modelMbs, type = 3)
summary(modelMbs)
plot(modelMbs)
emmeans(modelMbs, pairwise ~ Sex * Location, adjust = "tukey")




# Figure #4 Biogeochemical properties

# Biogeochemical boxplot multipanel


base_theme <- theme_classic(base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "lines"),
    legend.position = "none",
    plot.margin = unit(c(5, 5, 5, 30), "pt")
  )

# Significance annotations

# pH: ** at BTNP (CTC)
ph_star <- data.frame(
  Location = "CTC",
  y = max(biogeochemical$pH, na.rm = TRUE) * 1.08,
  label = "**"
)

p_s_m <- ggplot(soil_moisture, aes(x = Location, y = Moisture_VWC, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  labs(x = "Site", y = "Soil Moisture VWC %", fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male"),
    breaks = c("F", "M")
  ) +
  scale_x_discrete(labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")) +
  base_theme

p_p_h <- ggplot(biogeochemical, aes(x = Location, y = pH, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  geom_text(
    data = ph_star,
    aes(x = Location, y = y, label = label),
    inherit.aes = FALSE,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Site", y = "pH", fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male"),
    breaks = c("F", "M")
  ) +
  scale_x_discrete(labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")) +
  base_theme

p_t_n <- ggplot(org_matter, aes(x = Location, y = TotalN, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  labs(x = "Site", y = "Total Nitrogen-mg g⁻¹ soil)", fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male"),
    breaks = c("F", "M")
  ) +
  scale_x_discrete(labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")) +
  base_theme

p_m_b <- ggplot(microbial_mass, aes(x = Location, y = MBC, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85),
    alpha = 0.6, size = 2, shape = 21, stroke = 0.4, color = "black"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(clip = "off") +
  labs(x = "Site",
       y = "MBC (mg C g⁻¹ dry soil)",
       fill = "Sex") +
  scale_fill_manual(
    values = c("F" = "tomato", "M" = "darkcyan"),
    labels = c("F" = "Female", "M" = "Male"),
    breaks = c("F", "M")
  ) +
  scale_x_discrete(
    labels = c("AWi" = "LBJWC", "BH" = "PERL", "CTC" = "BTNP")
  ) +
  base_theme


# Function to format panels
format_panel <- function(p, ylab, show_x = FALSE) {
  p +
    labs(y = ylab) +
    theme(
      axis.title.x = if (show_x) element_text(size = 12) else element_blank(),
      axis.title.y = element_text(size = 12),
      axis.text.x = if (show_x) element_text(size = 12) else element_blank(),
      axis.text.y = element_text(size = 12)
    )
}

# Apply formatting
plot_s_m2 <- format_panel(p_s_m, "Soil moisture  VWC %", show_x = FALSE)
plot_p_h2 <- format_panel(p_p_h, "pH", show_x = FALSE)
plot_t_n2 <- format_panel(p_t_n, "Total Nitrogen (mg g⁻¹ soi)l", show_x = TRUE)
plot_m_b2 <- format_panel(p_m_b, "MBC (mg C g⁻¹ dry soil)", show_x = TRUE)

# Extract one legend (Female on top)
legend_box <- get_legend(
  p_s_m + theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0)
  )
)

# Stack plots
alpha_div <- plot_grid(
  plot_s_m2, plot_p_h2,
  plot_t_n2, plot_m_b2,
  ncol = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 1.00,
  align = "hv"
)

# Combine with legend on the right
plot_biog <- plot_grid(
  alpha_div,
  legend_box,
  ncol = 2,
  rel_widths = c(1, 0.15)
)

# Final plot
plot_biog


# save plots
ggsave("biog_figurenf.tiff",plot_biog, width = 9.5, height = 9, device='tiff', dpi=600)

ggsave("bio_figurefinal.png",plot_biog, width = 9.5, height = 10, device='png', dpi=600)




# Preparing data for structural equation modeling (SEM)



# ARTHROPOD FUNCTIONAL GUILDS; We will use this for SEM analyses

# Assign functional groups to arthropod taxa based on their diet 

functional_arthropods <- arthropod_wide %>%
mutate(Functional_Group = case_when(
 Taxa %in% c("Diplopoda", "Isopoda") ~ "Detritivores",
 Taxa %in% c("Entomobryidae", "Hypogastruridae", "Isotomidae", "Onychiuridae", 
         "Sminthuridae", "Oribatida", "Cryptophagidae") ~ "Fungivores",
 Taxa %in% c("Chilopoda", "Linyphiidae", "Lycosidae", "Mesostigmata", "Gnaphosidae",
           "Araneae_other", "Staphylinidae","Pseudoscorpions") ~ "Predators",
    Taxa %in% c("Hemiptera", "Thysanoptera", "Curculionidae") ~ "Herbivores",
    Taxa %in% c("Formicidae", "Prostigmata", "Blattodea") ~ "Omnivores",
 Taxa  %in% c("Arthropod_other", "Coleoptera_larval","Hymenoptera_other","Lepidoptera_larval","Coleoptera_other") ~ "Others"
  ))

functional_arthropods1<-subset(functional_arthropods, select = (-Taxa))

# Taxa with multiple no certain functional group remains as NA & they are removed here

functional_arthropods1_n<-na.omit(functional_arthropods1)

# Calculate abundance of each arthropod functional group
sum_data <- functional_arthropods1_n %>%
group_by(Log_ID, Functional_Group) %>%
summarize(total_value = sum(Countperweight, na.rm = TRUE), .groups = "drop")

relative_abundance_darthro <- sum_data %>%
  group_by(Log_ID) %>%
  mutate(relative_abundance = total_value / sum(total_value)) %>%
  ungroup()

relative_abundance_darthro_f <-subset(relative_abundance_darthro, select = -c(total_value))


dfun_art_wide <- relative_abundance_darthro_f  %>%
  group_by(Log_ID)  %>%
  pivot_wider(names_from = Functional_Group, values_from = relative_abundance, values_fill = 0)

# Recover Location and sex from Log_ID
dfun_art_func <- relative_abundance_darthro_f %>%
  separate(Log_ID,
           into = c("Location", "Sex"),
           sep = "-",
           remove = FALSE) %>%
  mutate(
    Location = str_remove(Location, "\\d+"),
    Sex = str_remove(Sex, "\\d+")
  )


# Question #3:What are the direct and indirect effects of plant sex on soil food webs?


# Structural Equation Modeling (SEM) ## update names of SEM

# Combine datasets to obtain Sex, biogeochemical, and food web trophic levels needed for SEMs

# Select the needed columns
soil_moisture_n<-soil_moisture %>%
  select(c(Log_ID,Location, Sex, Moisture_VWC))


biogeochemical_n<-biogeochemical %>%
  select(c(Log_ID,pH))

# Add pH data 
sem_d <- merge(x=soil_moisture_n,y=biogeochemical_n,
                  by=c("Log_ID"))

microbial_mass_n<-microbial_mass %>%
  select(c(Log_ID,MBC))

# Add microbial biomass carbon data
sem_da<-merge(x=sem_d,y=microbial_mass_n,
                by=c("Log_ID"))

org_matter_n<-org_matter %>%
  select(c(Log_ID,TotalN))

# Add Total Nitrogen data
sem_dat<-merge(x=sem_da,y=org_matter_n,
                by=c("Log_ID"))



### new data covariates

new_soilm<-new_soil_m %>%
  select(c(Log_ID,Moisture_VWC))

# combine data
sem_dat<-merge(x=new_soilm,y=biogeo_micro,
               by=c("Log_ID"))

sem_datn<-merge(x=sem_dat,y=new_org_matter,
               by=c("Log_ID"))

food_web<-dfun_art_wide %>%
  select(c(Log_ID,Detritivores,Fungivores,Predators))

# Add trophic level abundance data
sem_data<-merge(x=sem_datn,y=food_web,
                by=c("Log_ID"))


# Transform plant sex to levels
sem_data$Sex <- ifelse(sem_data$Sex == "F", 1, 0)



# Transform location to factors
sem_data$Location <- factor(sem_data$Location)



# SEM with location/site as random effect

# Figure #5 SEM - Note: Figure was created using Biorender

#SEM #1

modsem1<- lmer(MBC ~ pH + TotalN + Moisture_VWC + Sex + (1|Location), data = sem_data)
modsem1
anova(modsem1, ddf='Kenward-Roger')
plot(modsem1)

modsem2<- lmer(Fungivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + MBC, data = sem_data)
modsem2
anova(modsem2, ddf='Kenward-Roger')
plot(modsem2)

modsem3<- lmer(Detritivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location), data = sem_data)
modsem3
anova(modsem3, ddf='Kenward-Roger')
plot(modsem3)

modsem4<- lmer(Predators ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + MBC + Fungivores
               + Detritivores, data = sem_data)
modsem4
anova(modsem4, ddf='Kenward-Roger')
plot(modsem4)


sem_model1<- psem(modsem1,modsem2,modsem3,modsem4, data = sem_data)
summary(sem_model1, standardize = "scale")
plot(sem_model1)







# SUPPLEMENTARY DATA

# PERMANOVA PAIR COMPARISONS

# BACTERIA

# Pairwise comparison of multivariate, non-parametric data (community composition)
# to identify which specific groups differ from each other after a PERMANOVA.

PCoA_plot_B$group <- interaction(PCoA_plot_B$Sex, PCoA_plot_B$Location)

pairwise.adonis2(
  B.dis_b ~ group,
  data = PCoA_plot_B,
  permutations = 999
)

# FUNGI

PCoA_plot_F$group <- interaction(PCoA_plot_F$Sex, PCoA_plot_F$Location)

pairwise.adonis2(
  B.dist_F ~ group,
  data = PCoA_plot_F,
  permutations = 999
)


# ARTHROPOD

PCoA_plot$group <- interaction(PCoA_plot$Sex, PCoA_plot$Location)

pairwise.adonis2(
  Art_B_dist ~ group,
  data = PCoA_plot,
  permutations = 999
)



# PERMUTATION TEST FOR UNIQUE TAXA 


# Bacteria Permutation test


# make random number generation reproducible
set.seed(123) 

# Number of permutations
n_perm_b <- 999
# List to store permutation results
perm_results_b<- vector("list", n_perm_b)


for (i in 1:n_perm_b) {
# Randomly shuffle the Sex labels within each location
permuted_data_b <- filtered_datab %>%
group_by(Location)%>%
mutate(Sex = sample(Sex))%>%
ungroup()
# Compute summary for permuted dataset  
perm_summary_b <- get_summary_bacteria(permuted_data_b) %>%
# Label permutation iteration
mutate(perm = i)
# Store results in the list  
  perm_results_b[[i]] <- perm_summary_b
}

# Combine observed and permutation summaries into one dataset
all_summaries_b <- bind_rows(observed_summary_b, bind_rows(perm_results_b)) %>%
# Remove any rows with missing data
drop_na()

# Extract observed values by presence category

observed_vals_bf <- all_summaries_b %>%
  filter(perm == 0, presence == "female-only") %>%
  select(Location, n_taxa, presence)

observed_vals_bm <- all_summaries_b %>%
  filter(perm == 0, presence == "male-only") %>%
  select(Location, n_taxa, presence)

observed_vals_bs <- all_summaries_b %>%
  filter(perm == 0, presence == "shared") %>%
  select(Location, n_taxa,presence)

# Combine all observed values into one dataset
observed_vals_bacteria <- bind_rows(observed_vals_bf, observed_vals_bm, observed_vals_bs)

# Filter out permutation data for plotting null distribution
filtered_data_nb <- all_summaries_b %>% filter(perm != 0)

# Prepare parameters for histogram plotting

# Set consistent number of bins
bins <- 30

# Get global x limits and y-axis maximum count
x_limits <- range(filtered_data_nb$n_taxa, na.rm = TRUE)

# Estimate y max height from all data
y_max <- filtered_data_nb %>%
  group_by(Location, presence) %>%
  summarise(max_count = max(hist(n_taxa, breaks = bins, plot = FALSE)$counts), .groups = "drop") %>%
  pull(max_count) %>% max()

# Plot null distribution with observed values
combined_facet_plot_b <- ggplot(filtered_data_nb, aes(x = n_taxa)) +
  geom_histogram(bins = bins, fill = "gray", alpha = 0.7) +
  # Observed values
  geom_vline(data = observed_vals_bacteria, aes(xintercept = n_taxa), 
             color = "red", linetype = "dashed", linewidth = 1) +
  # Separate panels by location and presence category
  facet_wrap(
    Location ~ presence,
    labeller = labeller(
      Location = c(
        "Huntsville" = "PERL",
        "Turkey Creek" = "BTNP",
        "Wildflower" = "LBJWC"
      )
    )
  ) +
  coord_cartesian(xlim = x_limits, ylim = c(0, y_max)) +
  labs(x = "Number of taxa", y = "permutations") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
combined_facet_plot_b

ggsave(
  filename = "Figure_S1_null_distribution.png",
  plot = combined_facet_plot_b,
  width = 7, height = 6, units = "in",
  dpi = 600
)

ggsave(
  filename = "Figure_S1_null_distribution.tiff",
  plot = combined_facet_plot_b,
  width = 7, height = 6, units = "in",
  dpi = 600
)


# Compute 95% confidence intervals from permutations
quantile_summary_b<- filtered_data_nb %>%
  group_by(Location, presence) %>%
  summarise(
    lower_2.5 = quantile(n_taxa, 0.025, na.rm = TRUE),
    upper_97.5 = quantile(n_taxa, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Compare observed values to null confidence intervals
observed_vs_nullb <- observed_vals_bacteria %>%
  left_join(quantile_summary_b, by = c("Location", "presence")) %>%
  mutate(
    outside_CI = ifelse(n_taxa < lower_2.5 | n_taxa > upper_97.5, TRUE, FALSE)
  )



# Fungi Permutation test

# make random number generation reproducible
set.seed(123) 

# shuffle Sex labels within sites

# Number of permutations
n_perm <- 999
perm_results_f <- vector("list", n_perm)


for (i in 1:n_perm) {
  # Randomly permute sex labels within each site
  permuted_data_f <- filtered_dataf %>%
    group_by(Location)%>%
    mutate(Sex = sample(Sex))%>%
    ungroup()
  # Recalculate summary for permuted data
  perm_summary_f <- get_summary_fungi(permuted_data_f) %>%
    mutate(perm = i)
  # Store results  
  perm_results_f[[i]] <- perm_summary_f
}

# Combine observed and permuted results

all_summaries_f <- bind_rows(observed_summary_f, bind_rows(perm_results_f)) %>%
  drop_na()


# Extract observed values for each category

observed_vals_ff <- all_summaries_f %>%
  filter(perm == 0, presence == "female-only") %>%
  select(Location, n_taxa, presence)

observed_vals_fm <- all_summaries_f %>%
  filter(perm == 0, presence == "male-only") %>%
  select(Location, n_taxa, presence)

observed_vals_fs <- all_summaries_f %>%
  filter(perm == 0, presence == "shared") %>%
  select(Location, n_taxa,presence)

# Combine all observed values
observed_vals_fun <- bind_rows(observed_vals_ff, observed_vals_fm, observed_vals_fs)

# Keep only null (permuted) distributions
filtered_data_nf <- all_summaries_f %>% filter(perm != 0)

# Histogram for null distributions

# Set consistent number of bins
bins <- 30

# Get global x limits and y-axis maximum count
x_limits <- range(filtered_data_nf$n_taxa, na.rm = TRUE)

# Estimate y max height from all data
y_max <- filtered_data_nf %>%
  group_by(Location, presence) %>%
  summarise(max_count = max(hist(n_taxa, breaks = bins, plot = FALSE)$counts), .groups = "drop") %>%
  pull(max_count) %>% max()

# Plot null distributions with observed values
combined_facet_plot_f <- ggplot(filtered_data_nf, aes(x = n_taxa)) +
  # Null distribution
  geom_histogram(bins = bins, fill = "gray", alpha = 0.7) +
  # Observed value
  geom_vline(data = observed_vals_fun, aes(xintercept = n_taxa), 
             color = "red", linetype = "dashed", linewidth = 1) +
  # Facet by site and category
  facet_wrap(
    Location ~ presence,
    labeller = labeller(
      Location = c(
        "Huntsville" = "PERL",
        "Turkey Creek" = "BTNP",
        "Wildflower" = "LBJWC"
      )
    )
  ) +
  # Fix axis limits
  coord_cartesian(xlim = x_limits, ylim = c(0, y_max)) +
  # Axis labels
  labs(x = "Number of taxa", y = "permutations") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
combined_facet_plot_f

ggsave(
  filename = "Figure_f_S1_null_distribution.png",
  plot = combined_facet_plot_f,
  width = 7, height = 6, units = "in",
  dpi = 600
)

ggsave(
  filename = "Figure_f_S1_null_distribution.tiff",
  plot = combined_facet_plot_f,
  width = 7, height = 6, units = "in",
  dpi = 600
)


# Calculate 95% confidence intervals from null models
quantile_summary <- filtered_data_nf %>%
  group_by(Location, presence) %>%
  summarise(
    lower_2.5 = quantile(n_taxa, 0.025, na.rm = TRUE),
    upper_97.5 = quantile(n_taxa, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Compare observed values to null expectations
# Identify whether observed values fall outside the 95% CI
observed_vs_nullf <- observed_vals_fun %>%
  left_join(quantile_summary, by = c("Location", "presence")) %>%
  mutate(
    outside_CI = ifelse(n_taxa < lower_2.5 | n_taxa > upper_97.5, TRUE, FALSE)
  )


# Arthropod Permutation test 


# make random number generation reproducible
set.seed(123) 

# shuffle Sex labels within sites

# Number of permutations
n_perm <- 999
perm_results <- vector("list", n_perm)

for (i in 1:n_perm) {
  # Randomly permute sex labels within each site
  permuted_data <- filtered_data_art %>%
    group_by(Location) %>%
    mutate(Sex = sample(Sex)) %>%
    ungroup()
  # Recalculate summary for permuted data  
  perm_summary <- get_summary(permuted_data) %>%
    mutate(perm = i)
  # Store results   
  perm_results[[i]] <- perm_summary
}
# Combine observed and permuted results
all_summaries <- bind_rows(observed_summary, bind_rows(perm_results)) %>%
  drop_na()

# Extract observed values for each category

observed_vals_f <- all_summaries %>%
  filter(perm == 0, presence == "female-only") %>%
  select(Location, n_taxa, presence)

observed_vals_m <- all_summaries %>%
  filter(perm == 0, presence == "male-only") %>%
  select(Location, n_taxa, presence)

observed_vals_s <- all_summaries %>%
  filter(perm == 0, presence == "shared") %>%
  select(Location, n_taxa,presence)

# Combine all observed values
observed_vals <- bind_rows(observed_vals_f, observed_vals_m, observed_vals_s)

# Keep only null (permuted) distributions
filtered_data_n <- all_summaries %>% filter(perm != 0)

# Histogram for null distributions

# Set consistent number of bins
bins <- 30

# Get global x limits and y-axis maximum count
x_limits <- range(filtered_data_n$n_taxa, na.rm = TRUE)

# Estimate y max height from all data
y_max <- filtered_data_n %>%
  group_by(Location, presence) %>%
  summarise(max_count = max(hist(n_taxa, breaks = bins, plot = FALSE)$counts), .groups = "drop") %>%
  pull(max_count) %>% max()

# Plot null distributions with observed values

combined_facet_plot <- ggplot(filtered_data_n, aes(x = n_taxa)) +
  # Null distribution
  geom_histogram(bins = bins, fill = "gray", alpha = 0.7) +
  # Observed value
  geom_vline(data = observed_vals, aes(xintercept = n_taxa), 
             color = "red", linetype = "dashed", linewidth = 1) +
  # Facet by site and category
  facet_wrap(
    Location ~ presence,
    labeller = labeller(
      Location = c(
        "Huntsville" = "PERL",
        "Turkey_Creek" = "BTNP",
        "Wildflower" = "LBJWC"
      )
    )
  ) +
  # Fix axis limits
  coord_cartesian(xlim = x_limits, ylim = c(0, y_max)) +
  labs(x = "Number of taxa", y = "permutations") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))
combined_facet_plot


ggsave(
  filename = "Figure_art_S1_null_distribution.png",
  plot = combined_facet_plot,
  width = 7, height = 6, units = "in",
  dpi = 600
)

ggsave(
  filename = "Figure_art_S1_null_distribution.tiff",
  plot = combined_facet_plot,
  width = 7, height = 6, units = "in",
  dpi = 600
)

quantile_summary_ar <- filtered_data_n %>%
  group_by(Location, presence) %>%
  summarise(
    lower_2.5 = quantile(n_taxa, 0.025, na.rm = TRUE),
    upper_97.5 = quantile(n_taxa, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Checking if number of taxa happens for chance 
observed_vs_nullar <- observed_vals %>%
  left_join(quantile_summary_ar, by = c("Location", "presence")) %>%
  mutate(
    outside_CI = ifelse(n_taxa < lower_2.5 | n_taxa > upper_97.5, TRUE, FALSE)
  )



# Other candidate SEMs


# SEM #2 -  without detritivores


mod1<- lmer(MBC ~ pH + TotalN + Moisture_VWC + Sex + (1|Location), data = sem_data)
mod1
anova(mod1, ddf='Kenward-Roger')
plot(mod1)


mod2<- lmer(Fungivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + 
              MBC, data = sem_data)
mod2
anova(mod2, ddf='Kenward-Roger')
plot(mod2)

mod3<- lmer(Predators ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + Fungivores, data = sem_data)
mod3
anova(mod3, ddf='Kenward-Roger')
plot(mod3)


sem_model2<- psem(mod1,mod2,mod3, data = sem_data)
summary(sem_model2, standardize = "scale")
plot(sem_model2)




# SEM #3 - without fungivores

modA<- lmer(MBC ~ pH + TotalN + Moisture_VWC + Sex + (1|Location), data = sem_data)
modA
anova(modA, ddf='Kenward-Roger')
plot(modA)


modB<- lmer(Detritivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location), data = sem_data)
modB
anova(modB, ddf='Kenward-Roger')
plot(modB)

modC<- lmer(Predators ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + Detritivores
            + MBC, data = sem_data)
modC
anova(modC, ddf='Kenward-Roger')
plot(modC)


sem_model3<- psem(modA,modB,modC, data = sem_data)
summary(sem_model3, standardize = "scale")
plot(sem_model3)





#SEM #4 

modsem1n<- lmer(MBC ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems, data = sem_data)
modsem1n
anova(modsem1n, ddf='Kenward-Roger')
plot(modsem1n)

modsem2n<- lmer(Fungivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems + MBC, data = sem_data)
modsem2n
anova(modsem2n, ddf='Kenward-Roger')
plot(modsem2n)

modsem3n<- lmer(Detritivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems, data = sem_data)
modsem3n
anova(modsem3n, ddf='Kenward-Roger')
plot(modsem3n)

modsem4n<- lmer(Predators ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems + MBC + Fungivores
                + Detritivores, data = sem_data)
modsem4n
anova(modsem4n, ddf='Kenward-Roger')
plot(modsem4n)


sem_model4<- psem(modsem1n,modsem2n,modsem3n,modsem4n, data = sem_data)
summary(sem_model1n, standardize = "scale")
plot(sem_model1n)


# SEM #5


mod1n<- lmer(MBC ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems, data = sem_data)
mod1n
anova(mod1n, ddf='Kenward-Roger')
plot(mod1n)


mod2n<- lmer(Fungivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems +
              MBC, data = sem_data)
mod2n
anova(mod2n, ddf='Kenward-Roger')
plot(mod2n)

mod3n<- lmer(Predators ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems + Fungivores, data = sem_data)
mod3n
anova(mod3n, ddf='Kenward-Roger')
plot(mod3n)


sem_model5<- psem(mod1n,mod2n,mod3n, data = sem_data)
summary(sem_model2n, standardize = "scale")
plot(sem_model2n)




# SEM #6 

modAn<- lmer(MBC ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems, data = sem_data)
modAn
anova(modAn, ddf='Kenward-Roger')
plot(modAn)


modBn<- lmer(Detritivores ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems, data = sem_data)
modBn
anova(modBn, ddf='Kenward-Roger')
plot(modBn)

modCn<- lmer(Predators ~ pH + TotalN + Moisture_VWC + Sex + (1|Location) + N_stems + Detritivores
            + MBC, data = sem_data)
modCn
anova(modCn, ddf='Kenward-Roger')
plot(modCn)


sem_model6<- psem(modAn,modBn,modCn, data = sem_data)
summary(sem_model3n, standardize = "scale")
plot(sem_model3n)


# AIC selecting best SEM
AIC(sem_model1,sem_model2,sem_model3,sem_model4,sem_model5,sem_model6)
