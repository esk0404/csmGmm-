library(dplyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(data.table)
library(xtable)
library(devtools)
library(here)
library(tibble)
library(csmGmm)


# MANHATTAN PLOT (Fig 5)-----------------------------------------------------------------------------------

# for colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# plot manhattan function
plotManhattan <- function(plotRes, chrCounts, colValues, shapeValues, ylimits, legName) {
  # arrange data by chromosome
  plotRes <- plotRes %>% arrange(Chr)
  uniqueChrs <- sort(unique(plotRes$Chr))
  
  chrOffsets <- cumsum(chrCounts)
  truePos <- rep(NA, nrow(plotRes))
  counter <- 1
  for (chr_it in 1:length(uniqueChrs)) {
    tempChr <- uniqueChrs[chr_it]
    tempDat <- plotRes %>% filter(Chr == tempChr)
    # offset = sum of all previous chromosomes' max positions
    offsetVal <- ifelse(tempChr == 1, 0, chrOffsets[tempChr - 1])
    truePos[counter:(counter + nrow(tempDat) - 1)] <- offsetVal + tempDat$BP
    counter <- counter + nrow(tempDat)
  }
  
  # calculate x-axis positions
  # xBreaks <- sapply(1:22, function(i) {
  #start <- ifelse(i == 1, 0, chrOffsets[i - 1])
  #end <- chrOffsets[i]
  #(start + end) / 2
  #  })
  
  # xBreaksLabs <- ifelse((1:22) %% 2 == 0, "", as.character(1:22))
  
  # x-axis ticks at the end of each chromosome
  xBreaks <- chrOffsets[1:22]  
  xBreaksLabs <- ifelse((1:22) %% 2 == 0, "", as.character(1:22))
  
  plotDat <- plotRes %>% mutate(truePos = truePos)
  
  returnPlot <- ggplot(plotDat, aes(x=truePos, y=-log10(newLfdr), color=as.factor(cat), shape=as.factor(cat))) +
    geom_point() +
    xlab("Chromosome") + ylab("-log10(lfdr)") +
    #scale_color_manual(name="Group", values=c(gg_color_hue(3))) +
    scale_color_manual(name=legName, values=colValues) +
    scale_shape_manual(name=legName, values=shapeValues) +
    scale_x_continuous(name="Chr", breaks=xBreaks, labels=xBreaksLabs) +
    ylim(ylimits) +
    theme_cowplot() +
    theme(axis.text=element_text(size=18), axis.title=element_text(size=18),
          legend.title=element_text(size=18), legend.text=element_text(size=18)) +
    guides(colour = guide_legend(override.aes = list(size=4)))
  
  
  return(returnPlot)
}


#manDataTwo <- pleio_analysis %>% 
 # select(variant_id_iCOGS, chromosome, base_pair_location, lfdrResults) %>% 
  #mutate(cat = "iCOGS, VA MVP") %>%   # Category label for legend
  #mutate(Chr = as.numeric(chromosome), BP = as.numeric(base_pair_location)) %>%
  #arrange(lfdrResults) %>%
  #distinct(variant_id_iCOGS, .keep_all = TRUE) %>%
  #rename(newLfdr = lfdrResults)

manDataTwo <- pleio_analysis[, c("variant_id_iCOGS", "chromosome", "base_pair_location", "lfdrResults")]
manDataTwo$cat <- "iCOGS, VA MVP"
manDataTwo$Chr <- as.numeric(manDataTwo$chromosome)
manDataTwo$BP <- as.numeric(manDataTwo$base_pair_location)
manDataTwo <- manDataTwo[order(manDataTwo$lfdrResults), ]
manDataTwo <- manDataTwo[!duplicated(manDataTwo$variant_id_iCOGS), ]
manDataTwo$newLfdr <- manDataTwo$lfdrResults


chrCounts <- manDataTwo %>%
  group_by(Chr) %>%
  summarise(chrLength = max(BP, na.rm = TRUE)) %>%
  arrange(Chr) %>%
  pull(chrLength)

chrOffsets <- cumsum(chrCounts)

manPlotTwo <- plotManhattan(plotRes = manDataTwo, chrCounts,
                            colValues = gg_color_hue(1), shapeValues=c(16), ylimits=c(0, 11), legName="")


manPlotTwo

# BIVARIATE Z SCORE (Fig 6)--------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)

biv_data <- merged_data_1 %>%
  select(Z_iCOGS, Z_biobank_aligned) 

# Basic bivariate Z plot
z_plot <- ggplot(biv_data, aes(x = Z_iCOGS, y = Z_biobank_aligned)) +
  geom_point(color = "red", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Z-score (iCOGS)") +
  ylab("Z-score (VA MVP)") +
  theme_bw(base_size = 16)

z_plot

# MIAMI PLOT (Fig 7)------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(cowplot)


iCOGS_miami <- iCOGS_new %>%
  transmute(
    Chr = as.numeric(chromosome),
    BP  = as.numeric(base_pair_location),
    logp = -log10(p_value)
  ) %>%
filter(!is.na(Chr), Chr >= 1, Chr <= 22) 

rm(iCOGS_new)

VA_MVP_miami <- VA_MVP_new %>%
  transmute(
    Chr = as.numeric(chromosome),
    BP  = as.numeric(base_pair_location),
    logp = -log10(p_value)
  ) %>%
  filter(!is.na(Chr), Chr >= 1, Chr <= 22)

rm(VA_MVP_new)

chr_table <- bind_rows(
  iCOGS_miami %>% select(Chr, BP),
  VA_MVP_miami %>% select(Chr, BP)
) %>%
  group_by(Chr) %>%
  summarise(chrLength = max(BP, na.rm = TRUE)) %>%
  arrange(Chr)

chrOffsets <- cumsum(chr_table$chrLength)
names(chrOffsets) <- chr_table$Chr

add_true_pos <- function(df, chrOffsets) {
  df %>%
    arrange(Chr, BP) %>%
    mutate(
      truePos = BP + ifelse(Chr == 1, 0, chrOffsets[as.character(Chr - 1)])
    )
}

iCOGS_miami <- add_true_pos(iCOGS_miami, chrOffsets)
VA_MVP_miami <- add_true_pos(VA_MVP_miami, chrOffsets)

xBreaks <- chrOffsets[as.character(1:22)]
xBreaksLabs <- ifelse((1:22) %% 2 == 0, "", as.character(1:22))

miami_plot <- ggplot() +
  # Top: gcst1 (iCOGS)
  geom_point(
    data = iCOGS_miami,
    aes(x = truePos, y = logp),
    color = ifelse(iCOGS_miami$Chr %% 2 == 0, "grey70", "blue"),
    alpha = 0.7,
    size = 0.8
  ) +
  # Bottom: gcst3 (MVP)
  geom_point(
    data = VA_MVP_miami,
    aes(x = truePos, y = -logp),
    color = ifelse(VA_MVP_miami$Chr %% 2 == 0, "grey70", "red"),
    alpha = 0.7,
    size = 0.8
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = xBreaks,
    labels = xBreaksLabs
  ) +
  scale_y_continuous(
    name = "-log10(p-value)",
    limits = c(-30, 80),
    breaks = seq(-30, 80, 10)
  ) +
  theme_cowplot() +
  theme(
    axis.text  = element_text(size = 16),
    axis.title = element_text(size = 16)
  )

ggsave(
  filename = "new_miami_plot.png",
  plot = miami_plot,
  width = 14,
  height = 6,
  dpi = 300
)


miami_plot


# HYPOTHESIS PLOT------------------------------------------------------------------

dev.off()
library(dplyr)
library(tidyr)
library(ggplot2)

# Example 9 hypotheses
hypotheses <- expand.grid(UKBB = c(-1,0,1),
                          MVP   = c(-1,0,1)) %>%
  mutate(HypothesisID = row_number())

# Long format
hypotheses_long <- hypotheses %>%
  pivot_longer(cols = c(UKBB, MVP),
               names_to = "Study",
               values_to = "Effect") %>%
  mutate(Sign = case_when(
    Effect == 1 ~ "+",
    Effect == -1 ~ "-",
    TRUE ~ " "
  ))

# Identify blank and opposite rows
blank_rows <- hypotheses_long %>% filter(Sign == " ") %>% pull(HypothesisID)
opposite_rows <- hypotheses %>% filter((UKBB == 1 & MVP == -1) | (UKBB == -1 & MVP == 1)) %>% pull(HypothesisID)

# Remove ID row at bottom
combined <- hypotheses_long
combined$Study <- factor(combined$Study, levels = c("UKBB","MVP"))

# Create a separate data frame for top HypothesisID row (no color or border)
top_row <- data.frame(
  Study = "Configuration ID",
  HypothesisID = 1:9,
  Sign = as.character(1:9)
)

# Plot
ggplot() +
  geom_tile(data = combined,
            aes(y = Study,
                x = factor(HypothesisID, levels = 1:9),
                fill = case_when(
                  HypothesisID %in% blank_rows ~ "blank",
                  HypothesisID %in% opposite_rows ~ "opposite",
                  TRUE ~ "normal"
                )),
            color = "black") +
  
  # Text signs 
  geom_text(data = combined,
            aes(y = Study,
                x = factor(HypothesisID, levels = 1:9),
                label = Sign),
            size = 5, vjust = 0.5, hjust = 0.5) +
  
  # Top row with Hypothesis numbers (no tiles/borders)
  geom_text(data = top_row,
            aes(y = Study, x = factor(HypothesisID, levels = 1:9), label = Sign),
            size = 5, vjust = 0.5, hjust = 0.5) +

  
  # Fill colors for tiles
  scale_fill_manual(values = c("blank" = "plum", "opposite" = "pink", "normal" = "white"),
                    guide = "none") +
  
  # Layout
  scale_x_discrete(expand = c(0,0)) +
  coord_fixed(ratio = 1.2, clip = "off") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

#-------------------------------------------------------------------------------
