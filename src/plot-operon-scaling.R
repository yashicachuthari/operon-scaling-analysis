library(tidyverse)

input_file <- "processed_data/operon_scaling_results.csv"
output_dir <- "results/figures"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read_csv(input_file, show_col_types = FALSE)

plot_data <- data %>%
  mutate(
    operon_count = replace_na(operon_count, 0),
    replicon_type = str_to_lower(replicon_type)
  ) %>%
  filter(
    !is.na(length_bp),
    length_bp > 0,
    replicon_type %in% c("chromosome", "plasmid")
  )

print(plot_data %>% count(replicon_type))

p <- ggplot(plot_data,
            aes(x = length_bp,
                y = operon_count,
                color = replicon_type)) +
  geom_point(alpha = 0.25, size = 0.8) +
  scale_color_manual(
    values = c(
      chromosome = "#2563EB",
      plasmid = "#E76F51"
    )
  ) +
  labs(
    title = "Operon Scaling Across Bacterial Replicons",
    x = "Replicon length (bp)",
    y = "Number of operons",
    color = "Replicon type"
  ) +
  theme_classic(base_size = 13)

ggsave(
  "results/figures/operon_scaling_R.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  "results/figures/operon_scaling_R.pdf",
  plot = p,
  width = 9,
  height = 6
)

print("Plot saved successfully!")
