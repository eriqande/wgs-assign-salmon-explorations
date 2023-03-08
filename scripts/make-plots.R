


sample_info <- snakemake@input$sample_info
sample_order <- snakemake@input$sample_order
sim_output <- snakemake@input$sim_output
outfig <- snakemake@output$outfig

library(tidyverse)

meta <- read_tsv(sample_info)

bamlist_order <- read_tsv(sample_order, col_names = "vcf_name") %>%
  mutate(idx = 1:n())

log_likes <- read_table(
  sim_output,
  col_names = c("file", "idx", "klam_like", "sacto_like")
) %>%
  mutate(logl_for_klam = klam_like - sacto_like) %>%
  extract(
    file, 
    into = c("cov", "rep"),
    regex = "miss30/([FD0-9.]+X)/rep_([0-9]+)/mix",
    remove = FALSE
  )


LLMO <- log_likes %>%
  left_join(bamlist_order, by = "idx") %>%
  left_join(meta, by = "vcf_name")


for_plot <- LLMO %>%
  mutate(
    logl_to_correct_pop = case_when(
      ref_pop == "Klamath" ~ klam_like - sacto_like,
      ref_pop == "Sacramento" ~ sacto_like - klam_like,
      TRUE ~ NA_real_
    )
  ) %>%
  select(NMFS_DNA_ID, ref_pop, original_depth, cov, rep, logl_to_correct_pop) %>%
  mutate(
    cov_numeric = as.numeric(str_replace(cov, "X", "")),
    cov_numeric = ifelse(is.na(cov_numeric), original_depth, cov_numeric)
  ) %>%
  filter(cov_numeric <= original_depth) %>%
  arrange(ref_pop, original_depth) %>%
  mutate(samp_name = sprintf("%s: %s (%.2f)", ref_pop, NMFS_DNA_ID, original_depth)) %>%
  mutate(samp_name_f = factor(samp_name, levels = unique(samp_name))) %>%
  mutate(
    cov = recode(cov, FDX = "Full Depth")
  )

# some stuff for sorting the depths.  This is hard wired for our example
# so, would have to be modified for other depths...
cov_vals <- c("0.001X", "0.005X", "0.01X", "0.05X", "0.1X", "0.5X", "1.0X", "Full Depth")

fp2 <- for_plot %>%
  mutate(`Average \nRead Depth` = factor(cov, levels = cov_vals))

g <- ggplot(fp2, aes(y = samp_name_f, x = logl_to_correct_pop, fill = `Average \nRead Depth`)) + 
  geom_point(shape = 21) +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10), breaks = c(-5, -3, 0, 2, 5, 10, 25, 50, 100, 500, 2000, 5000, 50000)) +
  xlab("Log (base e) likelihood ratio for assignment to correct population") +
  ylab("Population: Sample_ID (original depth)") +
  scale_fill_manual(
    values = c(rainbow(n = length(levels(fp2$`Average \nRead Depth`)) - 1), "black")
  )
    
ggsave(g, filename = outfig, width = 9, height = 6)
