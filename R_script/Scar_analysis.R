### Scar analysis

### Load data and remove an sequencing error in the data
dataset <- read_tsv("scar_cell_whole.tsv") %>%
  mutate(scar = replace(scar, scar == "ACAGAGAGGGCAAAACCCCGCAAGGTGGGCTGCAG", "ACAGAGAGGGCAAACCCCGCAAGGTGGGCTGCAG"),
                              mat = replace(mat, mat == "MMMMMMIMMMMMMMMMMMMDDDM", "MMMMMMMMMMMMMMMMMMDDDMM"))


### Infer chimeric rate using cells with two scars in individual sites
sup_4 <- dataset %>%
  add_count(cell_id, site, scar, wt = umi_num) %>%
  distinct(cell_id, site, scar, .keep_all = T) %>%
  select(-umi_num, umi_num = n) %>%
  add_count(cell_id, site) %>%
  arrange(cell_id, site, desc(umi_num)) %>%
  ### retain cells with two or more scars in individual sites
  filter(n >= 2) %>% 
  mutate(rate = umi_num / t_umi) %>% group_by(cell_id, site) %>%
  # filter(umi_num<max(umi_num)) %>%
  mutate(rank = row_number()) %>% filter(rank != 1) %>%
  ungroup() %>%
  ggplot(aes(x = n, y = rate)) + geom_boxplot(aes(group = n)) +
  theme_minimal() +
  # annotate("segment", x = 1.5, xend = 14.5, y = 0.09, yend = 0.09, colour = "red") +
  scale_y_continuous(name = "Proportions of less frequent scars per site per cell",
                     labels = scales::percent, breaks = c(0, 0.09, 0.2, 0.3, 0.4, 0.5)) +
  scale_x_continuous(name = "Number of scars in the cell", breaks = c(2, 6, 10, 14))
#ggsave("Figure supplement 4.pdf", height = 5.874, width = 7.25, device=cairo_pdf)


### Utilize the 75th percentile (9%) as cut-off, the rest is considered as doublets/multiplets
Larv_scar <- dataset %>%
  add_count(cell_id, site, scar, wt = umi_num) %>% distinct(cell_id, site, scar, .keep_all = T) %>% select(-umi_num, umi_num = n) %>%
  mutate(rates = umi_num / t_umi) %>% filter(rates > 0.09) %>% add_count(site, scar, name = "cell_supp")


### Label scars
scar_lab <- Larv_scar %>% distinct(site, scar, cell_supp) %>% arrange(site, desc(cell_supp)) %>% group_by(site) %>%
  mutate(id = row_number()) %>% ungroup() %>% left_join(Larv_scar, by = c("scar", "site")) %>%
  mutate(lab = paste(site, id, sep = "#")) %>% select(site, id, lab, cell_supp = cell_supp.x, cell_id, scar, mat, umi_num, t_umi, rates)


### Remove possible doublets
### Remove scars detected in fewer than 3 cells 
### 8,737 cells and 512 scars left
doublets <- scar_lab %>% add_count(cell_id, site) %>% filter(n > 1) %>% distinct(cell_id) %>% pull()
cell_scar <- scar_lab %>% #filter(! cell_id %in% doublets) %>%
mutate(del = str_count(mat, "D"), ins = str_count(mat, "I")) %>% #select(cell_id,lab,del,ins) %>%
add_count(lab, name = "cell_num") %>% filter(cell_num >= 3) %>% arrange(cell_id)