# Required packages
# install.packages(c("tidyverse","biomaRt","STRINGdb","igraph","ggraph","ggplot2","scales"))
# BiocManager::install("STRINGdb")


library(tidyverse)
library(biomaRt)
library(STRINGdb)
library(igraph)
library(ggraph)
library(ggplot2)
library(scales)

# Parameters
human_file <- read.csv("/home/UpdownHumanD0vsD12.csv", header = TRUE)
mouse_file <- read.csv("/home/UpDownmMouseday0day6.csv", header = TRUE)

padj_cutoff <- 0.05
string_score_threshold <- 400
min_edges_for_plot <- 1

out_dir <- "net_outputs"
dir.create(out_dir, showWarnings = FALSE)

# Load data
human <- read_csv(human_file, col_types = cols())
mouse <- read_csv(mouse_file, col_types = cols())

# Ensure column names
human <- human_file %>% rename(gene = ID, logFC = logFC)
mouse <- mouse_file %>% rename(gene = ID, logFC = logFC)

# Optional filtering by adjusted p-value
if(!is.na(padj_cutoff) && "padj" %in% colnames(human)) human <- human %>% filter(padj <= padj_cutoff)
if(!is.na(padj_cutoff) && "padj" %in% colnames(mouse)) mouse <- mouse %>% filter(padj <= padj_cutoff)

# Split genes by regulation
human_up   <- human %>% filter(logFC > 0) %>% pull(gene) %>% unique()
human_down <- human %>% filter(logFC < 0) %>% pull(gene) %>% unique()
mouse_up   <- mouse %>% filter(logFC > 0) %>% pull(gene) %>% unique()
mouse_down <- mouse %>% filter(logFC < 0) %>% pull(gene) %>% unique()

# Ensembl ortholog mapping
hs <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

map_human_to_mouse <- function(human_genes) {
  if(length(human_genes)==0) return(tibble(human=character(), mouse=character()))
  res <- getLDS(attributes = c("hgnc_symbol"),
                filters = "hgnc_symbol",
                values = human_genes,
                mart = hs,
                attributesL = c("mgi_symbol"),
                martL = mm,
                uniqueRows = TRUE)
  colnames(res) <- c("human","mouse")
  as_tibble(res)
}

map_mouse_to_human <- function(mouse_genes) {
  if(length(mouse_genes)==0) return(tibble(mouse=character(), human=character()))
  res <- getLDS(attributes = c("mgi_symbol"),
                filters = "mgi_symbol",
                values = mouse_genes,
                mart = mm,
                attributesL = c("hgnc_symbol"),
                martL = hs,
                uniqueRows = TRUE)
  colnames(res) <- c("mouse","human")
  as_tibble(res)
}

# Load ortholog table
ortho <- read_csv("/home/humanMice.txt")

ortho <- ortho %>%
  rename(
    Mouse = `Mouse gene name`,
    Human = `Gene name`
  )

# Convert vectors to data frames
human_up   <- data.frame(gene = human_up)
human_down <- data.frame(gene = human_down)
mouse_up   <- data.frame(gene = mouse_up)
mouse_down <- data.frame(gene = mouse_down)

# Human → Mouse mapping
map_h2m_up <- human_up %>%
  inner_join(ortho, by = c("gene" = "Human")) %>%
  rename(mouse_gene = Mouse)

map_h2m_down <- human_down %>%
  inner_join(ortho, by = c("gene" = "Human")) %>%
  rename(mouse_gene = Mouse)

# Mouse → Human mapping
map_m2h_up <- mouse_up %>%
  inner_join(ortho, by = c("gene" = "Mouse")) %>%
  rename(human_gene = Human)

map_m2h_down <- mouse_down %>%
  inner_join(ortho, by = c("gene" = "Mouse")) %>%
  rename(human_gene = Human)

# Save mapping tables
write_csv(map_h2m_up, file.path(out_dir,"map_human_to_mouse_up.csv"))
write_csv(map_h2m_down, file.path(out_dir,"map_human_to_mouse_down.csv"))

# Initialize STRINGdb
string_hs <- STRINGdb$new(version="12.0", species=9606, score_threshold=string_score_threshold, input_directory="")
string_mm <- STRINGdb$new(version="12.0", species=10090, score_threshold=string_score_threshold, input_directory="")

# Extract STRING interactions for a gene set
get_string_network <- function(genes, stringdb_obj) {
  if(length(genes)==0) return(NULL)
  
  df <- tibble(gene = genes)
  mapped <- stringdb_obj$map(df, "gene", removeUnmappedRows = TRUE)
  if(nrow(mapped)==0) return(NULL)

  all_ids <- mapped$STRING_id
  edges <- stringdb_obj$get_interactions(all_ids)
  
  edges <- edges %>% filter(from %in% all_ids & to %in% all_ids)
  if(nrow(edges)==0) return(NULL)

  name_map <- mapped %>% select(STRING_id, gene) %>% distinct()

  edges <- edges %>%
    left_join(name_map, by=c("from"="STRING_id")) %>% rename(from_gene = gene) %>%
    left_join(name_map, by=c("to"="STRING_id")) %>% rename(to_gene = gene)

  edges
}

# Clean gene symbols
human_up_df <- data.frame(gene = human_up, stringsAsFactors = FALSE)
human_up_df$gene <- toupper(trimws(human_up_df$gene))

# STRING download URLs
url_aliases <- "https://stringdb-static.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
url_links   <- "https://stringdb-static.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz"

file_aliases <- "9606.protein.aliases.v12.0.txt.gz"
file_links   <- "9606.protein.links.full.v12.0.txt.gz"

# Download STRING data if missing
if(!file.exists(file_aliases)) download.file(url_aliases, file_aliases)
if(!file.exists(file_links)) download.file(url_links, file_links)

# Load STRING tables
aliases_hs <- readr::read_tsv(file_aliases) %>%
  rename(protein = `#string_protein_id`) %>%
  select(protein, alias) %>%
  mutate(alias = toupper(alias))

links_hs <- read_delim(file_links, delim = " ", col_names = TRUE, trim_ws = TRUE)

# Keep relevant columns
links_hs <- links_hs %>%
  select(protein1 = 1, protein2 = 2, combined_score = last_col())

links_hs <- links_hs %>%
  mutate(across(everything(), ~ sub("^9606\\.", "", .)))

aliases_hs <- aliases_hs %>%
  mutate(protein = sub("^9606\\.", "", protein),
         alias = toupper(alias))

# Prepare DEG tables
human_up_df   <- data.frame(gene = toupper(trimws(human_up)), stringsAsFactors = FALSE)
human_down_df <- data.frame(gene = toupper(trimws(human_down)), stringsAsFactors = FALSE)

# Map gene symbols to STRING proteins
mapped_up <- aliases_hs %>% filter(alias %in% human_up_df$gene)
mapped_down <- aliases_hs %>% filter(alias %in% human_down_df$gene)

# Extract interaction edges
edges_human_up <- links_hs %>%
  filter(protein1 %in% mapped_up$protein & protein2 %in% mapped_up$protein)

edges_human_down <- links_hs %>%
  filter(protein1 %in% mapped_down$protein & protein2 %in% mapped_down$protein)

# Build igraph networks
graph_human_up <- graph_from_data_frame(
  select(edges_human_up, protein1, protein2),
  directed = FALSE
)

graph_human_down <- graph_from_data_frame(
  select(edges_human_down, protein1, protein2),
  directed = FALSE
)

library(tidygraph)

graph_up_tbl <- as_tbl_graph(graph_human_up)
graph_down_tbl <- as_tbl_graph(graph_human_down)

# Plot STRING network
ggraph(graph_up_tbl, layout = "fr") +
  geom_edge_link(alpha = 0.3, color = "red") +
  geom_node_point(size = 3, color = "red") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  ggtitle("STRING Network – Human UP Genes")
