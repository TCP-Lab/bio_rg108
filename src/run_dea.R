
library(tidyverse)
requireNamespace("DeNumerator")
requireNamespace("DESeq2")

sanitize_names <- function(x) {
    x %>% str_replace_all(" ", "_")
}

#' Load and clean the input data
#' 
#' Paths are default for this project
load_data <- function(
        data_path = "./data/bioRG108_CountMatrix_genes_expected_count.tsv",
        metadata_path = "./data/in/Checklist_ENA-default sample checklist_1721814613698.tsv"
    ){
    raw_data <- read_tsv(data_path)
    gene_metadata <- raw_data[c("ENSEMBL", "SYMBOL", "GENENAME", "GENETYPE")]
    data <- raw_data %>% select(-all_of(c("SYMBOL", "GENENAME", "GENETYPE"))) %>%
        column_to_rownames("ENSEMBL")
    
    sample_metadata <- read_tsv(
        metadata_path,
        skip = 1, comment="#"
    )
    # We are interested in 2 conditions
    sample_metadata$media <- factor(
        ifelse(grepl("aMEM", sample_metadata$sample_alias), "amem", "mcdb131")
    )
    sample_metadata$cell_type <- as.factor(sample_metadata$cell_type %>% sanitize_names())
    sample_metadata <- sample_metadata %>% column_to_rownames("file_id")
    
    list(
        counts = data,
        metadata = list(
            gene = gene_metadata,
            sample = sample_metadata
        )
    )
}

data <- load_data()

# Check the colnames/rownames order since DeSeq2 is dumb
data$counts <- data$counts[, rownames(data$metadata$sample)]
assertthat::are_equal(colnames(data$counts), rownames(data$metadata$sample))

# We round the data when loading the matrix since Deseq doesn't like the
# RSEM partial abundances

dds <- DESeq2::DESeqDataSetFromMatrix(
    round(data$counts),
    data$metadata$sample,
    ~ cell_type + media + cell_type:media
)

# We have liftoff
dds_res <- DESeq2::DESeq(dds)

# Routine plots
dds_res %>% DESeq2::plotMA()
# Looks good to me!

# Enumerate the results
enum <- DeNumerator::denumerate(
    dds_res, alpha=0.05, fold_change=1,
    new_labels=list(
        "cell_type_Periodontal_ligament_stem_cells_vs_Dental_pulp_stem_cells" = "cell_type:ligament",
        "media_mcdb131_vs_amem" = "media:mcdb131",
        "cell_typePeriodontal_ligament_stem_cells.mediamcdb131" = "int:media_on_ligament"
    )
)

DeNumerator::plot_enumeration_frame(
    enum,
    title = "Denumeration of Biorg_108",
    exclude_all_negative = TRUE
)
