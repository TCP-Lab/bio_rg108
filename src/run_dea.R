
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
    
    # Remove the zero-counts. This has no effect on DESeq (see https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#Pre-filtering)
    # but it makes the object smaller and the computation (slightly) faster
    data <- data[rowSums(data) > 0,]
    
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
symbol_data <- read_csv("data/symbols.csv") %>% rename(
    ensg = "Gene stable ID",
    description = "Gene description",
    symbol = "Gene name",
    gene_type = "Gene type"
)

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

p <- DeNumerator::plot_enumeration_frame(
    enum,
    title = "Denumeration of Bio_rg108",
    exclude_all_negative = TRUE,
    category_renames = list(
        "cell_type.ligament" = "Ligament vs Pulp",
        "media.mcdb131" = "MCDB vs aMEM",
        "int.media_on_ligament" = "MCDB + Ligament"
    ),
    labels_y_nudge = 0.18
)

make_saver <- function(enum, symbol_data) {
    function(
        cell="zero", media="zero", inter="zero",
        file_path = file.path(
            "data", "out",
            paste0("enum_", "cell_", cell, "_media_", media, "_int_", inter, ".csv")
        )
    ) {
        selection <- row.names(enum)[
            enum$`cell_type:ligament` == cell & enum$`media:mcdb131` == media &
                enum$`int:media_on_ligament`== inter
        ]
        frame <- data.frame(ensg = selection) %>% merge(symbol_data, by="ensg", all.x = TRUE, all.y = FALSE)
        
        write_csv(frame, file_path)
        cat(paste0("Written ", file_path, "\n"))
    }
}

save_enum <- make_saver(enum, symbol_data)

pdf("./data/out/denumeration_plot.pdf", width = 8, height = 8)
print(p)
dev.off()

save_enum(inter = "negative")
save_enum(inter = "positive")

save_enum(cell = "negative")
save_enum(cell = "positive")

save_enum(media = "negative")
save_enum(media = "positive")

save_enum(media = "negative", inter = "positive")
save_enum(media = "positive", inter = "negative")
