requireNamespace("ggplot2")
requireNamespace("assertthat")

to_enumeration_vector <- function(
        data, alpha = 0.01, fold_change = 1,
        zero_lab = "zero", pos_lab = "positive", neg_lab = "negative",
        colname = "class"
) {
    # This makes just a single dframe with one col: the enumeration
    new_data <- data[, "baseMean", drop=FALSE]
    # This is just to not lose the rownames
    
    res <- factor(rep(zero_lab, nrow(data)), levels = c(zero_lab, pos_lab, neg_lab))
    
    res[data$log2FoldChange > fold_change & data$padj < alpha] <- pos_lab
    res[data$log2FoldChange < - fold_change & data$padj < alpha] <- neg_lab
    
    new_data[[colname]] <- res
    new_data$baseMean <- NULL
    
    return(new_data)
}

apply_to_results <- function(deseq_object, results, fun, ..., pass_names = FALSE) {
    res <- list()
    for (result in results) {
        print(result)
        res[[result]] <- fun(DESeq2::results(deseq_object, name = result), ...)
    }
    
    res
}

plot_enumeration_frame <- function(
    enum_data,
    zero_lab = "zero", pos_lab = "positive", neg_lab = "negative"
) {
    frequencies <- enum_data %>% as.data.frame() %>% table()
    frequencies %>% ggplot(aes(x = Freq, y = variable)) +
        geom_bar(stat = "identity") +
        ylab("") + xlab("Gene Count") +
        ggtitle(title) +
        theme_bw() +
        geom_text(stat="identity", aes(label = Freq), hjust = -0.1) +
        coord_cartesian(xlim = c(0, max(x$Freq) * 1.10))
}


denumerate <- function(computed_deseq_object, results = "nointercept") {
    # Parse the "results" parameter
    possible_results <- DESeq2::resultsNames(computed_deseq_object)
    if (results == "nointercept") {
        assertthat::assert_that(possible_results[1] == "Intercept")
        results <- possible_results[-1]
    } else if (results == "all") {
        results <- possible_results
    } else {
        assertthat::assert_that(
            all(results %in% possible_results)
        )
    }
    
    # Execute the enumerations and rename the internal frames, so we can merge
    enumerations <- apply_to_results(dds_res, results, to_enumeration_vector, pass_names = TRUE)
    # I merge just with cbind(), so I want to be extra sure the rows are all
    # in the same order
    test_row_order <- enumerations[[1]] %>% rownames()
    for (name in names(enumerations)) {
        colnames(enumerations[[name]]) <- name
        assertthat::are_equal(rownames(enumerations[[name]]), test_row_order)
    }
    enum_frame <- Reduce(cbind, enumerations)
    
    plot_enumeration_frame(enum_frame)
}

denumerate(dds_res)
