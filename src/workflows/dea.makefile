

## --- Decompress sources
./data/%: ./data/in/%.gz
	gunzip -cfv $< > $@

# If we have the input data as-is in the data/in folder, but we need it in
# data/ we can just copy it.
./data/%: ./data/in/%
	cp $< $@

ALL += ./data/out/done.flag
./data/out/done.flag : ./data/bioRG108_CountMatrix_genes_expected_count.tsv
	Rscript --vanilla ./src/run_dea.R
	touch $@

.ALL = ALL
