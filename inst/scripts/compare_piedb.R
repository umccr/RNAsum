require(dplyr)
require(readr)
require(here)
require(glue)
require(tibble)
require(tidyr)

dev <- readr::read_csv("~/UMCCR/data/wts/RNAsum/reference_update_comp/SBJ04426/dev/genes.expr.perc.csv")
pro <- readr::read_csv("~/UMCCR/data/wts/RNAsum/reference_update_comp/SBJ04426/pro/genes.expr.perc.csv")
dev <- readr::read_csv("~/UMCCR/data/wts/RNAsum/reference_update_comp/SBJ04187/dev/genes.expr.perc.csv")
pro <- readr::read_csv("~/UMCCR/data/wts/RNAsum/reference_update_comp/SBJ04187/pro/genes.expr.perc.csv")
cancer_genes <- readr::read_tsv("~/UMCCR/research/data/cancer_gene_list/somatic_panel-v24.03.0.tsv")

# now explore expression differences in reference and patient columns
# between dev and prod.
df <- dplyr::left_join(dev, pro, by = "Gene", suffix = c(".dev", ".pro")) |>
  dplyr::mutate(
    # Ref_equal = `KIRP (TCGA).dev` == `KIRP (TCGA).pro`,
    # Ref_equal = `PANCAN (TCGA).dev` == `PANCAN (TCGA).pro`,
    Ref_equal = `BRCA (TCGA).dev` == `BRCA (TCGA).pro`,
    Pat_equal = Patient.dev == Patient.pro,
    # Ref_diff = abs(`PANCAN (TCGA).dev` - `PANCAN (TCGA).pro`),
    Ref_diff = abs(`BRCA (TCGA).dev` - `BRCA (TCGA).pro`),
    Pat_diff = abs(Patient.dev - Patient.pro)
  ) |>
  dplyr::select(Gene, contains("BRCA"), contains("PANCAN"), Ref_diff, contains("Patient"), Pat_diff, contains("equal")) |>
  dplyr::filter(Pat_diff > 0) |>
  dplyr::filter(Gene %in% cancer_genes$ensembl_gene_symbol) |>
  dplyr::arrange(desc(Pat_diff)) |>
  dplyr::arrange(desc(Ref_diff)) |>
  datatable()

# plot Ref_diff values
hist(df$Ref_diff, breaks = 100)
