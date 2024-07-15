require(dplyr)
require(readr)
require(here)
require(glue)
require(tibble)
require(tidyr)
require(DT)

SBJ04426_dev <- readr::read_csv("../../data/wts/RNAsum/reference_update_comp/SBJ04426/dev/genes.expr.perc.csv")
SBJ04426_pro <- readr::read_csv("../../data/wts/RNAsum/reference_update_comp/SBJ04426/pro/genes.expr.perc.csv")
SBJ04187_dev <- readr::read_csv("../../data/wts/RNAsum/reference_update_comp/SBJ04187/dev/genes.expr.perc.csv")
SBJ04187_pro <- readr::read_csv("../../data/wts/RNAsum/reference_update_comp/SBJ04187/pro/genes.expr.perc.csv")
cancer_genes <- readr::read_tsv("../../research/data/cancer_gene_list/somatic_panel-v24.03.0.tsv")

# now explore expression differences in reference and patient columns
# between dev and prod.
SBJ04426_df <- dplyr::left_join(SBJ04426_dev, SBJ04426_pro, by = "Gene", suffix = c(".dev", ".pro")) |>
  dplyr::mutate(
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

SBJ04187_df <- dplyr::left_join(SBJ04187_dev, SBJ04187_pro, by = "Gene", suffix = c(".dev", ".pro")) |>
  dplyr::mutate(
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
hist(SBJ04426_df[[1]]$data$Ref_diff, breaks = 100)


