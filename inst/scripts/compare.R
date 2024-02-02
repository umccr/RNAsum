require(dplyr)
require(readr)
require(dracarys)
require(here)
require(glue)
require(tibble)
require(tidyr)

# gds results per subject
o <- here("nogit/comparisons")
tokpro <- Sys.getenv("ICA_ACCESS_TOKEN_PRO")
tokdev <- Sys.getenv("ICA_ACCESS_TOKEN_DEV")
s <- tibble::tribble(
  ~sbj, ~pro, ~dev,
  "SBJ04468", "20231212b3b22850/SBJ04468__L2301504/PRJ231328.results", "output/SBJ04468__L2301504/PRJ231328.results",
  "SBJ04469", "20231209eeca1fe7/SBJ04469__L2301505/MDX230556.results", "output/SBJ04469__L2301505/MDX230556.results",
  "SBJ04578", "2024011447748eb1/SBJ04578__L2400043/MDX230546.results", "output_test/SBJ04578__L2400043/MDX230546.results"
) |>
  mutate(
    pro = glue("gds://production/analysis_data/{sbj}/rnasum/{pro}"),
    dev = glue("gds://development/sehrish/rnasum_0.5.0/{sbj}/{dev}")
  ) |>
  pivot_longer(cols = c("pro", "dev"), names_to = "namespace", values_to = "gdsdir") |>
  mutate(
    outdir = file.path(o, sbj, namespace),
    token = ifelse(namespace == "pro", tokpro, tokdev)
  )

# download results
regex <- tibble::tribble(
  ~regex, ~fun,
  "^genes\\.expr\\.perc\\.html$", "foo1",
  "^genes\\.expr\\.z\\.html$", "foo2",
)
# try one first
# dracarys::dr_gds_download(gdsdir = s$gdsdir[1], outdir = s$outdir[1], token = s$token[1], page_size = 200, dryrun = TRUE, regexes = regex)
x <- s |>
  rowwise() |>
  mutate(
    dl = list(dracarys::dr_gds_download(gdsdir = gdsdir, outdir = outdir, token = token, page_size = 200, dryrun = FALSE, regexes = regex))
  ) |>
  ungroup()

# after download, open html and export csv
d <- s |>
  mutate(
    exprp = file.path(outdir, "perc.csv"),
    exprz = file.path(outdir, "z.csv")
  ) |>
  select(sbj, namespace, exprp, exprz) |>
  pivot_longer(cols = c("exprp", "exprz"), names_to = "ftype", values_to = "fpath") |>
  pivot_wider(names_from = "namespace", values_from = "fpath")


i <- 6
dev <- readr::read_csv(d$dev[i])
pro <- readr::read_csv(d$pro[i])

# now explore expression differences in reference and patient columns
# between dev and prod.
dplyr::left_join(dev, pro, by = "Gene", suffix = c(".dev", ".pro")) |>
  dplyr::mutate(
    # Ref_equal = `KIRP (TCGA).dev` == `KIRP (TCGA).pro`,
    # Ref_equal = `PANCAN (TCGA).dev` == `PANCAN (TCGA).pro`,
    # Ref_equal = `PAAD (TCGA).dev` == `PAAD (TCGA).pro`,
    Pat_equal = Patient.dev == Patient.pro,
    # Ref_diff = abs(`PANCAN (TCGA).dev` - `PANCAN (TCGA).pro`),
    # Ref_diff = abs(`PAAD (TCGA).dev` - `PAAD (TCGA).pro`),
    Pat_diff = abs(Patient.dev - Patient.pro)
  ) |>
  dplyr::select(Gene, contains("KIRP"), contains("PANCAN"), contains("Patient"), contains("diff"), everything()) |>
  dplyr::filter(Pat_diff > 0) |>
  dplyr::arrange(desc(Pat_diff)) |>
  dplyr::arrange(desc(Ref_diff))
