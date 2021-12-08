get_sumstats_GWAScatalog <- function(study) {
  # Avoid useless queries to GWAScatalog
  stopifnot(
    stringr::str_detect(
      study,
      stringr::regex("\\AGCST[0-9]{6}\\z",
                     multiline = TRUE
      )
    )
  )

  message("Querying GWAScatalog...")
  # Get everything first
  query <-
    gwasrapidd::get_associations(
      study_id = study,
      verbose = TRUE
    )
  message("Done.")

  data_rsid <-
    dplyr::full_join(
      query@associations,
      query@risk_alleles,
      by = "association_id"
    ) %>%
    # TODO add chr and pos (they're in GRCh38)
    dplyr::select(
      variant_id,
      risk_allele,
      risk_frequency,
      beta_number,
      beta_direction,
      pvalue
    ) %>%
    dplyr::mutate(
      beta = dplyr::if_else(
        beta_direction == "decrease", # if beta decreases effect
        -1 * beta_number, # then make it negative
        beta_number # else keep it positive
      ),
      other_allele = dplyr::case_when(
        risk_allele == "A" ~ "T",
        risk_allele == "T" ~ "A",
        risk_allele == "G" ~ "C",
        risk_allele == "C" ~ "G"
      )
    ) %>%
    dplyr::select(-c("beta_direction", "beta_number"))

  message("Getting GRCh37/hg19 data from myVariant.info...")
  data_rsid %$%
    # Submit query to myvariant
    myvariant::queryVariants(
      qterms = variant_id,
      scopes = "dbsnp.rsid",
      fields = c("dbsnp.chrom", "vcf.position")
    ) %>%
    # Change obj type from DataFrame to tibble
    tibble::as_tibble() %>%
    # Select relevant columns and rename them
    dplyr::select(query, dbsnp.chrom, vcf.position) %>%
    dplyr::rename(
      chr = dbsnp.chrom,
      position = vcf.position
    ) %>%
    dplyr::full_join(
      data_rsid,
      by = c("query" = "variant_id")
    ) %>%
    dplyr::rename(rsid = query)
}
