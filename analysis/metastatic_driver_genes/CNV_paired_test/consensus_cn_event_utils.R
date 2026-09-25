create_consensus_dataframe <- function(cn_event_calls) {
  events <- cn_event_calls[!duplicated(cn_event_calls[c("seqnames", "start", "end")]), c("seqnames", "start", "end", "width")]
  events$is_relative_amplification_mean_clonal <- FALSE
  events$is_absolute_amplification_clonal <- FALSE
  events$is_relative_amplification_mean_subclonal <- FALSE
  events$is_absolute_amplification_subclonal <- FALSE
  events$is_relative_amplification_mean_subclonal_prop <- NA
  events$is_absolute_amplification_subclonal_prop <- NA
  events$is_relative_gain_mean_clonal <- FALSE
  events$is_absolute_gain_clonal <- FALSE
  events$is_relative_gain_mean_subclonal <- FALSE
  events$is_absolute_gain_subclonal <- FALSE
  events$is_relative_gain_mean_subclonal_prop <- NA
  events$is_absolute_gain_subclonal_prop <- NA
  events$is_relative_loss_mean_clonal <- FALSE
  events$is_absolute_loss_clonal <- FALSE
  events$is_relative_loss_mean_subclonal <- FALSE
  events$is_absolute_loss_subclonal <- FALSE
  events$is_relative_loss_mean_subclonal_prop <- NA
  events$is_absolute_loss_subclonal_prop <- NA
  events$is_relative_amplification_ttest_clonal <- FALSE
  events$is_relative_amplification_ttest_subclonal <- FALSE
  events$is_relative_amplification_ttest_subclonal_prop <- NA
  events$is_relative_gain_ttest_clonal <- FALSE
  events$is_relative_gain_ttest_subclonal <- FALSE
  events$is_relative_gain_ttest_subclonal_prop <- NA
  events$is_relative_loss_ttest_clonal <- FALSE
  events$is_relative_loss_ttest_subclonal <- FALSE
  events$is_relative_loss_ttest_subclonal_prop <- NA
  events$is_LOH_clonal <- FALSE
  events$is_LOH_subclonal <- FALSE
  events$is_LOH_subclonal_prop <- NA
  events$is_homozygous_deletion_clonal <- FALSE
  events$is_homozygous_deletion_subclonal <- FALSE
  events$is_homozygous_deletion_subclonal_prop <- NA
  events$is_ai_clonal <- FALSE
  events$is_ai_subclonal <- FALSE
  events$is_ai_subclonal_prop <- NA
  events$is_msai <- FALSE
  events
}
calc_consensus <- function(event_type, sample_cn_events, events_df, num_samples) {
  frequency_table <- sample_cn_events[, c(event_type, "seqnames", "start")] %>% dplyr::count(sample_cn_events[[event_type]], sample_cn_events$seqnames, sample_cn_events$start)
  for (j in seq_len(nrow(frequency_table))) {
    if ((frequency_table[j, "event"] == TRUE) & (frequency_table[j, "frequency"] == num_samples) & (!is.na(frequency_table[j, "event"]))) {
      clonalitycall <- paste0(event_type, "_clonal")
      events_df[which((as.character(events_df$seqnames) == as.character(frequency_table$chrom[j])) &
        (events_df$start == frequency_table$start[j])), clonalitycall] <- TRUE
    } else if ((frequency_table[j, "event"] == TRUE) & (frequency_table[j, "frequency"] != num_samples) & (!is.na(frequency_table[j, "event"]))) {
      clonalitycall <- paste0(event_type, "_subclonal")
      events_df[which((as.character(events_df$seqnames) == as.character(frequency_table$chrom[j])) &
        (events_df$start == frequency_table$start[j])), clonalitycall] <- TRUE
      subclonalityprop <- paste0(event_type, "_subclonal_prop")
      events_df[which((as.character(events_df$seqnames) == as.character(frequency_table$chrom[j])) &
        (events_df$start == frequency_table$start[j])), subclonalityprop] <- ((frequency_table[j, "frequency"]) / num_samples)
    }
  }
  events_df
}
calc_consensus_parallel <- function(modifiable_events_df, phased_segs_w_cn_events, cn_event_columns_df) {
  cn_event_system_vec <- colnames(cn_event_columns_df)
  msai_events_df <- phased_segs_w_cn_events %>%
    dplyr::mutate(ai_and_maj_allele_is_A = paste0(.[["is_ai"]], "_", .[["mirrored_vs_ref"]])) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarise(
      ai_all_ref = ifelse(!any(ai_and_maj_allele_is_A %in% "TRUE_TRUE") && any(ai_and_maj_allele_is_A %in% "TRUE_FALSE"), TRUE, FALSE),
      ai_all_alt = ifelse(any(ai_and_maj_allele_is_A %in% "TRUE_TRUE") && !any(ai_and_maj_allele_is_A %in% "TRUE_FALSE"), TRUE, FALSE),
      is_msai = ifelse(any(ai_and_maj_allele_is_A %in% "TRUE_TRUE") && any(ai_and_maj_allele_is_A %in% "TRUE_FALSE"), TRUE, FALSE)
    )
  modifiable_events_df$is_msai <- msai_events_df$is_msai
  modifiable_events_df$ai_all_ref <- msai_events_df$ai_all_ref
  modifiable_events_df$ai_all_alt <- msai_events_df$ai_all_alt
  parallel_LOH_events_df <- phased_segs_w_cn_events %>%
    dplyr::mutate(LOH_maj_cn = paste0(.[["is_LOH"]], "_", .[["mirrored_vs_ref"]])) %>%
    dplyr::group_by(seqnames, start, end) %>%
    dplyr::summarise(
      LOH_all_ref = ifelse(!any(LOH_maj_cn %in% "TRUE_TRUE") && any(LOH_maj_cn %in% "TRUE_FALSE"), TRUE, FALSE),
      LOH_all_alt = ifelse(any(LOH_maj_cn %in% "TRUE_TRUE") && !any(LOH_maj_cn %in% "TRUE_FALSE"), TRUE, FALSE),
      is_LOH_parallel = ifelse(any(LOH_maj_cn %in% "TRUE_TRUE") && any(LOH_maj_cn %in% "TRUE_FALSE"), TRUE, FALSE)
    )
  modifiable_events_df$is_LOH_parallel <- parallel_LOH_events_df$is_LOH_parallel
  modifiable_events_df$LOH_all_ref <- parallel_LOH_events_df$LOH_all_ref
  modifiable_events_df$LOH_all_alt <- parallel_LOH_events_df$LOH_all_alt
  for (cn_event_system in cn_event_system_vec) {
    cn_events_vec <- cn_event_columns_df[, cn_event_system]
    parallel_summary_events <- phased_segs_w_cn_events %>%
      dplyr::mutate(amp_maj_cn = paste0(.[["is_ai"]], "_", .[[cn_events_vec["amplification"]]], "_", .[["mirrored_vs_ref"]])) %>%
      dplyr::mutate(gain_maj_cn = paste0(.[["is_ai"]], "_", .[[cn_events_vec["gain"]]], "_", .[["mirrored_vs_ref"]])) %>%
      dplyr::mutate(loss_maj_cn = paste0(.[["is_ai"]], "_", .[[cn_events_vec["loss"]]], "_", .[["mirrored_vs_ref"]])) %>%
      dplyr::mutate(strict_gain_maj_cn = paste0(.[["is_ai"]], "_", .[[cn_events_vec["gain"]]], "_", .[[cn_events_vec["amplification"]]], "_", .[["mirrored_vs_ref"]])) %>%
      dplyr::group_by(seqnames, start, end) %>%
      dplyr::summarise(
        amp_all_ai_ref = ifelse(!any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        gain_all_ai_ref = ifelse(!any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        loss_all_ai_ref = ifelse(!any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        strict_gain_all_ai_ref = ifelse(!any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && !any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        amp_all_ai_alt = ifelse(any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && !any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        gain_all_ai_alt = ifelse(any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && !any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        loss_all_ai_alt = ifelse(any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && !any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        strict_gain_all_ai_alt = ifelse(any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && !any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && !any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        amp_ai_ref_bal = ifelse(any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && !any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        gain_ai_ref_bal = ifelse(any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && !any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        loss_ai_ref_bal = ifelse(any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && !any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        strict_gain_ai_ref_bal = ifelse(any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && !any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        amp_ai_alt_bal = ifelse(!any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        gain_ai_alt_bal = ifelse(!any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        loss_ai_alt_bal = ifelse(!any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        strict_gain_ai_alt_bal = ifelse(!any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        amp_all_bal = ifelse(!any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && !any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        gain_all_bal = ifelse(!any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && !any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        loss_all_bal = ifelse(!any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && !any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        strict_gain_all_bal = ifelse(!any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && !any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        parallel_and_bal_amp = ifelse(any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        parallel_and_bal_gain = ifelse(any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        parallel_and_bal_loss = ifelse(any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        parallel_and_bal_strict_gain = ifelse(any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        parallel_no_bal_amp = ifelse(any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && any("TRUE_TRUE_FALSE" %in% amp_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% amp_maj_cn), TRUE, FALSE),
        parallel_no_bal_gain = ifelse(any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && any("TRUE_TRUE_FALSE" %in% gain_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% gain_maj_cn), TRUE, FALSE),
        parallel_no_bal_loss = ifelse(any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && any("TRUE_TRUE_FALSE" %in% loss_maj_cn) && !any(c("FALSE_TRUE_FALSE", "NA_TRUE_FALSE") %in% loss_maj_cn), TRUE, FALSE),
        parallel_no_bal_strict_gain = ifelse(any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn) && !any("FALSE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE),
        parallel_amp = ifelse(any("TRUE_TRUE_TRUE" %in% amp_maj_cn) && any("TRUE_TRUE_FALSE" %in% amp_maj_cn), TRUE, FALSE),
        parallel_gain = ifelse(any("TRUE_TRUE_TRUE" %in% gain_maj_cn) && any("TRUE_TRUE_FALSE" %in% gain_maj_cn), TRUE, FALSE),
        parallel_loss = ifelse(any("TRUE_TRUE_TRUE" %in% loss_maj_cn) && any("TRUE_TRUE_FALSE" %in% loss_maj_cn), TRUE, FALSE),
        parallel_strict_gain = ifelse(any("TRUE_TRUE_FALSE_TRUE" %in% strict_gain_maj_cn) && any("TRUE_TRUE_FALSE_FALSE" %in% strict_gain_maj_cn), TRUE, FALSE)
      ) %>%
      dplyr::ungroup()
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_all_ai_ref")] <- parallel_summary_events$amp_all_ai_ref
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_all_ai_ref")] <- parallel_summary_events$gain_all_ai_ref
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_all_ai_ref")] <- parallel_summary_events$loss_all_ai_ref
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_all_ai_ref")] <- parallel_summary_events$strict_gain_all_ai_ref
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_all_ai_alt")] <- parallel_summary_events$amp_all_ai_alt
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_all_ai_alt")] <- parallel_summary_events$gain_all_ai_alt
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_all_ai_alt")] <- parallel_summary_events$loss_all_ai_alt
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_all_ai_alt")] <- parallel_summary_events$strict_gain_all_ai_alt
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_ai_ref_bal")] <- parallel_summary_events$amp_ai_ref_bal
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_ai_ref_bal")] <- parallel_summary_events$gain_ai_ref_bal
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_ai_ref_bal")] <- parallel_summary_events$loss_ai_ref_bal
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_ai_ref_bal")] <- parallel_summary_events$strict_gain_ai_ref_bal
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_ai_alt_bal")] <- parallel_summary_events$amp_ai_alt_bal
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_ai_alt_bal")] <- parallel_summary_events$gain_ai_alt_bal
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_ai_alt_bal")] <- parallel_summary_events$loss_ai_alt_bal
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "ai_alt_bal")] <- parallel_summary_events$strict_gain_ai_alt_bal
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_all_bal")] <- parallel_summary_events$amp_all_bal
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_all_bal")] <- parallel_summary_events$gain_all_bal
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_all_bal")] <- parallel_summary_events$loss_all_bal
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_all_bal")] <- parallel_summary_events$strict_gain_all_bal
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_parallel_and_bal")] <- parallel_summary_events$parallel_and_bal_amp
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_parallel_and_bal")] <- parallel_summary_events$parallel_and_bal_gain
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_parallel_and_bal")] <- parallel_summary_events$parallel_and_bal_loss
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_parallel_and_bal")] <- parallel_summary_events$parallel_and_bal_strict_gain
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_parallel_no_bal")] <- parallel_summary_events$parallel_no_bal_amp
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_parallel_no_bal")] <- parallel_summary_events$parallel_no_bal_gain
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_parallel_no_bal")] <- parallel_summary_events$parallel_no_bal_loss
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_parallel_no_bal")] <- parallel_summary_events$parallel_no_bal_strict_gain
    modifiable_events_df[, paste0(cn_events_vec["amplification"], "_parallel")] <- parallel_summary_events$parallel_amp
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_parallel")] <- parallel_summary_events$parallel_gain
    modifiable_events_df[, paste0(cn_events_vec["loss"], "_parallel")] <- parallel_summary_events$parallel_loss
    modifiable_events_df[, paste0(cn_events_vec["gain"], "_strict", "_parallel")] <- parallel_summary_events$parallel_strict_gain
  }
  modifiable_events_df
}
initialise_msai_columns <- function(events) {
  events$is_relative_amplification_mean_clonal_nomsai <- FALSE
  events$is_absolute_amplification_clonal_nomsai <- FALSE
  events$is_relative_amplification_mean_subclonal_nomsai <- FALSE
  events$is_relative_amplification_mean_subclonal_msai <- FALSE
  events$is_absolute_amplification_subclonal_nomsai <- FALSE
  events$is_absolute_amplification_subclonal_msai <- FALSE
  events$is_relative_gain_mean_clonal_nomsai <- FALSE
  events$is_absolute_gain_clonal_nomsai <- FALSE
  events$is_relative_gain_mean_subclonal_nomsai <- FALSE
  events$is_relative_gain_mean_subclonal_msai <- FALSE
  events$is_absolute_gain_subclonal_nomsai <- FALSE
  events$is_absolute_gain_subclonal_msai <- FALSE
  events$is_relative_loss_mean_clonal_nomsai <- FALSE
  events$is_absolute_loss_clonal_nomsai <- FALSE
  events$is_relative_loss_mean_subclonal_nomsai <- FALSE
  events$is_relative_loss_mean_subclonal_msai <- FALSE
  events$is_absolute_loss_subclonal_nomsai <- FALSE
  events$is_absolute_loss_subclonal_msai <- FALSE
  events$is_relative_amplification_ttest_clonal_nomsai <- FALSE
  events$is_relative_amplification_ttest_subclonal_nomsai <- FALSE
  events$is_relative_amplification_ttest_subclonal_msai <- FALSE
  events$is_relative_gain_ttest_clonal_nomsai <- FALSE
  events$is_relative_gain_ttest_subclonal_nomsai <- FALSE
  events$is_relative_gain_ttest_subclonal_msai <- FALSE
  events$is_relative_loss_ttest_clonal_nomsai <- FALSE
  events$is_relative_loss_ttest_subclonal_nomsai <- FALSE
  events$is_relative_loss_ttest_subclonal_msai <- FALSE
  events$is_LOH_clonal_nomsai <- FALSE
  events$is_LOH_subclonal_nomsai <- FALSE
  events$is_LOH_subclonal_msai <- FALSE
  events$is_ai_clonal_nomsai <- FALSE
  events$is_ai_subclonal_nomsai <- FALSE
  events
}
accounting_for_msai <- function(x, events_df) {
  clonalcall <- paste0(x, "_clonal")
  subclonalcall <- paste0(x, "_subclonal")
  clonalnomsaicall <- paste0(x, "_clonal_nomsai")
  subclonalnomsaicall <- paste0(x, "_subclonal_nomsai")
  subclonalmsaicall <- paste0(x, "_subclonal_msai")
  for (j in seq_len(nrow(events_df))) {
    if ((events_df$is_msai[j] == FALSE) & (events_df[, clonalcall][j] == TRUE) & (!is.na(events_df[, clonalcall][j]))) {
      events_df[, clonalnomsaicall][j] <- TRUE
    } else if (x != "is_ai") {
      if ((events_df[, paste0(x, "_parallel")][j] == FALSE) & (events_df[, subclonalcall][j] == TRUE) & (!is.na(events_df[, subclonalcall][j]))) {
        events_df[, subclonalnomsaicall][j] <- TRUE
      }
    } else if ((x == "is_ai") & (events_df[, "is_msai"][j] == FALSE) & (events_df[, subclonalcall][j] == TRUE) & (!is.na(events_df[, subclonalcall][j]))) {
      events_df[, subclonalnomsaicall][j] <- TRUE
    }
  }
  events_df
}
combine_loss_and_LOH_calls <- function(events) {
  events$is_relative_anyloss_mean_clonal <- FALSE
  events$is_relative_anyloss_mean_clonal_nomsai <- FALSE
  events$is_relative_anyloss_mean_subclonal <- FALSE
  events$is_relative_anyloss_mean_subclonal_nomsai <- FALSE
  events$is_relative_anyloss_mean_subclonal_msai <- FALSE
  events$is_relative_anyloss_ttest_clonal <- FALSE
  events$is_relative_anyloss_ttest_clonal_nomsai <- FALSE
  events$is_relative_anyloss_ttest_subclonal <- FALSE
  events$is_relative_anyloss_ttest_subclonal_nomsai <- FALSE
  events$is_relative_anyloss_ttest_subclonal_msai <- FALSE
  events$is_absolute_anyloss_clonal <- FALSE
  events$is_absolute_anyloss_clonal_nomsai <- FALSE
  events$is_absolute_anyloss_subclonal <- FALSE
  events$is_absolute_anyloss_subclonal_nomsai <- FALSE
  events$is_absolute_anyloss_subclonal_msai <- FALSE
  events$is_relative_anyloss_mean_parallel <- FALSE
  events$is_relative_anyloss_ttest_parallel <- FALSE
  events$is_absolute_anyloss_parallel <- FALSE
  for (i in seq_len(nrow(events))) {
    events$is_relative_anyloss_mean_clonal[i] <- ifelse((events$is_relative_loss_mean_clonal[i] == TRUE) | (events$is_LOH_clonal[i] == TRUE), TRUE, FALSE)
    events$is_relative_anyloss_mean_clonal_nomsai[i] <- ifelse((events$is_relative_loss_mean_clonal_nomsai[i] == TRUE) | (events$is_LOH_clonal_nomsai[i] == TRUE), TRUE, FALSE)
    events$is_relative_anyloss_mean_subclonal[i] <- ifelse((events$is_relative_loss_mean_clonal[i] == FALSE) & (events$is_LOH_clonal[i] == FALSE) & ((events$is_relative_loss_mean_subclonal[i] == TRUE) | (events$is_LOH_subclonal[i] == TRUE)), TRUE, FALSE)
    events$is_relative_anyloss_mean_subclonal_nomsai[i] <- ifelse((events$is_relative_loss_mean_clonal_nomsai[i] == FALSE) & (events$is_LOH_clonal_nomsai[i] == FALSE) & ((events$is_relative_loss_mean_subclonal_nomsai[i] == TRUE) | (events$is_LOH_subclonal_nomsai[i] == TRUE)) & ((events$is_relative_loss_mean_subclonal_msai[i] == FALSE) | (events$is_LOH_subclonal_msai[i] == FALSE)), TRUE, FALSE)
    events$is_relative_anyloss_mean_subclonal_msai[i] <- ifelse((events$is_relative_loss_mean_clonal_nomsai[i] == FALSE) & (events$is_LOH_clonal_nomsai[i] == FALSE) & ((events$is_relative_loss_mean_subclonal_msai[i] == TRUE) | (events$is_LOH_subclonal_msai[i] == TRUE)), TRUE, FALSE)
    events$is_relative_anyloss_ttest_clonal[i] <- ifelse((events$is_relative_loss_ttest_clonal[i] == TRUE) | (events$is_LOH_clonal[i] == TRUE), TRUE, FALSE)
    events$is_relative_anyloss_ttest_clonal_nomsai[i] <- ifelse((events$is_relative_loss_ttest_clonal_nomsai[i] == TRUE) | (events$is_LOH_clonal_nomsai[i] == TRUE), TRUE, FALSE)
    events$is_relative_anyloss_ttest_subclonal[i] <- ifelse((events$is_relative_loss_ttest_clonal[i] == FALSE) & (events$is_LOH_clonal[i] == FALSE) & ((events$is_relative_loss_ttest_subclonal[i] == TRUE) | (events$is_LOH_subclonal[i] == TRUE)), TRUE, FALSE)
    events$is_relative_anyloss_ttest_subclonal_nomsai[i] <- ifelse((events$is_relative_loss_ttest_clonal_nomsai[i] == FALSE) & (events$is_LOH_clonal_nomsai[i] == FALSE) & ((events$is_relative_loss_ttest_subclonal_nomsai[i] == TRUE) | (events$is_LOH_subclonal_nomsai[i] == TRUE)) & ((events$is_relative_loss_ttest_subclonal_msai[i] == FALSE) | (events$is_LOH_subclonal_msai[i] == FALSE)), TRUE, FALSE)
    events$is_relative_anyloss_ttest_subclonal_msai[i] <- ifelse((events$is_relative_loss_ttest_clonal_nomsai[i] == FALSE) & (events$is_LOH_clonal_nomsai[i] == FALSE) & ((events$is_relative_loss_ttest_subclonal_msai[i] == TRUE) | (events$is_LOH_subclonal_msai[i] == TRUE)), TRUE, FALSE)
    events$is_absolute_anyloss_clonal[i] <- ifelse((events$is_absolute_loss_clonal[i] == TRUE) | (events$is_LOH_clonal[i] == TRUE), TRUE, FALSE)
    events$is_absolute_anyloss_clonal_nomsai[i] <- ifelse((events$is_absolute_loss_clonal_nomsai[i] == TRUE) | (events$is_LOH_clonal_nomsai[i] == TRUE), TRUE, FALSE)
    events$is_absolute_anyloss_subclonal[i] <- ifelse((events$is_absolute_loss_clonal[i] == FALSE) & (events$is_LOH_clonal[i] == FALSE) & ((events$is_absolute_loss_subclonal[i] == TRUE) | (events$is_LOH_subclonal[i] == TRUE)), TRUE, FALSE)
    events$is_absolute_anyloss_subclonal_nomsai[i] <- ifelse((events$is_absolute_loss_clonal_nomsai[i] == FALSE) & (events$is_LOH_clonal_nomsai[i] == FALSE) & ((events$is_absolute_loss_subclonal_nomsai[i] == TRUE) | (events$is_LOH_subclonal_nomsai[i] == TRUE)) & ((events$is_absolute_loss_subclonal_msai[i] == FALSE) | (events$is_LOH_subclonal_msai[i] == FALSE)), TRUE, FALSE)
    events$is_absolute_anyloss_subclonal_msai[i] <- ifelse((events$is_absolute_loss_clonal_nomsai[i] == FALSE) & (events$is_LOH_clonal_nomsai[i] == FALSE) & ((events$is_absolute_loss_subclonal_msai[i] == TRUE) | (events$is_LOH_subclonal_msai[i] == TRUE)), TRUE, FALSE)
    events$is_relative_anyloss_mean_parallel[i] <- ifelse(((events$is_relative_loss_mean_parallel[i] == TRUE) | (events$is_LOH_parallel[i] == TRUE)), TRUE, FALSE)
    events$is_relative_anyloss_ttest_parallel[i] <- ifelse(((events$is_relative_loss_ttest_parallel[i] == TRUE) | (events$is_LOH_parallel[i] == TRUE)), TRUE, FALSE)
    events$is_absolute_anyloss_parallel[i] <- ifelse(((events$is_absolute_loss_parallel[i] == TRUE) | (events$is_LOH_parallel[i] == TRUE)), TRUE, FALSE)
  }
  events
}
