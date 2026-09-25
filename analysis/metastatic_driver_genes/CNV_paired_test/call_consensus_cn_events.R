call_consensus_cn_events <- function(phased_segs_cn_event_calls, sample_number) {
  consensus_events <- create_consensus_dataframe(phased_segs_cn_event_calls)
  events_vector <- c(
    "is_relative_amplification_mean", "is_relative_gain_mean", "is_relative_loss_mean",
    "is_relative_amplification_ttest", "is_relative_gain_ttest", "is_relative_loss_ttest",
    "is_absolute_amplification", "is_absolute_gain", "is_absolute_loss",
    "is_ai", "is_LOH", "is_homozygous_deletion"
  )
  for (event_type_entry in events_vector) {
    consensus_events <- calc_consensus(
      event_type = event_type_entry,
      sample_cn_events = phased_segs_cn_event_calls,
      events_df = consensus_events,
      num_samples = sample_number
    )
  }
  cn_event_category_cols_df <- data.frame(
    relative_mean_test = c(
      "is_relative_amplification_mean",
      "is_relative_gain_mean",
      "is_relative_loss_mean"
    ),
    relative_ttest = c(
      "is_relative_amplification_ttest",
      "is_relative_gain_ttest",
      "is_relative_loss_ttest"
    ),
    absolute = c(
      "is_absolute_amplification",
      "is_absolute_gain",
      "is_absolute_loss"
    ),
    stringsAsFactors = FALSE
  )
  consensus_events <- calc_consensus_parallel(modifiable_events_df = consensus_events, phased_segs_cn_event_calls, cn_event_category_cols_df)
  consensus_events <- initialise_msai_columns(consensus_events)
  events_vector_msai <- c(
    "is_relative_amplification_mean", "is_relative_gain_mean", "is_relative_loss_mean",
    "is_relative_amplification_ttest", "is_relative_gain_ttest", "is_relative_loss_ttest",
    "is_absolute_amplification", "is_absolute_gain", "is_absolute_loss",
    "is_ai", "is_LOH"
  )
  for (i in events_vector_msai) {
    consensus_events <- accounting_for_msai(x = i, events_df = consensus_events)
  }
  consensus_events <- combine_loss_and_LOH_calls(consensus_events)
  consensus_events
}
