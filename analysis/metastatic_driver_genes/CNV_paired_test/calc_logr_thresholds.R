calc_logr_thresholds <- function(segs, sample_data) {
  segs_df <- as.data.frame(segs)
  segs_df_subset <- segs_df[, c("group_name", "seqnames", "start", "end", "width", "cn_a", "cn_b", "is_ai", "mirrored_vs_ref", "is_LOH")]
  sample_data_subset <- sample_data[, c("sample_id", "purity", "ploidy")]
  segs_df_subset <- merge(sample_data_subset, segs_df_subset, by.x = "sample_id", by.y = "group_name", all.y = TRUE)
  segs_df_subset$amplification_logRthreshold <- log2((2 * (1 - segs_df_subset$purity) + (segs_df_subset$purity * segs_df_subset$ploidy * 2)) / (2 * (1 - segs_df_subset$purity) + (segs_df_subset$purity * segs_df_subset$ploidy)))
  segs_df_subset$gain_logRthreshold <- log2((2 * (1 - segs_df_subset$purity) + (segs_df_subset$purity * segs_df_subset$ploidy * 1.25)) / (2 * (1 - segs_df_subset$purity) + (segs_df_subset$purity * segs_df_subset$ploidy)))
  segs_df_subset$loss_logRthreshold <- log2((2 * (1 - segs_df_subset$purity) + (segs_df_subset$purity * segs_df_subset$ploidy * 0.75)) / (2 * (1 - segs_df_subset$purity) + (segs_df_subset$purity * segs_df_subset$ploidy)))
  segs_df_subset
}
