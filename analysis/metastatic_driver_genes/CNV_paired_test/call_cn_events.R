call_cn_events <- function(segs, snps, sample_data) {
  segs_with_thresholds <- calc_logr_thresholds(segs, sample_data)
  for (i in seq_len(nrow(segs_with_thresholds))) {
    snps_sample <- as.data.frame(snps[[segs_with_thresholds$sample_id[i]]])[, c("seqnames", "pos", "logr")]
    segment_snps <- snps_sample[which((as.character(snps_sample$seqnames) == as.character(segs_with_thresholds$seqnames[i])) &
      (snps_sample$pos >= segs_with_thresholds$start[i]) & (snps_sample$pos <= segs_with_thresholds$end[i])), ]
    segs_with_thresholds$meanlogR[i] <- mean(segment_snps$logr)
    segs_with_thresholds$is_relative_amplification_mean[i] <- (segs_with_thresholds$meanlogR[i] >= segs_with_thresholds$amplification_logRthreshold[i])
    segs_with_thresholds$is_relative_gain_mean[i] <- (segs_with_thresholds$meanlogR[i] >= segs_with_thresholds$gain_logRthreshold[i])
    segs_with_thresholds$is_relative_loss_mean[i] <- (segs_with_thresholds$loss_logRthreshold[i] >= segs_with_thresholds$meanlogR[i])
    if (nrow(segment_snps) > 1) {
      relative_amplification_ttest_pval <- t_test_with_error_handling(segment_snps$logr, mu = segs_with_thresholds$amplification_logRthreshold[i], alternative = "greater")$p.value
      relative_gain_ttest_pval <- t_test_with_error_handling(segment_snps$logr, mu = segs_with_thresholds$gain_logRthreshold[i], alternative = "greater")$p.value
      relative_loss_ttest_pval <- t_test_with_error_handling(segment_snps$logr, mu = segs_with_thresholds$loss_logRthreshold[i], alternative = "less")$p.value
      segs_with_thresholds$is_relative_amplification_ttest[i] <- (relative_amplification_ttest_pval < 0.01)
      segs_with_thresholds$is_relative_gain_ttest[i] <- (relative_gain_ttest_pval < 0.01)
      segs_with_thresholds$is_relative_loss_ttest[i] <- (relative_loss_ttest_pval < 0.01)
    } else {
      segs_with_thresholds$is_relative_amplification_ttest[i] <- NA
      segs_with_thresholds$is_relative_gain_ttest[i] <- NA
      segs_with_thresholds$is_relative_loss_ttest[i] <- NA
    }
    segs_with_thresholds$is_absolute_amplification[i] <- ((segs_with_thresholds$cn_a[i] + segs_with_thresholds$cn_b[i]) >= 4)
    segs_with_thresholds$is_absolute_gain[i] <- ((segs_with_thresholds$cn_a[i] + segs_with_thresholds$cn_b[i]) >= 2.5)
    segs_with_thresholds$is_absolute_loss[i] <- (1.5 >= (segs_with_thresholds$cn_a[i] + segs_with_thresholds$cn_b[i]))
    segs_with_thresholds$is_homozygous_deletion[i] <- ifelse((round(segs_with_thresholds$cn_b[i]) <= 0) & (round(segs_with_thresholds$cn_a[i]) <= 0), TRUE, FALSE)
    segs_with_thresholds$negative_cn_called[i] <- ifelse((segs_with_thresholds$cn_a[i] < 0 | segs_with_thresholds$cn_b[i] < 0), TRUE, FALSE)
  }
  return(segs_with_thresholds)
}
