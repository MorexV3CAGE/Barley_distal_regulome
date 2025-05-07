# NETCAGE directionality analysis
library(rtracklayer)

# Load NETCAGE dataset
data <- read.table("path/to/NETCAGE/data.bed", header=TRUE, stringsAsFactors=FALSE)

# Load bigwig files
bw_plus <- readBigwig("path/to/NETCAGE/bigwig_plus_strand.bw")
bw_minus <- readBigwig("path/to/NETCAGE/bigwig_minus_strand.bw")
data <- data[,c(1:6, 7, 11)]
colnames(data) <- c("seqnames", "start", "end", "geneID", "TPM", "strand", "domTSS", "domTPM")

# Initialize data frames to store results
both_signals <- data.frame()
plus_signal_only <- data.frame()
minus_signal_only <- data.frame()

# Iterate over each candidate
for (i in 1:nrow(data)) {
  enhancer_row <- data[i, ]
  
  # Check for signal in bw_plus
  plus_signal <- check_signal(enhancer_row, bw_plus)
  
  # Check for signal in bw_minus
  minus_signal <- check_signal(enhancer_row, bw_minus)
  
  # Categorize rows
  if (plus_signal && minus_signal) {
    enhancer_row$domTSS_plus <- get_dominant_TSS(enhancer_row, bw_plus)
    enhancer_row$domTSS_minus <- get_dominant_TSS(enhancer_row, bw_minus)
    if(enhancer_row$domTSS_plus <= enhancer_row$domTSS_minus){
      enhancer_row$start <- enhancer_row$domTSS_plus
      enhancer_row$end <- enhancer_row$domTSS_minus
    }
    else{
      enhancer_row$start <- enhancer_row$domTSS_minus
      enhancer_row$end <- enhancer_row$domTSS_plus
    }
    both_signals <- rbind(both_signals, enhancer_row)
  } else if (plus_signal) {
    enhancer_row$domTSS_plus <- get_dominant_TSS(enhancer_row, bw_plus)
    plus_signal_only <- rbind(plus_signal_only, enhancer_row)
  } else if (minus_signal) {
    enhancer_row$domTSS_minus <- get_dominant_TSS(enhancer_row, bw_minus)
    minus_signal_only <- rbind(minus_signal_only, enhancer_row)
  }
}

write.table(both_signals, file="output/path/bidirectional.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(plus_signal_only, file="output/path/unidirectional_plus.bed", quote=F, sep="\t", row.names=F, col.names=F)
write.table(minus_signal_only, file="output/path/unidirectional_minus.bed", quote=F, sep="\t", row.names=F, col.names=F)

# Iterate over each bidirectional row and check if the signals are balanced
both_signals_balanced <- data.frame()

for (i in 1:nrow(both_signals)) {
  enhancer_row <- both_signals[i, ]

  if (check_signal_balance(enhancer_row, bw_plus, bw_minus)) {
    both_signals_balanced <- rbind(both_signals_balanced, enhancer_row)
  }
  
}

write.table(both_signals_balanced, 
            file="output/path/bidirectional_balanced.bed", quote=F, sep="\t", row.names=F, col.names=F)


#---------------------------------------------------------------------------
# Used functions
#---------------------------------------------------------------------------
# Check if there is plus and minus signal for each row of our dataset
check_signal <- function(enhancer_row, signal_df) {
  enhancer_domTSS <- enhancer_row$domTSS
  
  # Check if there is any signal within the radius of 600bp from domTSS
  signal_within_enhancer <- signal_df[signal_df$seqnames == enhancer_row$seqnames &
                                        ((signal_df$start >= enhancer_domTSS-600 & signal_df$start <= enhancer_domTSS+600)), ]
  # Filter similar to CAGEr
  if (nrow(signal_within_enhancer) == 1) {
    # If there is only one row, keep it if abs(score) > 5
    if (abs(signal_within_enhancer$score) > 5) {
      signal_within_enhancer <- signal_within_enhancer
    } else {
      signal_within_enhancer <- data.frame()
    }
  } else {
    # If there are two or more rows, keep them if the sum of abs(score) > 0.1
    if (sum(abs(signal_within_enhancer$score)) > 0.1) {
      signal_within_enhancer <- signal_within_enhancer
    } else {
      signal_within_enhancer <- data.frame()
    }
  }
  
  return(nrow(signal_within_enhancer) > 0)
}

# Get dominant TSS position (highest TPM) from the bigwig file
get_dominant_TSS <- function(enhancer_row, signal_df) {
  enhancer_domTSS <- enhancer_row$domTSS
  
  # Subset the signals within the region
  signal_subset <- signal_df[signal_df$seqnames == enhancer_row$seqnames &
                               ((signal_df$start >= enhancer_domTSS-600 & signal_df$start <= enhancer_domTSS+600)), ]
  
  # Get the dominant TSS
  if (nrow(signal_subset) > 0) {
    max_tpm_row <- signal_subset[which.max(abs(signal_subset$score)), ]
    return(c(max_tpm_row$start))
  } else {
    return(c(NA))
  }
}


# Function to check if the signals on both strands are balanced
check_signal_balance <- function(enhancer_row, bw_plus, bw_minus) {
  enhancer_start <- enhancer_row$start
  enhancer_end <- enhancer_row$end
  
  # Extract the signals within the bidirectional region (and +/- 50 bp around)
  signal_plus <- bw_plus[bw_plus$seqnames == enhancer_row$seqnames &
                           bw_plus$start >= (enhancer_start - 50) &
                           bw_plus$end <= (enhancer_end + 50), ]
  
  signal_minus <- bw_minus[bw_minus$seqnames == enhancer_row$seqnames &
                             bw_minus$start >= (enhancer_start - 50) &
                             bw_minus$end <= (enhancer_end + 50), ]
  
  # Sum the scores (signal values) for plus and minus strands
  total_plus_signal <- sum(signal_plus$score)
  total_minus_signal <- sum(abs(signal_minus$score))
  
  # Check for balance
  log_fc <- log2(total_plus_signal / total_minus_signal)
  if(log_fc > 0){
    return(log_fc < 1)
  }
  else{
    return(log_fc > -1)
  }
}



