#########################################################
### FOR KEEPS FOR KEEPS FOR KEEPS KOR KEEPS FOR KEEPS ###
#########################################################
compute_y_intercept <- function(x1, filt = TRUE) {
  if (is.null(x1)) { return(NULL) }
  slope <- sapply(1:length(x1), function(z) { sum <- lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]]))); return(sum$fit[[1]]) })
  return(as.vector(slope))
}

compute_x_intercept <- function(x1, filt = TRUE) {
  filt_val <- 0
  read_size <- 90
  
  if (is.null(x1)) { return(NULL) }
  
  slope <- sapply(1:length(x1), function(z) { sum <- lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]]))); return(sum$fit[[length(x1[[z]])]]) })
  return(as.vector(slope))
}

get_chr_coverage <- function(x, chr, cov_data) {
  chrs1 <- cov_data[[x]][seqnames(cov_data[[x]]) == chr]
  covs1 <- coverage(chrs1, weight = chrs1$score)[[chr]]
  return(covs1)
}

compute_slope <- function(x1, filt = TRUE) {            
  if (is.null(x1)) { return(NULL) }
  slope <- sapply(1:length(x1), function(z) { sum <- summary(lm(as.numeric(x1[[z]]) ~ seq(1:length(x1[[z]])))); return(sum[4][[1]][[2]]) })
  return(slope)                  
}

filter_gene <- function(x, chr, strand_sense = '') {
  if (strand_sense == 'minus') {
    gff_exon <- x[(x$chr == chr) & (x$strand == "-"),]
  }
  if (strand_sense == 'plus') {
    gff_exon <- x[(x$chr == chr) & (x$strand == "+"),]
  }
  return(gff_exon)
}

paralel_import_bed_graphs <- function(bedgraph_filepaths, cluster) {
  imported_bed_graphs <- parLapply(cluster, bedgraph_filepaths, function(d) {
    imported_bed_graph <- import.bedGraph(d)
    #seqlengths(imported_bed_graph) <- chr_size[levels(seqnames(imported_bed_graph))]                       
    return(imported_bed_graph)
    
  }) 
  return(imported_bed_graphs)
}

get_cov_strand_spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, cl) { 
  
  # Get annotation
  gff_filt_fwd <- filter_gene(x = gff_b, chr = chr, chr_size_m = chrsize, strand_sense = "plus")
  gff_filt_rev <- filter_gene(x = gff_b, chr = chr, chr_size_m = chrsize, strand_sense = "minus")
  
  # Split by chromosome
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get_chr_coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get_chr_coverage, chr, cov_data_rev)
  
  # Extract genomic region of interest
  size_of_annot_fwd <- nrow(gff_filt_fwd)
  size_of_annot_rev <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x) {                                                       
    mcmapply(function(s, e) { window(cov_data_chr_fwd[[x]], s, e) }, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
  })
  
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x) {                                                       
    mcmapply(function(s, e) { rev(window(cov_data_chr_rev[[x]], s, e)) }, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
  })
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <- lapply(1:size_of_annot_fwd, function(y) lapply(r_fw, function(x) x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:size_of_annot_rev, function(y) lapply(r_rev, function(x) x[[y]]))
  
  gen_data <- c(r_fw, r_rev)
  return(gen_data)
}

run_chr_y_intercept_strand_spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl) { 
  
# Get annotation
gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")

# Split by chromosome
cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get_chr_coverage, chr, cov_data_fwd)        
cov_data_chr_rev <- lapply(1:length(cov_data_rev), get_chr_coverage, chr, cov_data_rev)
# Extract genomic region of interest
size_of_annot_fwd <- nrow(gff_filt_fwd)
size_of_annot_rev <- nrow(gff_filt_rev)

r_fw <- lapply(1:length(cov_data_chr_fwd), function(x) {                                                       
  mcmapply(function(s, e) { window(cov_data_chr_fwd[[x]], s, e) }, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
})

r_rev <- lapply(1:length(cov_data_chr_rev), function(x) {                                                       
  mcmapply(function(s, e) { rev(window(cov_data_chr_rev[[x]], s, e)) }, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
})

if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <- lapply(1:size_of_annot_fwd, function(y) lapply(r_fw, function(x) x[[y]]))

if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:size_of_annot_rev, function(y) lapply(r_rev, function(x) x[[y]]))

gen_data <- c(r_fw, r_rev)
gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)

cat("Get slope \n")
gen_slope <- mclapply(gen_data, compute_y_intercept, mc.cores=12)

gen_slope[gen_slope == "NULL"] <- NA

cat(paste(chr, "\n"))
gen_slope <- do.call(rbind, gen_slope)
gen_slope <- cbind(gff_filt, gen_slope)        
return(gen_slope)
}

run_chr_x_intercept_strand_spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, repli, chrsize, direc, cl) { 
  # Get annotation
  gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
  gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
  
  # Split by chromosome
  cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get_chr_coverage, chr, cov_data_fwd)        
  cov_data_chr_rev <- lapply(1:length(cov_data_rev), get_chr_coverage, chr, cov_data_rev)
  
  # Extract genomic region of interest
  size_of_annot_fwd <- nrow(gff_filt_fwd)
  size_of_annot_rev <- nrow(gff_filt_rev)
  
  r_fw <- lapply(1:length(cov_data_chr_fwd), function(x) {                                                       
    mcmapply(function(s, e) { window(cov_data_chr_fwd[[x]], s, e) }, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
  })
  
  r_rev <- lapply(1:length(cov_data_chr_rev), function(x) {                                                       
    mcmapply(function(s, e) { rev(window(cov_data_chr_rev[[x]], s, e)) }, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
  })
  
  if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <- lapply(1:size_of_annot_fwd, function(y) lapply(r_fw, function(x) x[[y]]))
  
  if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:size_of_annot_rev, function(y) lapply(r_rev, function(x) x[[y]]))
  
  gen_data <- c(r_fw, r_rev)
  gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
  
  # Get slope
  cat("Get slope \n")
  gen_slope <- mclapply(gen_data, compute_x_intercept, mc.cores=12)
  cat(paste(format(object.size(gen_slope), units = "MB"), "\n"))
  
  cat("alphab")
  gen_slope[gen_slope == "NULL"] <- NA
  cat(paste(chr, "\n"))
  gen_slope <- do.call(rbind, gen_slope)
  gen_slope <- cbind(gff_filt, gen_slope)        
  return(gen_slope)
}

run_chr_slope_strand_spe <- function(chr, cov_data_fwd, cov_data_rev, gff_b, chrsize, direc, cl) { 
    # Get annotation
    gff_filt_fwd <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "plus")
    gff_filt_rev <- filter_gene(gff_b, gen_elem, chr, chrsize, direc, "minus")
    
    # Split by chromosome
    cov_data_chr_fwd <- lapply(1:length(cov_data_fwd), get_chr_coverage, chr, cov_data_fwd)        
    cov_data_chr_rev <- lapply(1:length(cov_data_rev), get_chr_coverage, chr, cov_data_rev)
    
    # Extract genomic region of interest
    size_of_annot_fwd <- nrow(gff_filt_fwd)
    size_of_annot_rev <- nrow(gff_filt_rev)
    
    r_fw <- lapply(1:length(cov_data_chr_fwd), function(x) {                                                       
      mcmapply(function(s, e) { window(cov_data_chr_fwd[[x]], s, e) }, gff_filt_fwd$start, gff_filt_fwd$end, mc.cores = 12)                                      
    })
    
    r_rev <- lapply(1:length(cov_data_chr_rev), function(x) {                                                       
      mcmapply(function(s, e) { rev(window(cov_data_chr_rev[[x]], s, e)) }, gff_filt_rev$start, gff_filt_rev$end, mc.cores = 12)       
    })
    
    if (length(r_fw[[1]]) == 0) r_fw <- NULL else r_fw <- lapply(1:size_of_annot_fwd, function(y) lapply(r_fw, function(x) x[[y]]))
    if (length(r_rev[[1]]) == 0) r_rev <- NULL else r_rev <- lapply(1:size_of_annot_rev, function(y) lapply(r_rev, function(x) x[[y]]))
    
    gen_data <- c(r_fw, r_rev)
    gff_filt <- rbind(gff_filt_fwd, gff_filt_rev)
    
    # Get slope
    cat("Get slope \n")
    gen_slope <- mclapply(gen_data, compute_slope, mc.cores = 12)
    cat(paste(format(object.size(gen_slope), units = "MB"), "\n"))
    
    cat("alphab")
    gen_slope[gen_slope == "NULL"] <- NA
    cat(paste(chr, "\n"))
    gen_slope <- do.call(rbind, gen_slope)
    gen_slope <- cbind(gff_filt, gen_slope)        
    return(gen_slope)
  }
  