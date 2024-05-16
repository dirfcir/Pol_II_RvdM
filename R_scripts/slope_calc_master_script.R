calculate_slope_master_funct <- function(species                 = args[1], 
                                         sj_tab_folder           = args[2], 
                                         gtf_file                = args[3], 
                                         outfile_defined_introns = args[4], 
                                         outfile_slope_calc      = args[5]){
  define_introns(outfile_defined_introns = outfile_defined_introns, gtf_file = gtf_file)
  sj_annot()
}


