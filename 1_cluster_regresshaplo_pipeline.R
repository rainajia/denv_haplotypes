## RegressHaplo_NMC_passage vs non-passage 
## Last version: 28-06-2021
######################################
# Preparation 
######################################
# set up regressHaplo
#install.packages(c("igraph", "plyr", "dplyr"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("Biostrings","Rsamtools","GenomicAlignments", "rmutil"))
#install.packages("devtools")
#install.packages("Rtools")

#library("igraph")
#library("plyr")
#library("dplyr")
#library("rmutil")
#library("Biostrings")
#library("Rsamtools")
#library("GenomicAlignments")
#library("devtools")
#devtools::install_github("SLeviyang/RegressHaplo")
library("RegressHaplo")

# end of regressHaplo set-up 
#install.packages("doParallel")
#install.packages("seqinr")
#install.packages("vcfR")
library("seqinr")
library("doParallel")
library("vcfR")


# module load pkg-config-0.29.2-gcc-6.2.0-we4glmw
# module load R/3.6
# module load gcc/9

## n=62
list_indiv = c(1207,1223,1238,1276,1278,1280,1282,1284,1294,1307,1310,1312,1315,1324,1325,1338,1364,1378,1388,2430,2432,2433,2434,2436,2442,2443,2444,2447,2451,2452,2453,2454,2455,2456,2457,2460,2462,2463,2464,2467,2468,2469,2470,2471,2472,2477,2478,2480,2481,2482,2483,2909,2910,2911,2912,2913,2914,2929,2930,2931,2935,3023) # 
length(list_indiv)

list_bam = list.files(path='/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/' , pattern ="*.bam", all.files = FALSE,
                      full.names = FALSE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE) ## Careful : takes everything with "trim_sorted.bam" is it, including *.bam.bai files !
length(list_bam)
index = seq(1, length(list_bam), 2)
list_bam = list_bam[index]
list_bam


list_vcf = list.files(path='/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/',pattern = "*.vcf", all.files = FALSE,
                      full.names = FALSE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#index = seq(1, length(list_vcf), 3)
#list_vcf = list_vcf[index]
list_vcf

list_ref = list.files(path='/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/',pattern = "*.fasta", all.files = FALSE,
                      full.names = FALSE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#print(list_ref)
index = seq(1, length(list_ref), 7)
list_ref = list_ref[index]
list_ref

## For loop for the VCF filtering and RegressHaplo pipeline
# Parallelization
cores <- 30
#cluster <- makeCluster(cores)
registerDoParallel(cores)

#.packages = c("Rcpp", "tibble", "dplyr", "RcppArmadillo")  

foreach(i = 1:30) %dopar% { #length(list_indiv)
  
    ## STEP 1 
    ## Read and filter own variant calling vcf files
    print(paste0('Indiv: ', i, '/', length(list_indiv)))
    
    bam_file = paste0('/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/', list_bam[i])
    out_dir = paste0('/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/run2_snv03/', list_indiv[i], "/") 
    
    print(paste0("bam file is",bam_file))
    print(paste0("Outdir file is", out_dir))
    
    ## own variant calling
    vcf = read.vcfR(paste0('/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/',list_vcf[i]), verbose = TRUE)
    info = INFO2df(vcf)   # For lofreq 
    vcf = cbind(data.frame(vcf@fix),info)    #fix field as a dataframe
    vcf$ind_id = list_indiv[i]
    
    print(paste0("vcf file created is ",list_vcf[i])) 
  
    # Filter SNVs by the DP4 field
    DP4 = data.frame("DP"=vcf$DP, "DP4"=vcf$DP4)  
    DP4_tmp = matrix(0, nrow = nrow(DP4), ncol = 4)
    
    for (j in (1:length(DP4$DP4))){
      DP4_tmp[j,] = as.numeric(unlist(strsplit(as.character(DP4$DP4[j]), ',')))[1:4]
    }
    
    # Filtering steps 
    vcf$DP4_REF = DP4_tmp[,1]+DP4_tmp[,2] #number of reads for the REF
    vcf$DP4_ALT = DP4_tmp[,3]+DP4_tmp[,4] #number of reads for the ALT (snp)
    vcf$DP4_ALT_F = DP4_tmp[,3] #number of reads for the REF
    vcf$DP4_ALT_R = DP4_tmp[,4] #number of reads for the ALT
    vcf$DP4_F = vcf$DP4_ALT_F+ DP4_tmp[,1] #number of F reads
    vcf$DP4_R = vcf$DP4_ALT_R + DP4_tmp[,2]#number of R reads
    vcf$sum = vcf$DP4_REF + vcf$DP4_ALT #total number of reads
    vcf$ALT_SB_F = vcf$DP4_ALT_F/vcf$DP4_ALT
    vcf$ALT_SB_R = vcf$DP4_ALT_R/vcf$DP4_ALT
  
    
    vcf = vcf[which(vcf$DP4_ALT/vcf$sum >=0.03),] # filtering on the frequency, no variant should have freq >50% in this case. 
    vcf = vcf[which(vcf$DP >= 1000),] # hard filtering on mean read depth;use  1st Qu.if caller is LoFreq
    vcf = vcf[which(vcf$SB < 60 & (vcf$ALT_SB_F<0.80 || vcf$ALT_SB_R<0.80)),] # strand bias 
    
    variant_number_remian = nrow(vcf)
    print(paste0(list_indiv[i], " has ", variant_number_remian," variants"))

    #####
    if(nrow(vcf) == 0){
      ref = getSequence(seqinr::read.fasta(paste0("/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/",list_ref[i])))
      setwd(out_dir)
      seqinr::write.fasta(ref[[1]], names = paste0(list_indiv[i],'_hap_1_freq_1'), file.out = paste0(list_indiv[i], '_hap_final.fasta'))
    }else if(nrow(vcf) > 0){     # if no variants left after filtering, save as fasta in RegressHaplo outout dir
      write.csv(file = paste0(out_dir,list_indiv[i],"_own_snv_filtered.csv"), x = vcf , row.names = FALSE)

    ## STEP 2
    ## compute the variant_call.csv file as an input/ into RegressHaplo
    bam_to_variant_calls.pipeline(bam_file, out_dir, start_pos=NULL, end_pos=NULL, sig=0.01, heavy_tail=T)
    
    ## Compare own variant calling with RegressHaplo variant caller 
    variant_caller = get_variant_calls.pipeline(out_dir) #returns the data.frame written in pileup.csv and variant_calls.csv
    choose = NULL
    for (k in vcf$POS){
      choose = c(choose, which(variant_caller$pos == k))
    }
    
    if(length(choose) == 0){
      ref = getSequence(seqinr::read.fasta("/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/",list_ref[i]))
      setwd(out_dir)
      write.fasta(ref[[1]], names = paste0(list_indiv[i],'_hap_1_freq_1'), file.out = paste0(list_indiv[i], '_hap_final.fasta'))
   }else if(length(choose) > 0){
    
    ## use own variant calling results
    variant_caller_filtered = variant_caller[choose,]
    variant_caller_filtered = variant_caller_filtered[which(variant_caller_filtered$pos!=10080 & variant_caller_filtered$pos!=10128 & variant_caller_filtered$pos!=10137 & variant_caller_filtered$pos!=2866),] ## exclude these positions
  
    ## check variant number change 
    row_own = nrow(variant_caller)
    row_filtered = nrow(variant_caller_filtered)
    print(paste0(list_indiv[i],"- Before filtering:", row_own, " After filtering:", row_filtered))
   
    if(nrow(variant_caller_filtered) == 0){     # if no variants left after filtering, save as fasta in RegressHaplo outout dir
      ref = getSequence(seqinr::read.fasta(paste0("/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/",list_ref[i])))
      setwd(out_dir)
      seqinr::write.fasta(ref[[1]], names = paste0(list_indiv[i],'_hap_1_freq_1'), file.out = paste0(list_indiv[i], '_hap_final.fasta'))
    }else if(nrow(variant_caller_filtered ) > 0){
      write.csv(file = paste0(out_dir, "variant_calls.csv"), x = variant_caller_filtered, row.names = FALSE)
    
    ## STEP 3
    ## continue regress haplo
    variant_calls_to_read_table.pipeline(bam_file, out_dir, use_raw_read_table=F, sig=0.01) #- outputs raw_read_table.csv and read_table.csv. ## SUPER LONG
    read_table_to_loci.pipeline(out_dir, max_num_haplotypes=1000) #- outputs loci.csv.
    loci_to_haplotypes.pipeline(out_dir, max_num_haplotypes=1000) #- outputs h.csv.
    haplotypes_to_parameters.pipeline(out_dir) #- outputs y.csv, P.csv
    parameters_to_solutions.pipeline(out_dir, num_trials=2400, rho_vals = c(0.1,0.5,1,2,3.5,5,7,10)) #- outputs solutions.csv ## SUPER LONG (take from 0,1 to 10 )
    solutions_to_haplotypes.pipeline(out_dir, K_val = NULL) #- outputs final_haplo.csv
    haplotypes_to_fasta.pipeline(bam_file, out_dir) #- outputs final_haplo.fasta
    }
   }
  }
    
    ## rename outputS

    ref = getSequence(read.fasta(paste0('/home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/regresshaplo/input/',list_ref[i])))
    setwd(out_dir)
    variants = read.csv('final_haplo.csv', header = F, colClasses = 'character')
    
    pos = get_variable_positions.pipeline(out_dir)
    
    names = paste0(list_indiv[i], '_hap_', 1:length(variants[,1]), '_freq_', as.character(round(as.numeric(as.character(variants[,1])), digits = 3)))
    
    seq = as.list(rep(NA, length(variants[,1])))
    for (j in 1:length(variants[,1])){
      tmp = ref[[1]]
      tmp[pos] = tolower(variants[j,2:length(variants[1,])])
      seq[[j]] = tmp
    }
    
    data = data.frame(seq[[1]])
    for (k in 2:length(variants[,1])){
      data = cbind(data, seq[[k]])
    }
    seqinr::write.fasta(data, names = names, file.out = paste0(list_indiv[i], '_hap_final.fasta'))
    
}









