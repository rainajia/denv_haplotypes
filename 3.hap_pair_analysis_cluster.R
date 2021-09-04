# Haplotype pairs and relevant analysis  
# pair n = 25
# 02-03-2021
# Check word documents for Henrik's comment 
###################################################
# Needs datagrames created in :
#     1.check_hap.R
####################################################
# dataframe catalogue:
# hap_analysis  - contains indiv. haplotype frequency 
# shared_haplo_index - haplotype pairs with full id  - created in "1.check_haps"
# df_new - indiv.id, hap frequency, total hap no. 
#####################################################

# Analysis 1. For each individual in a pair, what proportion of its haplotypes are shared? 
# use df_new and shared_haplo_index dataframes as the starting point
# make a dataframe "hap_pair"

# task: asign a"id-hapID" to each individual, match with the pairs to fine which "id-hapID" shares haplotypes and with which "id-hap"
head(df_new)

hap_pair = df_new
hap_pair$unique_id  = paste0(hap_pair$id, "-", hap_pair$hap_index)    # replace unique_id column to "id-hapID" format (id = individual id, hapID= index of haplotype, i.e. haplotype 1 has hapID = 1) 
####################################################################################
####################################################################################
# 1. data cleaning to prepare the right dataframe
# re-arrange the shared_haplo_index dataframe 
####################################################################################
####################################################################################
head(shared_haplo_index)
#view(shared_haplo_index)
shared_hap_tmp = shared_haplo_index[,c(1,4)]
#head(shared_hap_tmp)

shared_hap_tmp[c("ID_A","delete1","hap_ind_A","delete2", "freq_A")] = t(data.frame(str_split(shared_hap_tmp$newColName, "_" ))) # solutions for how to create multiple columns from a text column with a delimiter 

shared_hap_tmp[c("ID_B","delete3","hap_ind_B","delete4", "freq_B")] = t(data.frame(str_split(shared_hap_tmp$name_col, "_" ))) 


shared_hap_tmp = shared_hap_tmp[,c(1,2,3,5,7,8,10,12)] #clean

shared_hap_tmp$pair_A_id = paste0(shared_hap_tmp$ID_A,"-",shared_hap_tmp$hap_ind_A) # create a unique id for the haplo (id + haplo index) from the individual (pair A
shared_hap_tmp$pair_B_id = paste0(shared_hap_tmp$ID_B,"-",shared_hap_tmp$hap_ind_B) # create a unique id for the haplo (id + haplo index) from the individual pair B

head(shared_hap_tmp)
####################################################################################
####################################################################################
#  Analysis 1: 
#  Qs: For each individual in a pair, what proportion of its haplotypes are shared?
####################################################################################
####################################################################################
# calculate the no. of times each hap id (id+hap index) is shared with another pair
#head(shared_hap_tmp)    #dataframe needed
#head(hap_pair)    # dataframe needed
shared_hap_tmp$hap_id = paste0(shared_hap_tmp$pair_A_id,"-",shared_hap_tmp$pair_B_id)
head(shared_hap_tmp)
nrow(shared_hap_tmp)
## match total hap no. to A and B, respectively 
#view(hap_analysis)    # from "1.check_hap.R"

## merge(dat1, dat2, by.x = 2, by.y = 0, all.x = TRUE)
shared_hap_tmp2 = merge(shared_hap_tmp,hap_analysis[,c(1,2)], by.x = "ID_A", by.y="id", all.x = TRUE)   # add total hap number for A
nrow(shared_hap_tmp2)   
colnames(shared_hap_tmp2)[12] = "hap_no_A"   # clean up the column names

## Then fine the total no. of hap for B
shared_hap_tmp2=merge(shared_hap_tmp2,hap_analysis[,c(1,2)], by.x = "ID_B", by.y="id", all.x = TRUE)   # add total hap number for B
colnames(shared_hap_tmp2)[13] = "hap_no_B"      # clean up the column names


## clean the freq data in the table (the original freq from the tree has "xxx.1" suffix to indicate duplicate of the same haplotypes in the haplotype pair table)
#head(df_new)   # all the original hap_freq data is in df_new from "1.check_hap"
shared_hap_tmp2 = subset(shared_hap_tmp2, select = -c(freq_A,freq_B))    # delete the original freq column with wrong format 

## use df_new to match the frequency 
df_new$ID_hap_ID=paste0(df_new$id,"-",df_new$hap_index)    # create an ID to match with shared_hap

shared_hap = merge(shared_hap_tmp2, df_new[,c(6,4)], by.x = "pair_A_id", by.y="ID_hap_ID", all.x = TRUE)  # re-match hap_freq from df_new to shared_hap
index_A = ncol(shared_hap)
colnames(shared_hap)[index_A] = "hap_freq_A"

shared_hap = merge(shared_hap, df_new[,c(6,4)], by.x = "pair_B_id", by.y="ID_hap_ID", all.x = TRUE) 
index_B = ncol(shared_hap)
colnames(shared_hap)[index_B] = "hap_freq_B"

## calculate the max. no. of hap to share between A and B
shared_hap = transform(shared_hap,max_hap_to_share = pmin(shared_hap$hap_no_A,shared_hap$hap_no_B))

view(shared_hap)
colnames(shared_hap)

####################################################################################
# key step: 
# for each individual in a pair, find the proportion of its haplotype
####################################################################################
nrow(shared_hap)

shared_hap$aux_pair = paste0(shared_hap$ID_A,"-",shared_hap$ID_B)
shared_hap$aux_pair_reverse = paste0(shared_hap$ID_B,"-",shared_hap$ID_A)

shared_hap$aux_pair[1]

for (i in (1:nrow(shared_hap))){
  
  # search for the existence of a pair of A and B (with different haplotypes)

  shared_hap$aux_pair_freq[i] = nrow(shared_hap[which(shared_hap$aux_pair == shared_hap$aux_pair[i] | shared_hap$aux_pair_reverse == shared_hap$aux_pair[i]  
                                                   | shared_hap$aux_pair == shared_hap$aux_pair_reverse[i] | shared_hap$aux_pair_reverse == shared_hap$aux_pair_reverse[i] ),])
}

# final step: 
# calculate the proportion of haplotype shared
shared_hap$porportion_of_shared = round(as.numeric(shared_hap$aux_pair_freq)/as.numeric(shared_hap$max_hap_to_share),2)

shared_hap
setwd('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/Tree/hap_analysis_snv01_final/')
write.csv(shared_hap, file = "porportion_hap_shared_between_pairs.csv") 


####################################################################################
####################################################################################
#  Analysis 2: 
#  Qs: Are the proportion of each haplotype correlated between linked pairs? 
#  (calculate by means of absolute difference)
####################################################################################
####################################################################################
# step 1:
# calculate the absolute difference

shared_hap$abs_diff = abs(as.numeric(shared_hap$hap_freq_A)-as.numeric(shared_hap$hap_freq_B))

head(shared_hap[,c("hap_freq_A","hap_freq_B","abs_diff")])

# step 2:
# look at the summary and histogram of the absolute difference between pairs
summary(shared_hap$abs_diff)

shared_hap_median = quantile(shared_hap$abs_diff, probs = 0.5)
shared_hap_mean = mean(shared_hap$abs_diff)

#shared_hap= shared_hap[which(shared_hap$ID_A != "1388"),] # get rid of unwanted pairs
#nrow(shared_hap)

# present the histogram of results
ggplot(shared_hap, aes(x=abs_diff)) + geom_histogram(binwitch =0.5, colour = "grey") +
  labs(title = "Absolute freq differences between haplotype pairs", x="Absolute difference between haplotype pairs", y = "Frequency") 
      #+ geom_text(aes(label=hap_freq), vjust=-0.5, color="black", position = position_dodge(1.5), size=3.5)


summary(shared_hap$abs_diff)


setwd('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/Tree/hap_analysis_snv03')#CHECK
write.csv(shared_hap, "absolute_freq_diff_between_pairs.csv")
####################################################################################
## compute the null expectations 
## i.e. the mean absolute difference if randomly pick one pair of haplotypes from all possible
## haplotype pairs regardless of their genetic distance
####################################################################################
#class(distance_matrix)  # use the distance matrix created in "1.check_hap"

## step 1:
## create a dataframe of all pairs from the distance matrix 
all_pairs = as.data.frame(as.table(distance_matrix))

# nrow(all_pair)     # check
#summary(all_pair$Freq)   # check 
#colnames(all_pair)
colnames(all_pairs) = c("unique_id_A", "unique_id_B","distance")

## split the haplotype frequency into individual columns
all_pairs[c("ID_A","delete1","hap_index_A","delete2", "freq_A")] = t(data.frame(str_split(all_pairs$unique_id_A, "_" ))) # solutions for how to create multiple columns from a text column with a delimiter 
all_pairs[c("ID_B","delete3","hap_index_B","delete4", "freq_B")] = t(data.frame(str_split(all_pairs$unique_id_B, "_" ))) 
all_pairs = subset(all_pairs, select = -c(delete1,delete2,delete3,delete4)) # clean up - delete the "delete" columns
#head(all_pairs)
nrow(all_pairs)    #14112

## step 2:
## bootstrap the data to calculate condifence intervals 

## set the parameters
n_boot = 100    #number of bootstrap iterations (the more the better, 200 is generally a good number)
boot_vect = rep(NA, n_boot) #set up a vector to store the results from each bootstrap iteration

for(i in 1:n_boot){ # loop through the iterations
  tmp = sample(1:nrow(all_pairs), 1) ## randomly draw one row
  boot_vect[i] = abs(as.numeric(all_pairs$freq_A[tmp]) - as.numeric(all_pairs$freq_B[tmp]))
}

#summary(boot_vect)

## Compute the mean and cis 
boot.median = quantile(boot_vect, probs = c(0.5), na.rm = T) # median
boot.mean = mean(boot_vect)

## check results
boot.median
boot.mean


## t-test to compare two means
shared_hap_t=shared_hap[which(shared_hap$abs_diff!=0),"abs_diff"] # don't include 0 in mean calculation
mean(shared_hap_t)

null_exp_t=boot_vect

t.test(shared_hap_t,null_exp_t)

#results of t-test; boot_n=200
#mean of x mean of y 
#0.3161893 0.3418450
#p-value = 0.3229
#95%CI = -0.07  0.02

#results of t-test; boot_n=5000
#mean of x mean of y 
#0.3161893 0.3462114 
#p-value = 0.101
#95%CI
#-0.065  0.005

# check plot
plot(boot_vect)   # completely random

####################################################################################
## output result
####################################################################################

## Analysis 2:output 

table = format(data.frame(boot.median,boot.mean,shared_hap_median,shared_hap_mean),digit=3)
colnames(table) = c("median_freq_null", "mean_freq_null", "median_full_pair", "mean_hap_pair")
#head(table)
write.csv(table, file ="C:/Users/micro/OneDrive - University of Cambridge/Documents/MPhil_project/trial_n67/haplotype_analysis_n62/hap_freq.csv")

# possible function to compute the minimum distance between two nodes 
#install.packages("igraph")
#?graph_from_data_frame


# Analysis 1:output table for the porportion of hap shraed 
colnames(shared_hap)
#head(shared_hap)
out = shared_hap[,c(2,1,6,14,18)]
nrow(out)


write.csv(out, file ="C:/Users/micro/OneDrive - University of Cambridge/Documents/MPhil_project/trial_n67/haplotype_analysis_n62/shared_hap_porportion.csv")



