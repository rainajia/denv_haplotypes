## check all variants fed into regresshaplo
library('vcfR')
library('tidyverse')



## n=10 =-rect samples only
## n=20
list_indiv = c(1207,1223,1238,1276,1278,1280,1282,1284,1294,1307,1310,1312,1315,1324,1325,1338,1364,1378,1388,2430,2432,2433,2434,2436,2442,2443,2444,2447,2451,2452,2453,2454,2455,2456,2457,2460,2462,2463,2464,2467,2468,2469,2470,2471,2472,2477,2478,2480,2481,2482,2483,2909,2910,2911,2912,2913,2914,2929,2930,2931,2935,3023) # 20053 20054 20055 20056 20059 20060 20061 20062 2006$length(list_indiv)
length(list_indiv)

##################################################
## load csv
###################################################

input_dir = 'C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/regresshaplo_out_snv03/' # CHECK snv05 or ssnv03
setwd('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/regresshaplo_out_snv03') # CHECK snv05 or ssnv03

for( i in 1:length(list_indiv)){
  if(i==1){
    check = paste0(input_dir,list_indiv[i],'/final_haplo.csv')
    if (file.exists(check)=='TRUE'){
      df = as.data.frame(read.csv(paste0(input_dir,list_indiv[i],'/variant_calls.csv'),header=TRUE))
      df = df[,1:5]
      df$id = list_indiv[i]
    }else{
      df = data.frame(pos=0, A=0, C=0, G=0,T=0,id=list_indiv[i])
    }
  } else if (i>1){
    check = paste0(input_dir,list_indiv[i],'/final_haplo.csv')
    if (file.exists(check)!='TRUE'){
      df_tmp = data.frame(pos=0, A=0, C=0, G=0,T=0,id=list_indiv[i])
      df = rbind(df,df_tmp)
    }else{
    df_tmp = as.data.frame(read.csv(paste0(input_dir,list_indiv[i],'/variant_calls.csv'),header=TRUE))
    df_tmp = df_tmp[,1:5]
    df_tmp$id = list_indiv[i]
    df = rbind(df,df_tmp)
    }
  }
}


unique(df$id)


write.csv(df, "all_snv_to_regress.csv")

#write.csv(file = "snv_input_to_regresshaplo.csv", x=m)

#dir_snv = df  # view the number of SNVs
#iso_snv = df  # view the numbero f SNVs

table(df$id)   # count the number of haplotypes 

head(df)
nrow(df)

## check SNV numbers 
## check sample without snv
df_no_snv = df[which(df$pos == 0),]
df_no_snv = data.frame(table(df_no_snv$id))
write.csv(df_no_snv, "df_no_snv.csv")

#check sample with SNV
df_with_snv = df[which(df$pos > 0),]
df_with_snv = data.frame(table(df_with_snv$id))
write.csv(df_with_snv, "df_with_snv.csv")





