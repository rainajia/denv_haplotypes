## check all variants fed into regresshaplo
library('vcfR')
library('tidyverse')

setwd('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/2nd_round_vcf')
setwd('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/1st_round_vcf')

## n=62
list_indiv = c(1207,1223,1238,1276,1278,1280,1282,1284,1294,1307,1310,1312,1315,1324,1325,1338,1364,1378,1388,2430,2432,2433,2434,2436,2442,2443,2444,2447,2451,2452,2453,2454,2455,2456,2457,2460,2462,2463,2464,2467,2468,2469,2470,2471,2472,2477,2478,2480,2481,2482,2483,2909,2910,2911,2912,2913,2914,2929,2930,2931,2935,3023)
length(list_indiv)

list_vcf = list.files(path='./', pattern = "*.vcf", all.files = FALSE,
                      full.names = FALSE, recursive = FALSE,
                      ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)  # NOTES: use the vcf before trimming or read frequency > 0.05 for pseudo-consensus
list_vcf

## load all vcf to a dataframe
## samples with no SNV will occupy one row with indiv. ID presetn in the firsit column
for (i in 1:length(list_indiv)){      
  
  print(paste0('Indiv: ', list_vcf[i]))
  if(i==1){
    vcf_df = read.vcfR(paste0('./', list_vcf[1])) 
    info = INFO2df(vcf_df)
    nrow=nrow(info)
    if(nrow==0){
      fix = as.data.frame(vcf_df@fix)
      fix[1,1] = ""
      rm(info)
      info= data.frame("DP"="","AF"="","SB"="","DP4"="","INDEL"="","CONSVAR"="","HRUN"="")
      vcf_df = cbind(fix,info) #vcf_df=vcf_df[,c(1:2,4:15)]  # get rid of the original ID column and assign a new ID col as the list_indiv id
      vcf_df =cbind(ID=list_indiv[1], vcf_df)
      vcf_df = vcf_df[,c("ID","CHROM","POS", "REF","ALT","QUAL","FILTER","INFO","DP","AF", "SB","DP4","INDEL","CONSVAR","HRUN")]
    }else if(nrow>0){
     vcf_df = cbind(data.frame(ID=list_indiv[i], vcf_df@fix),info) #fix field as a dataframe
     vcf_df = vcf_df[,c("ID","CHROM","POS", "REF","ALT","QUAL","FILTER","INFO","DP","AF", "SB","DP4","INDEL","CONSVAR","HRUN")]
    }
  }else if (i>1){
    vcf_tmp = read.vcfR(paste0('./', list_vcf[i])) 
    info = INFO2df(vcf_tmp)
    nrow=nrow(info)
    if(nrow==0){
      fix = as.data.frame(vcf_tmp@fix)
      fix[1,1] = ""
      rm(info)
      info= data.frame("DP"="","AF"="","SB"="","DP4"="","INDEL"="","CONSVAR"="","HRUN"="")
      vcf_tmp = cbind(fix,info)
      vcf_tmp2 =cbind(ID=list_indiv[i], vcf_tmp)
      vcf_tmp2 = vcf_tmp2[,c("ID","CHROM","POS", "REF","ALT","QUAL","FILTER","INFO","DP","AF", "SB","DP4","INDEL","CONSVAR","HRUN")]
      vcf_df=rbind(vcf_df, vcf_tmp2)
    }else if (nrow>0){
    vcf_tmp = cbind(data.frame(ID = list_indiv[i], vcf_tmp@fix),info) #fix field as a dataframe
    vcf_tmp = vcf_tmp[,c("ID","CHROM","POS", "REF","ALT","QUAL","FILTER","INFO","DP","AF", "SB","DP4","INDEL","CONSVAR","HRUN")]
    vcf_df = rbind(vcf_df, vcf_tmp)
    }
 }
}


write.csv(vcf_df,"./all_1st_round_snv.csv")

#######################################################################################################

vcf_df$ID = as.character(vcf_df$ID)
vcf_df$POS = as.numeric(vcf_df$POS)
vcf_df$DP = as.numeric(vcf_df$DP)
vcf_df$AF = as.numeric(vcf_df$AF)
vcf_df$QUAL = as.numeric(vcf_df$QUAL)

vcf_df_analysis = vcf_df[which(vcf_df$DP>=1000),] ## double-check on depth filter

minor_snv_pos <- ggplot(vcf_df_analysis, aes(x=POS, y=AF)) +
                 geom_point(color="turquoise3")+xlim(0,10100)+ylim(0,0.5) + 
                 geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
                 labs(x="Genome position", y="SNV frequency ")


minor_snv_DP <- ggplot(vcf_df_analysis, aes(x=DP, y=AF)) +
           geom_point(color="turquoise3")+xlim(0,10100)+ylim(0,0.5) + 
           geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=0.5) +
           labs(x="Depth", y="SNV frequency ")

ggarrange(minor_snv_pos,minor_snv_DP,
          common.legend = TRUE, legend = "right",
          ncol = 2, nrow = 1)









######################################################################################
## other descriptive anlaysis
##################################################################################

summary(vcf_df$DP)



ggplot(vcf_df, aes(x=ID, y=DP)) +
  geom_boxplot()

ggplot(vcf_df, aes(DP)) +
  geom_density(fill = "lightblue")

DP_cutoff = data.frame(
  "0.05" = quantile(vcf_df$DP,0.05,na.rm=TRUE),
  "0.95" = quantile(vcf_df$DP,0.95,na.rm=TRUE),
  "0.01" = quantile(vcf_df$DP,0.01,na.rm=TRUE),
  "0.99" = quantile(vcf_df$DP,0.99,na.rm=TRUE),
  "0.005" = quantile(vcf_df$DP,0.005,na.rm=TRUE),
  "0.995" = quantile(vcf_df$DP,0.995,na.rm=TRUE))
DP_cutoff
summary(vcf_df$DP)

DP_cutoff = vcf_df[which(vcf_df$DP<1000),]
table(DP_cutoff$POS)
DP_cutoff$AF


## average depth by position 
dp_by_pos = aggregate(DP~POS,data=vcf_df,FUN=mean)
dp_by_pos 
#write.csv(dp_by_pos, "./dp_by_pos_iso.csv")   ## CHECK THIS WHEN RERUN

## check qual 
summary(vcf_df$QUAL)
ggplot(vcf_df, aes(QUAL)) +
  geom_density(fill = "lightblue")


summary(vcf_df$AF)
ggplot(vcf_df, aes(x=POS, y=AF, colour=ID)) +
  geom_point()+xlim(0,11000)+ylim(0,1) + geom_hline(yintercept=0.05, linetype="dashed", 
                                                     color = "red", size=1)







## check QUAL
summary(vcf_df$QUAL)

##############################################################################
## check snv greater than 50%

vcf_greater_than_50 = vcf_df[which(vcf_df$FREQ>0.5 & vcf_df$DP >= 1000),]
table(vcf_greater_than_50$POS)

vcf_analysis =vcf_df[which(vcf_df$DP >= 1000 & vcf_df$FREQ>=0.05),]   

#####################################################################
## visualisations
####################################################################

## check SB
colnames(vcf_analysis)
ggplot(vcf_analysis, aes(x=DP, y=SB)) +
  geom_point() + theme(legend.title = element_blank())+geom_hline(yintercept=100, linetype="dashed", 
                                                                  size=0.7)

ggplot(vcf_analysis, aes(x=FREQ, y=SB)) +
  geom_point() + theme(legend.title = element_blank())+geom_hline(yintercept=100, linetype="dashed", 
                                                                  size=0.7)
quantile(vcf_df$SB, probs=0.95)

## check distribution of SNV along POS

ggplot(vcf_analysis, aes(x=POS, y=FREQ)) +
  geom_point()#+xlim(0,30000)+ylim(0,1) + geom_hline(yintercept=0.05, linetype="dashed", 
                                                    #color = "red", size=1)

## find SNV pos frequencies 
check = as.data.frame(table(vcf_analysis$POS))
check = check[order(check[,2],decreasing = TRUE),]
head(check)

ggplot(check, aes(x=check[,1], y=check[,2])) +
  geom_point()


check[which(check$Var1!=10080 & check$Var1!=1012 & check$Var1!=10137& check$Var1!=2866),]

check

write.csv(check,"C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/2nd_round_vcf/check_snv_input_to_regresshaplo.csv")





##################################################
## load csv
###################################################

input_dir = 'C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/passage_vs_non-passage/May_updates/run4_results/'


list_indiv = c(20053,20054,20055,20056,20059,20060,20061,20064,20065,20069)  # 20053 20054 20055 20056 20059 20060 20061 20062 20064 20065 20069 20999 21000 21001 21002 21003 21004 21005 21006 21007 21008
list_indiv
length(list_indiv)

list_indiv = c(20999,21000,21001,21002,21003,21004,21005,21006,21007,21008) 

list_indiv=c(20053,20054,20055,20056,20059,20060,20061,20064,20065,20069,20999,21000,21001,21002,21003,21004,21005,21006,21007,21008)


for( i in 1:20){
  if(i==1){
    check = paste0(input_dir,list_indiv[i],'/final_haplo.csv')
    if (file.exists(check)!='TRUE'){
      next
    }
    df = as.data.frame(read.csv(paste0(input_dir,list_indiv[i],'/final_haplo.csv'),header=FALSE))
    df = df[,1:2]
    df$id = list_indiv[i]
  } else if (i>1){
    check = paste0(input_dir,list_indiv[i],'/final_haplo.csv')
    if (file.exists(check)!='TRUE'){
      next
    }
    df_tmp = as.data.frame(read.csv(paste0(input_dir,list_indiv[i],'/final_haplo.csv'),header=FALSE))
    df_tmp = df_tmp[,1:2]
    df_tmp$id = list_indiv[i]
    df = rbind(df,df_tmp)
  }
}

#dir_snv = df  # view the number of SNVs
#iso_snv = df  # view the numbero f SNVs

table(df$id)   # count the number of haplotypes 

head(df)
nrow(df)

write.csv(df, "./all_snv_to_regress.csv")
