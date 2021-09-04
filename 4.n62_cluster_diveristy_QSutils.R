#BiocManager::install("QSutils")
library(QSutils)
library(seqinr)
library(doParallel)



###########################################
## edit the name of fasta files 
###########################################
## sample IDs
list_indiv = c(1207,1223,1238,1276,1278,1280,1282,1284,1294,1307,1310,1312,1315,1324,1325,1338,1364,1378,1388,2430,2432,2433,2434,2436,2442,2443,2444,2447,2451,2452,2453,2454,2455,2456,2457,2460,2462,2463,2464,2467,2468,2469,2470,2471,2472,2477,2478,2480,2481,2482,2483,2909,2910,2911,2912,2913,2914,2929,2930,2931,2935,3023) # 20053 20054 20055 20056 20059 20060 20061 20062 2006$length(list_indiv)
length(list_indiv)


for (i in 1:length(list_indiv)){

 input_dir = 'C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/regresshaplo_out_snv03/' # CHECK snv03 or ssnv03

 seq = getSequence(read.fasta(paste0(input_dir,"/",list_indiv[i],"/",list_indiv[i],"_hap_final.fasta")))

final_haplo = paste0(input_dir,"/",list_indiv[i],"/",'final_haplo.csv')

if(file.exists(final_haplo)=='TRUE'){
 variants = read.csv(paste0(input_dir,"/",list_indiv[i],"/",'final_haplo.csv'), header = F, colClasses = 'character')
 
 names = paste0(list_indiv[i], '.haplo',1:length(variants[,1]), '|10735|', as.character(round(as.numeric(as.character(variants[,1])), digits = 2)))

 setwd(input_dir)
 seqinr::write.fasta(seq, names = names, file.out =paste0(list_indiv[i], '_qstuils.fasta')) 
 
}else{
 name = paste0(list_indiv[i], '.haplo1|10735|1') 
 setwd(input_dir)
 seqinr::write.fasta(seq, names = names, file.out = paste0(list_indiv[i], '_qstuils.fasta'))
}
}


#########################################################
## entropy calculation 
#########################################################

setwd('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/regresshaplo_out_snv03') # CHE

m = matrix("", nrow =length(list_indiv), ncol = 6)
rownames(m) = list_indiv
colnames(m) = c("no_of_haplotypes","no_of_mutation","Shannon_entropy","Shannon_entropy_variance","mutation_frequency","nucleotide_diveristy")

for(i in 1:length(list_indiv)){
lst <- ReadAmplSeqs(paste0('C:/Users/micro/OneDrive - University of Cambridge/Documents/Mphil_project/June_cluster/regresshaplo_out_snv03/', list_indiv[i],'_qstuils.fasta'),type="DNA")

no_hap = length(lst$hseqs)
no_mut = TotalMutations(lst$hseqs)
#Shannon entropy
shannon = Shannon(lst$nr)
shannon_var = ShannonVar(lst$nr)
# compute p-distance for mutation freqeuncy and nucleotide diveristy estimation
dst = DNA.dist(lst$hseqs, model = "raw", gamma = FALSE, pairwise.deletion = FALSE)
# haplotype-based mutation freuqency
mutation_freq=MutationFreq(dst) 
# haplotype-based nucleotide diversity 
nucleotide_diversity = NucleotideDiversity(dst)

m[i,1]=no_hap
m[i,2]=no_mut
m[i,3]=shannon
m[i,4]=shannon_var
m[i,5]=mutation_freq
m[i,6]=nucleotide_diversity 
}

m
data = data.frame(m)
data$ID <- rownames(data)
colnames(data)[1]


write.csv(file = "diveristy_estimation_table.csv", x=data)

##comparison between direct and passage 
dir_list = data[1:10,]
iso_list=data[11:20,]

dir_list = dir_list[order(dir_list$pair_id,decreasing=FALSE),]
iso_list = iso_list[order(iso_list$pair_id,decreasing=FALSE),]

colnames(dir_list)

## no. of haplotypes  
x= dir_list[,'no_of_haplotypes']
y= iso_list[,'no_of_haplotypes']

x= as.numeric(x)
y= as.numeric(y)

#wilcoxon test
wilcox.test(x, y, paired = TRUE, alternative = "two.sided")

# F-test
var.test(x, y, alternative = "two.sided")

## Shannon entropy
x= dir_list[,'Shannon_entropy']
y= iso_list[,'Shannon_entropy']

x= as.numeric(x)
y= as.numeric(y)

wilcox.test(x, y, paired = TRUE, alternative = "two.sided")

# F-test
var.test(x, y, alternative = "two.sided")


## nucleotide diveristy
x= dir_list[,'nucleotide_diveristy']
y= iso_list[,'nucleotide_diveristy']

x= as.numeric(x)
y= as.numeric(y)

wilcox.test(x, y, paired = TRUE, alternative = "two.sided")

## Mutation freuqnecy by haplotype
x= dir_list[,'mutation_frequency']
y= iso_list[,'mutation_frequency']

x= as.numeric(x)
y= as.numeric(y)

wilcox.test(x, y, paired = TRUE, alternative = "two.sided")


## SNV variability (use input SNV)
#snv003
x = c(8,1,9,30,3,5,5,3,11)  # take out the outlier 
y = c(8,0,9,0,2,2,7,3,7)

#snv005
x = c(5,0,9,4,0,1,3,7,4,2,5)
y = c(3,0,6,5,0,2,2,2,5,3,3)



x= as.numeric(x)
y= as.numeric(y)

var.test(x, y, alternative = "two.sided")

wilcox.test(x, y, paired = TRUE, alternative = "two.sided")
