
cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq

#!/bin/bash
for i in *_R1.fastq
do
 
# trim by min. quality 30 and min. length 50
time cutadapt -q 30,30 --minimum-length 50   \
         -o $(basename ${i} _R1.fastq)_R1_trim1.fastq \
         -p $(basename ${i} _R1.fastq)_R2_trim1.fastq\
         $(basename ${i} _R1.fastq)_R1.fastq \
         $(basename ${i} _R1.fastq)_R2.fastq \
         >> cutadapt_trim1.txt
done 


#!/bin/bash
for i in *_R1_trim1.fastq
do
 
# trim 20 bases from each ends of the read; R1
cutadapt -u 5 -u -5  \
         -o $(basename ${i} _R1_trim1.fastq)_R1_trim2.fastq \
         $(basename ${i} _R1_trim1.fastq)_R1_trim1.fastq \
         >> cutadapt_trim2_R1.txt

# trim 20 bases from each ends of the read; R2


cutadapt -u 5 -u -5  \
         -o $(basename ${i} _R1_trim1.fastq)_R2_trim2.fastq \
         $(basename ${i} _R1_trim1.fastq)_R2_trim1.fastq \
         >> cutadapt_trim2_R2.txt

done


## 1st round of alignment 

cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq
#!/bin/bash
for i in *_R1_trim1.fastq
do
## Create SAM file
# code formate: bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sa
   
   ref=Den1_HM469967_Thailand_2007.fasta         #KM279601.fasta
   fq1=$(basename ${i} _R1_trim1.fastq)_R1_trim2.fastq \
   fq2=$(basename ${i} _R1_trim1.fastq)_R2_trim2.fastq \

   echo "ref is $ref, fastq files are: $fq1 and $fq2" \ 

   # indx the fasta file
  bwa mem -t 8 $ref $fq1 $fq2 > $(basename ${i} _R1_trim1.fastq).sam     #takes time


## Create BAM file
samtools view -Sb -h $(basename ${i} _R1_trim1.fastq).sam  > $(basename ${i} _R1_trim1.fastq).bam

## Sort BAM files by coordidates
samtools sort $(basename ${i} _R1_trim1.fastq).bam -o $(basename ${i} _R1_trim1.fastq)_sorted.bam

## Index sorted BAM
samtools index $(basename ${i} _R1_trim1.fastq)_sorted.bam
 
done 

#######################################################################################################
## total no. of reads, coverage, and mapping quality (MAPQ) check
#######################################################################################################
## total reads
#!/bin/bash
for i in *_sorted.bam
do 
#echo $i
samtools view -c $i #>> total_read_of_bam.txt
done

## check length of the samples
#!/bin/bash
for i in *sorted.bam
do 
samtools view -H $i | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'
done

## calculating average read depth including only the covered reads
#!/bin/bash
for i in *sorted.bam
do 
#echo $i
samtools depth $i | awk '{sum+=$3} END { print sum/10735}' #>> average_coverage.txt
# 10176 = total number of bases in the sample; if "NR", then only the covered reads is taken into account
done 

## mapping quality (check the number of reads)
#!/bin/bash
for i in *_sorted.bam
do 
samtools view -q30 -c $i
done 


#######################################################################################################
## filter out mapq <20
#######################################################################################################
cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq
#!/bin/bash
for i in *_sorted.bam
do 
samtools view -b -q 20 $i > $(basename ${i} _sorted.bam)_sorted_mapq.bam
samtools index $(basename ${i} _sorted.bam)_sorted_mapq.bam
done 

#######################################################################################################
## remove singletons 
#######################################################################################################
#!/bin/bash
for i in *_sorted_mapq.bam
do 
samtools view -b -F 0x8 $i > $(basename ${i} .bam)_nosingle.bam
samtools index $(basename ${i} .bam)_nosingle.bam

## create uniquely mapped reads bam 
#######################################################################################################
samtools view -b -F 4 $(basename ${i} .bam)_nosingle.bam > $(basename ${i} .bam)_nosingle_unique.bam

samtools index $(basename ${i} .bam)_nosingle_unique.bam

done 

######################################################################################################
## deduplicate after removing singletons and low quality. - DO NOT APPLY THIS STEP TO THE NMC DATA
#######################################################################################################
#!/bin/bash
for i in *_unique.bam
do 
## Deduplicate 
java -jar /home/yrj21/rds/hpc-work/privatemodules/picard/build/libs/picard-2.25.0-5-ga2f44ae-SNAPSHOT-all.jar MarkDuplicates \
           -I $i \
           -O $(basename ${i} .bam)_dedup.bam \
           -M marked_dup_metrics.txt  \
           --REMOVE_DUPLICATES true
done 


## check final bam files 
#################################################
#!/bin/bash
for i in *_dedup.bam
do 
#echo $i
samtools view -c $i #>> total_read_of_bam.txt
done

#!/bin/bash
for i in *_dedup.bam
do 
#echo $i
samtools depth $i | awk '{sum+=$3} END { print sum/10735}' #>> average_coverage.txt
# 10176 = total number of bases in the sample; if "NR", then only the covered reads is taken into account
done 


######################################################################################################
## Variant calling and filtering for consensus fasta
#######################################################################################################
cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq
#!/bin/bash
for i in *_unique.bam
do 
## Sort BAM files by coordidates
samtools sort $i -o $(basename ${i} .bam)_sorted.bam

## Index sorted BAM
samtools index $(basename ${i} .bam)_sorted.bam
## variant call using lofreq 
/home/yrj21/rds/hpc-work/privatemodules/lofreq_star-2.1.2/bin/lofreq call --verbose -f Den1_HM469967_Thailand_2007.fasta   \
-o $(basename ${i} .bam).vcf $(basename ${i} .bam)_sorted.bam  

done 

##########################################################################################################
## index vcf for filtering and align with consensus
##########################################################################################################
#!/bin/bash
for i in *_unique.vcf
do 
## index trimmed vcf
/home/yrj21/rds/hpc-work/privatemodules/tabix-0.2.6/bgzip -c $(basename ${i} .vcf).vcf > $(basename ${i} .vcf).vcf.gz

/home/yrj21/rds/hpc-work/privatemodules/tabix-0.2.6/tabix -p vcf $(basename ${i} .vcf).vcf.gz

## filter vcf
/home/yrj21/rds/hpc-work/privatemodules/bcftools-1.9/bcftools filter --exclude 'DP<1000 || (DP4[2]+DP4[3])/DP<0.5' $(basename ${i} .vcf).vcf.gz \
 > $(basename ${i} .vcf)_DP250.vcf

 ## construct individual consensus seuquence
###########################################
## index trimmed vcf
/home/yrj21/rds/hpc-work/privatemodules/tabix-0.2.6/bgzip -c $(basename ${i} .vcf)_DP250.vcf > $(basename ${i} .vcf)_DP250.vcf.gz

/home/yrj21/rds/hpc-work/privatemodules/tabix-0.2.6/tabix -p vcf $(basename ${i} .vcf)_DP250.vcf.gz

## consensus
/home/yrj21/rds/hpc-work/privatemodules/bcftools-1.9/bcftools consensus \
          -f Den1_HM469967_Thailand_2007.fasta $(basename ${i} .vcf)_DP250.vcf.gz  \  # check reference
          -o $(basename ${i} .vcf)_ind.fasta
done 




 ########################################################################################################
 ## convert downsample by 0.5 bam files to fastq files
 ########################################################################################################
cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq
#!/bin/bash
for i in *_unique.bam
do 
 
# sort paired read alignment .bam file (sort by name -n)
samtools sort -n $i -o $(basename ${i} .bam)_sorted.bam

# save fastq reads in separate R1 and R2 files
time samtools fastq -@ 8 \
    -1 ./align_consensus_fastq/$(basename ${i} .bam)_R1.fastq.gz \
    -2 ./align_consensus_fastq/$(basename ${i} .bam)_R2.fastq.gz \
    -0 ./align_consensus_fastq/null -n -s ./align_consensus_fastq/singleton \
    $(basename ${i} .bam)_sorted.bam
done 


#######################################################################################################
## align with the individual consensus
#######################################################################################################
cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq/align_consensus_fastq

#!/bin/bash
for i in *_R1.fastq.gz
do
## Create SAM file
# code formate: bwa mem ref_genome.fasta input_file_R1.fastq input_file_R2.fastq > output.sa
   
   ref=../$(basename ${i} _R1.fastq.gz)_ind.fasta    #REMEMBER TO CHECK HERE

   fq1=$(basename ${i} _R1.fastq.gz)_R1.fastq.gz\
   fq2=$(basename ${i} _R1.fastq.gz)_R2.fastq.gz \

   echo "ref is $ref, fastq files are: $fq1 and $fq2" \ 

   # indx the fasta file
    bwa index $ref 
   
time   bwa mem -t 10 $ref $fq1 $fq2 > ../$(basename ${i} _R1.fastq.gz)_ind.sam     #takes time

## sam to bam
## Create BAM file
samtools view -Sb ../$(basename ${i} _R1.fastq.gz)_ind.sam  > ../$(basename ${i} _R1.fastq.gz)_ind.bam 

## Sort BAM files by coordidates
samtools sort ../$(basename ${i} _R1.fastq.gz)_ind.bam -o ../$(basename ${i} _R1.fastq.gz)_ind_sorted.bam 

## Index sorted BAM
samtools index ../$(basename ${i} _R1.fastq.gz)_ind_sorted.bam 
done 

## 2nd variant calling 
cd /home/yrj21/rds/rds-hs743-arbodynamic/Raina/cluster_june/fastq
#!/bin/bash
for i in *_ind_sorted.bam
do
samtools faidx $(basename ${i} _ind_sorted.bam)_ind.fasta 
/home/yrj21/rds/hpc-work/privatemodules/lofreq_star-2.1.2/bin/lofreq call --verbose -f $(basename ${i} _ind_sorted.bam)_ind.fasta  \
-o /$(basename ${i} _ind_sorted.bam)_ind.vcf $i
done 

