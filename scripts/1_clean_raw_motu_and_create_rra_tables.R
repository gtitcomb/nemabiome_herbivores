
# 1. Clean Raw mOTU Data and Create RRA Table
# 6 October 2021
# Georgia Titcomb

# Use Ctrl+Shift+F10 to refresh whole workspace
##################

library(microDecon)
library(tidyverse)
library(vegan)

# navigate to parent folder before loading this library
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")
library(here) 

##################

# read in data
reads = read.csv(here("data/raw_motu_table_mpala_nematodes.csv"))

# determine clustering level
length(levels(as.factor(reads$cluster_sumatra98)))
elbow_data = data.frame(sum90 = length(levels(as.factor(reads$cluster_sumatra90))),
                        sum95 = length(levels(as.factor(reads$cluster_sumatra95))),
                        sum96 = length(levels(as.factor(reads$cluster_sumatra96))),
                        sum97 = length(levels(as.factor(reads$cluster_sumatra97))),
                        sum98 = length(levels(as.factor(reads$cluster_sumatra98))),
                        sum99 = length(levels(as.factor(reads$cluster_sumatra99))))

clust_elbow = elbow_data %>% 
  pivot_longer(sum90:sum99, names_to="Clustering_Level", values_to="N_mOTUs") %>%
  mutate(clust_level = as.numeric(sub("sum","",Clustering_Level))) %>% 
  ggplot(aes(x=clust_level, y=N_mOTUs))+
  geom_point()+
  geom_line()+
  theme_bw(base_size=14)+
  annotate(geom="text", label=paste(length(reads$id), "Unique Sequences"), x=92.5, y=elbow_data$sum99)+
  scale_x_continuous(breaks=c(90,95,96,97,98,99))+
  labs(x="Clustering Level", y="Number of mOTUs in dataset")
clust_elbow

ggsave(here("plots/1_clust_elbow.png"), dpi=300, height=5, width=5, device="png")


# combine with plate plan data
plate_plan = read.csv(here("data/plate_plan_mpala_nematodes.csv"))
head(plate_plan)

##################################################
#### Digression
# calculate sequence distances for forward and reverse
# taglist = sub(":","",plate_plan$Tag_Combo)
# barcode_distance = data.frame(tags1 = taglist,tags2=taglist, distance=0)
# barcode_distance = barcode_distance %>% pivot_wider(names_from = tags2, values_from=distance)
# barcode_distance = barcode_distance %>% 
#   pivot_longer(-tags1) 
# seqDist(barcode_distance$tags1[3],barcode_distance$name[3])
# 
# library(alakazam)
# for(i in 1:length(barcode_distance$value)){
#   barcode_distance$value[i] = seqDist(barcode_distance$tags1[i], barcode_distance$name[i])
# }
# 
# ggplot(barcode_distance, aes(x=tags1, y=name))+
#   geom_tile(aes(fill=value))
# 
# 
# lower_distance = barcode_distance %>% 
#   filter(value<5) %>% 
#   filter(value>0)
# 
# plate_plan$Tag_Combo2 = sub(":","",plate_plan$Tag_Combo)
# lower_distance$sample1 = plate_plan[match(lower_distance$tags1, plate_plan$Tag_Combo2),]$Sample
# lower_distance$sample2 = plate_plan[match(lower_distance$name, plate_plan$Tag_Combo2),]$Sample
# 
# blank_ld = lower_distance[grep("Blank", lower_distance$sample1),] 
# levels(as.factor(blank_ld$sample1))
######################################################

# transpose samples
names(reads) = sub("sample.","",names(reads))

# any mismatches
names(reads)[which((names(reads) %in% plate_plan$Sample2)==F)]
plate_plan$Sample2[which((plate_plan$Sample %in% names(reads))==F)]

# examine distribution of total reads
read_dist = as.data.frame(colSums(reads[,28:411]))
read_dist$Sample = row.names(read_dist)
names(read_dist)[1]="Reads"

# join to plate plan
read_dist = left_join(read_dist, plate_plan)
read_dist %>% 
  filter(is.na(Plate))

# assign sample type information
read_dist$Sample_type = ""
read_dist$Sample_type[grep("Blank",read_dist$Sample)]="Blank"
read_dist$Sample_type[grep("POS",read_dist$Sample)]="Positive"
read_dist$Sample_type[grep("NEG",read_dist$Sample)]="Negative"
read_dist$Sample_type[which(read_dist$Sample_type=="")]="Sample"

# quick boxplot
ggplot(read_dist, aes(x=Sample_type, y=Reads))+
  geom_boxplot()+
  facet_wrap(~Plate)+
  scale_y_log10()

# there is one blank in plate one that has caused problems
# exclude that here
read_info  = read_dist %>% 
  filter(Sample!="Blank02") %>% 
  group_by(Sample_type) %>% 
  summarize_at(vars(Reads), funs(mean, median, min, max, n()))

write.table(read_info, here("docs/1_read_info.txt"), sep="\t", row.names=F)

#reads = reads[,-(which(names(reads)=="Blank02"))]

########################################################
# Filtering steps:
# 1. All Samples with <1000 reads
# read_dist_1 = read_dist %>% 
#   filter(Reads >= 1000)
# read_dist_1 %>% group_by(Sample_type) %>% summarize(n())
# read_dist %>% group_by(Sample_type) %>% summarize(n())
# # 37 samples had fewer than 1000 reads


# 1. Cluster
# Calculate the sum of reads for each cluster
clust_reads = reads %>% 
  group_by(cluster_sumatra98) %>% 
  summarize_at(vars(Blank01:MRC_WAT_205), funs(sum))

# step1 read tracking
a_n_motu = clust_reads %>% group_by(cluster_sumatra98) %>% nrow()
a_n_samp = length(grep("MRC", names(clust_reads)))

#organizing vector
col_order = plate_plan %>% 
  arrange(Species) %>% 
  dplyr::select(Sample) #%>% 
  #filter(Sample != "Blank02")
col_order

# put negs with blanks
grep("NEG",col_order$Sample)
col_order = col_order[c(1:47,324,325,326, 48:323,327:384),]

clust_reads = clust_reads[,c(1,match(col_order,names(clust_reads)))]


# 2. Use decontaminator
read_to_decontam = clust_reads %>% 
  rename(OTU_ID = cluster_sumatra98)
read_to_decontam$OTU_ID = paste("mOTU",read_to_decontam$OTU_ID, sep="_")

# get groups
sp_counts = plate_plan %>% 
  arrange(Species) %>% 
  filter(Species != "Blank") %>% 
  filter(Species != "Negative") %>% 
  group_by(Species) %>% 
  summarize(n=n())

# total number of samples (includes 325 samples and 8 positives)
sum(sp_counts$n)

# run microDecon
set.seed(0); decontaminated = decon(data = read_to_decontam, numb.blanks = 51, numb.ind =333, taxa=F, runs=2)

# average read depth for samples, after subtractions
mean(colSums(decontaminated$decon.table[,-c(1:2)]))

b_n_motu = a_n_motu - dim(decontaminated$OTUs.removed)[1]
b_n_samp = length(grep("MRC", names(decontaminated$decon.table)))

# transpose
clust_reads2 = decontaminated$decon.table %>% 
  pivot_longer(cols=MRC_BUF_105:MRC_WAR_103, names_to="Sample", values_to="Reads") %>% 
  pivot_wider(names_from=OTU_ID, values_from=Reads, -Mean.blank)
dim(clust_reads2) # ensure 333 rows

# remove empty columns and unassigned sequence motu
clust_reads2 = clust_reads2[,-(which(colSums(clust_reads2[,-1])==0)+1)]
clust_reads2 = clust_reads2[,-(which(names(clust_reads2)=="mOTU_NA"))]
which(colSums(clust_reads2[,-1])==0) # check

# motus lost from positive sample and NA removed
c_n_motu = b_n_motu - 2
c_n_samp = b_n_samp - length(which(rowSums(clust_reads2[,-1])==0))

length(grep("POS",clust_reads2$Sample))
# 333 - 8
length(grep("mOTU", names(clust_reads2)))


# calculate reads again
clust_reads2$reads = rowSums(clust_reads2[,-1])

# exclude <1000
clust_reads3 = clust_reads2 %>% 
  filter(reads>=1000)

d_n_motu = c_n_motu - length(which(colSums(clust_reads3[,-1])==0))
d_n_samp = length(grep("MRC", clust_reads3$Sample))


###### Calculate mOTU loss as a function of RRA threshold ######

rra_sensitivity = function(data, startcol, endcol, thresholds){
  # data = data frame
  # startcol = name of first column
  # endcol = name of end column
  # thresholds = vector of thresholds to try
  
  store = data.frame(thresh = rep(0,length(thresholds)), n_mOTUs = rep(0, length(thresholds)))
  for(i in 1:length(thresholds)){
    mutated = data %>% 
      mutate_at(vars(startcol:endcol), funs(./reads)) %>% 
      mutate_at(vars(startcol:endcol), funs(ifelse(. < thresholds[i], 0,.))) %>% 
      mutate_at(vars(startcol:endcol), funs(ifelse(.>0,1,0)))
    new_data = dplyr::select(mutated, startcol:endcol) * dplyr::select(data, startcol:endcol)
    
    
    store$n_mOTUs[i] = (dim(mutated)[2]-2) - length(which(colSums(new_data)==0))
    store$thresh[i] = thresholds[i]
  }
  return(store)
}

sens = rra_sensitivity(clust_reads3, "mOTU_1", "mOTU_568", thresholds = seq(from = 0, to=0.03, by=0.001))

rra_elbow = ggplot(sens, aes(x=thresh, y=n_mOTUs))+
  geom_point()+
  geom_path()+
  theme_bw(base_size=14)+
  labs(x="RRA Threshold", y="Number of mOTUs in dataset")+
  geom_vline(xintercept = c(0.002, 0.02), linetype="dotted", col="red", size=1)
rra_elbow
ggsave(here("plots/1_rra_elbow.png"), dpi=300, height=5, width=5, device="png")


#######
# remove RRA lower threshold

# drop mOTUs representing less than .2% of reads in each sample
decision_include = clust_reads3 %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(./reads)) %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(ifelse(. < 0.002, 0,.))) %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(ifelse(.>0,1,0)))

clust_reads_n1 = dplyr::select(decision_include, mOTU_1:mOTU_568) * dplyr::select(clust_reads3, mOTU_1:mOTU_568)

# calculate n mOTUs excluded
e_n_motu1 = d_n_motu - length(which(colSums(clust_reads_n1)==0))
e_n_samp1 = d_n_samp - length(which(rowSums(clust_reads_n1) ==0))

clust_reads_n1 = clust_reads_n1[,-(which(colSums(clust_reads_n1)==0))]

# add sample back and calculate new reads
clust_reads_n1$reads = rowSums(clust_reads_n1)
clust_reads_n1$Sample = clust_reads3$Sample


# drop mOTUs representing less than 2% of reads in each sample
decision_include = clust_reads3 %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(./reads)) %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(ifelse(. < 0.02, 0,.))) %>% 
  mutate_at(vars(mOTU_1:mOTU_568), funs(ifelse(.>0,1,0)))

clust_reads_n2 = dplyr::select(decision_include, mOTU_1:mOTU_568) * dplyr::select(clust_reads3, mOTU_1:mOTU_568)

# calculate n mOTUs excluded
e_n_motu2 = d_n_motu - length(which(colSums(clust_reads_n2)==0))
e_n_samp2 = d_n_samp - length(which(rowSums(clust_reads_n2) ==0))

# remove empty columns
clust_reads_n2 = clust_reads_n2[,-(which(colSums(clust_reads_n2)==0))]
dim(clust_reads_n2)

# add sample back and calculate new reads
clust_reads_n2$reads = rowSums(clust_reads_n2)
clust_reads_n2$Sample = clust_reads3$Sample


# rarefy to 1000 reads
# Used repeated rarefaction:
# clust_n1
min(clust_reads_n1$reads)
set.seed(0); rare_stack = replicate(100, rrarefy(clust_reads_n1[,1:(dim(clust_reads_n1)[2]-2)], min(clust_reads_n1$reads)))
# average through stack
mean_rare1 = as.data.frame(apply(rare_stack, 1:2, mean))

# calculate RRA by dividing by minimum
clust_reads_n1 = mean_rare1 %>% 
  mutate_at(vars(mOTU_1:mOTU_488), funs(./min(clust_reads_n1$reads)))
clust_reads_n1$Sample = clust_reads3$Sample

# table 1 is our lower bound
table_1 = left_join(clust_reads_n1, read_dist, by="Sample")
head(table_1)
dim(table_1)
end1 = dim(table_1)[2]-14


# clust_n2
min(clust_reads_n2$reads)
set.seed(0); rare_stack = replicate(100, rrarefy(clust_reads_n2[,1:(dim(clust_reads_n2)[2]-2)], min(clust_reads_n2$reads)))
# average through stack
mean_rare2 = as.data.frame(apply(rare_stack, 1:2, mean))

# calculate RRA by dividing by minimum
clust_reads_n2 = mean_rare2 %>% 
  mutate_at(vars(mOTU_1:mOTU_94), funs(./min(clust_reads_n2$reads)))
clust_reads_n2$Sample = clust_reads3$Sample

# table 1 is our lower bound
table_2 = left_join(clust_reads_n2, read_dist, by="Sample")
head(table_2)
dim(table_2)
end2 = dim(table_2)[2]-14


#quick mds visual
nmds1 = metaMDS(table_1[,1:end1])
plot(nmds1)
levels(as.factor(table_1$Species))
mds_data1 = data.frame(x1 = nmds1$points[,1], x2=nmds1$points[,2], sample_type=table_1$Sample_type, species=table_1$Species, sample=table_1$Sample, reads=table_1$Reads)
colors = c( "violet", "orange", "violet","blue","magenta", "blue2", "green", "orange", "blue", "magenta",
           "blue", "green", "magenta", "blue", "blue", "blue", "magenta", "red", "green", "magenta")
ggplot(mds_data1, aes(x=x1, y=x2))+
  geom_point(aes(col=species), alpha=0.2)+
  ggConvexHull::geom_convexhull(aes(fill=species), alpha=0.2)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  theme_bw()

# 
nmds2 = metaMDS(table_2[,1:end2])
plot(nmds2)
mds_data2 = data.frame(x1 = nmds2$points[,1], x2=nmds2$points[,2], sample_type=table_2$Sample_type, species=table_2$Species, sample=table_2$Sample, reads=table_2$Reads)
ggplot(mds_data2, aes(x=x1, y=x2))+
  geom_point(aes(col=species), alpha=0.2)+
  ggConvexHull::geom_convexhull(aes(fill=species), alpha=0.2)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  theme_bw()

#table_2[which(mds_data2$x1 < -5),]

# exclude remaining controls
table_1 = table_1 %>% 
  filter(Sample_type=="Sample")
table_2 = table_2 %>% 
  filter(Sample_type=="Sample")

# exclude spurious hippo, waterbuck, and hybrid zebra
table_1 = table_1 %>% 
  filter(Species != "Hybrid Zebra") %>% 
  filter(Species != "Hybrid zebra") %>% 
  filter(Sample != "MRC_17_HIP_93")
table_2 = table_2 %>%
  filter(Species != "Hybrid Zebra") %>% 
  filter(Species != "Hybrid zebra") %>% 
  filter(Sample != "MRC_17_HIP_93")


# double check empty columns
(which(colSums(table_1[,1:end1])==0))
#table_1 = table_1[,-(which(colSums(table_1[,1:end1])==0))]
dim(table_1)

f_n_motu1 = dim(table_1[,1:end1])[2]
f_n_samp1 = dim(table_1[,1:end1])[1]

table_1 = table_1 %>% 
  rename(Filtered_Reads = Reads)

# double check empty columns
(which(colSums(table_2[,1:end2])==0)) # this should be none

dim(table_2)

f_n_motu2 = dim(table_2[,1:end2])[2]
f_n_samp2 = dim(table_2[,1:end2])[1]

table_2 = table_2 %>% 
  rename(Filtered_Reads = Reads)

# there are 15 metadata columns, so this should be 110 and 94 mOTUs respectively

##############################
# save the cleaned data tables
write.csv(dplyr::select(table_1, mOTU_1:Filtered_Reads), here("data/RRA_table_1.csv"), row.names=F)
write.csv(dplyr::select(table_2, mOTU_1:Filtered_Reads), here("data/RRA_table_2.csv"), row.names=F)

# sample and otu filtering summary
filt_summary = data.frame(step = c("initial","decontam","less_1000","POS_NA_mOTU","RRA","extras"),
           n_motu1 = c(a_n_motu, b_n_motu, c_n_motu, d_n_motu, e_n_motu1, f_n_motu1),
           n_motu2 = c(a_n_motu, b_n_motu, c_n_motu, d_n_motu, e_n_motu2, f_n_motu2),
           n_samp1 = c(a_n_samp, b_n_samp, c_n_samp, d_n_samp, e_n_samp1, f_n_samp1),
           n_samp2 = c(a_n_samp, b_n_samp, c_n_samp, d_n_samp, e_n_samp2, f_n_samp2))

write.table(filt_summary, here("docs/1_filt_summary.txt"), sep="\t", row.names=F)  
