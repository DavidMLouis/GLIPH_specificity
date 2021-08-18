#!/usr/bin/env Rscript

# Author: David M. Louis
# Contact: dmlouis@stanford.edui
# Contact: dmlouis87@gmail.com

############################### Libraries ###############################

library(ggplot2)
library(plyr)

############################### Arguments ###############################

# expected format of input file is
cmd_args = commandArgs(trailingOnly = TRUE);
# midread_table_file = cmd_args[6] ???????????

input_file=read.csv(cmd_args, header=T)
# tableLength=length(midread_table[,1])
# print(cmd_args) #prints file name
print(paste0("input length: ", length(rownames(input_file))))
#midread_matrix=data.matrix(midread_table)

############################################################################################
# Analysis Functions
############################################################################################

#clean up GLIPH output
GLIPH_cleaner <- function(GLIPH_df) {
  GLIPH_df[GLIPH_df==""] <- NA
  data <- GLIPH_df[rowSums(is.na(GLIPH_df)) != ncol(GLIPH_df), ]
  data <- data[!grepl("single", data$pattern), ]
  if (sum(is.na(data$Sample)) > 0) {
    print("Samples containes NA. Removing all rows. Add missing values if row data is desired in outcome")
    data <- data[!is.na(data$Sample),]
  }
  return(data)
}

#split data
data_splitter <- function(df_input) {
  split_data <- df_input$Sample
  # split_data <- gsub("_heart", "", split_data)
  split_data <- gsub("^.*_", "", split_data)
  return(split_data)
}

#creating specificity sheets
#inputs:dataframe #2nd counted column #2nd columns merge files (needs 2)
#ex: specificty_table(df, "split_column", "split_col_item_1", "split_col_item_2")
specificity_table <- function(clean_input, count_col, count_data_1, count_data_2){
  #creating specificity tables
  data_count <- count(clean_input, vars = c("pattern", count_col))
  spec_table <- merge(subset(data_count, disease == count_data_1), subset(data_count, disease == count_data_2), by = 'pattern', all = T)
  colnames(spec_table)[colnames(spec_table) == 'freq.x'] <- count_data_1
  colnames(spec_table)[colnames(spec_table) == 'freq.y'] <- count_data_2
  spec_table[[paste0(count_col,'.x')]] <- NULL
  spec_table[[paste0(count_col,'.y')]] <- NULL
  spec_table[is.na(spec_table)] <- 0
  #advanced frequency
  spec_table[[paste0(count_data_1, "Freq")]] <- spec_table[[count_data_1]]/sum(spec_table[[count_data_1]])*100
  spec_table[[paste0(count_data_2, "Freq")]] <- spec_table[[count_data_2]]/sum(spec_table[[count_data_2]])*100
  #creating advanced specificity
  spec_table[[paste0(count_data_1, "Spec")]] <- spec_table[[count_data_1]]/(rowSums(spec_table[,c(count_data_1, count_data_2)]))*100
  spec_table[[paste0(count_data_2, "Spec")]] <- spec_table[[count_data_2]]/(rowSums(spec_table[,c(count_data_1, count_data_2)]))*100
  spec_table[['sum']] <- rowSums(spec_table[,c(count_data_1, count_data_2)])
  return(spec_table)
}


# test <- Reduce(function(x,y) merge(x,y,by="pattern",all=TRUE),
#                list(subset(stem_01_data, cell_type == "CD4"), 
#                     subset(stem_01_data, cell_type == "nonTR1"), 
#                     subset(stem_01_data, cell_type == "Tallo"), 
#                     subset(stem_01_data, cell_type == "Tallo10"), 
#                     subset(stem_01_data, cell_type == "TR1")))





spec_Fig <- function(spec_df){
  mid=50
  spec_plot <- ggplot(environment = environment()) +
    geom_point(data = spec_df, aes(
      x = spec_df[,4],
      y = spec_df[,5],
      fill = spec_df[,6],
      size=spec_df[['sum']]),
      pch=21, alpha=1, position = position_jitter(w = max(c(spec_df[,4], spec_df[,5]) + sd(c(spec_df[,4], spec_df[,5])))*.02, h = max(c(spec_df[,4], spec_df[,5]) + sd(c(spec_df[,4], spec_df[,5])))*.02)) +
    scale_size_continuous(range = c(1, 15), guide = "none") +
    scale_fill_gradient2(midpoint=mid, low="blue", mid="white",high="red", name = colnames(spec_df[6])) +
    # geom_point(data=samp ,aes(MSfreq, Healthyfreq, fill=MSspec), pch=21, size=samp$sum, alpha=1, position = position_jitter(w = 0.002, h = 0.002)) +
    # scale_color_gradient(low = 'green', high = 'black') +
    xlim(NA,max(c(spec_df[,4], spec_df[,5]) + sd(c(spec_df[,4], spec_df[,5])))) +
    ylim(NA,max(c(spec_df[,4], spec_df[,5]) + sd(c(spec_df[,4], spec_df[,5])))) +
    # ggtitle("MS CD4") +
    xlab(colnames(spec_df[4])) +
    ylab(colnames(spec_df[5])) +
    theme(panel.background = element_blank(), axis.line = element_line(color = 'black'),panel.border = element_rect(colour = "black", fill=NA, size=2),
          axis.text = element_text(size = 20))
  return(spec_plot)
}

#Fstat Closer to 1 the more equal the data is, closer to 0 the farther away they are
#differences within group is larger than difference between groups 
correlation_figure <- function(spec_table) {
  # data correlations with fstat
  plot(spec_table[,2] ~ spec_table[,3], main="Data Correlation",  xlab=colnames(spec_table[2]), ylab = colnames(spec_table[3])) #shows randomness of data (maybe make it jitter)
  #getting fstatistic(?) of observed data
  # summary(lm(data_specificity[,2] ~ data_specificity[,3]))$fstat[1]
  mtext(paste0("Fstat: ", summary(lm(spec_table[,2] ~ spec_table[,3]))$fstat[1]))
  corr_plot <- recordPlot()
  return(corr_plot)
}

permutations <- function(spec_data){
  simulation_counts <- 1000
  spec_data$rand_add <- NA
  spec_data$rand_disease_spec <- NA
  nreps <- simulation_counts
  disease_simulations <- numeric(nreps)
  for(i in 1:nreps){
    disease_simulations[i] <- mean(
      ifelse(
        sample(
          c(colnames(spec_data[2]), colnames(spec_data[3])),
          nrow(spec_data), replace = T
        ) == colnames(spec_data)[2], spec_data[,6], spec_data[,7]
      )
    )
  }
  return(disease_simulations)
}

# permutation_plot(data_specificity, permutation_simulation)
permutation_plot <- function(OG_file, perm_data) {
  simulation_counts <- 1000
  observe1 <- mean(OG_file[,6]) #input[2]
  observe2 <- mean(OG_file[,7]) #input[3]
  plot(density(perm_data), main = "Permutation Distribution of Specificity", 
       # xlab = NA,
       xlab = paste0("Simulated p-value = ",
                     if (mean(OG_file[,6]) > mean(OG_file[,7])) {
                       sum(perm_data > mean(OG_file[,6]))/simulation_counts
                     } else {
                       sum(perm_data > mean(OG_file[,7]))/simulation_counts
                     }),
       xlim=c(
         min(perm_data, observe1, observe2), 
         max(perm_data, observe1, observe2)))
  abline(v=observe1)
  abline(v=observe2)
  text(observe1, max(density(perm_data)$y)/2, colnames(OG_file[2]), srt=90, pos = 3) #srt=90, pos = 3
  text(observe2, max(density(perm_data)$y)/2, colnames(OG_file[3]), srt=90, pos = 3) #srt=90, pos = 3
  perm_plot <- recordPlot()
  return(perm_plot)
}

#Function ends

#pvalue
# pf(summary(lm(data_specificity[,4] ~ data_specificity[,5]))$fstat[1], summary(lm(data_specificity[,4] ~ data_specificity[,5]))$fstat[2], summary(lm(data_specificity[,4] ~ data_specificity[,5]))$fstat[3], lower.tail = F)

############################### Analysis  #################################

clean_input <- GLIPH_cleaner(input_file)
clean_input$disease <- data_splitter(clean_input)
data_specificity  <- specificity_table(clean_input, "disease", unique(clean_input$disease)[1], unique(clean_input$disease)[2])
bubble_specificity <- spec_Fig(data_specificity)
# correlation_plot <- correlation_figure(data_specificity) #moved to outputs
permutation_simulation <- permutations(data_specificity)

############################### Outputs #################################
#creating directory
mainDir <- "specificity_analysis"
dir.create(file.path(mainDir), showWarnings = FALSE)
setwd(file.path(mainDir))

############## exporting data and figures
# stat table
write.csv(data_specificity, file = paste0("specificity_table_", gsub("^.*/", "", cmd_args)), quote = F, row.names = F)

#bubble plot
png("specificity_plot.png",height=600,width=800)
bubble_specificity
dev.off()

#correlation plot
png("correlation_plot.png",height=600,width=800)
# correlation_plot
correlation_figure(data_specificity)
dev.off()

#correlation plot
png("disease_simulation_hist.png",height=600,width=800)
hist(permutation_simulation)
dev.off()

#permutation plot
png("permutation_plot.png",height=600,width=800)
permutation_plot(data_specificity, permutation_simulation)
dev.off()

















