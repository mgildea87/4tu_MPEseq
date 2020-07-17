
library(minpack.lm)
library(nlstools)
library(nlme)

# Load in data
setwd('//cbsumbgfs1.biohpc.cornell.edu/storage/MBG-LAB-Pleiss/mag456/4-thio U labeling/Time course 5.2019/data/Pipeline_V9/')
SI_table = read.delim('Combined_SI.txt', header = T)
row.names(SI_table) = SI_table[,1]
SI_table = SI_table[,2:31]

Total_counts_table = read.delim('Combined_total_counts.txt', header = T)
row.names(Total_counts_table) = Total_counts_table[,1]
Total_counts_table = Total_counts_table[,2:31]

Total_unspliced_table = read.delim('Combined_concordant_unspliced_counts.txt', header = T)
row.names(Total_unspliced_table) = Total_unspliced_table[,1]
Total_unspliced_table = Total_unspliced_table[,2:31]


Total_pre1st_step_table = read.delim('pre1st_step_adjusted.txt', header = T)
row.names(Total_pre1st_step_table) = Total_pre1st_step_table[,1]
Total_pre1st_step_table = Total_pre1st_step_table[,2:31]

Total_branched_table = read.delim('branched_adjusted.txt', header = T)
row.names(Total_branched_table) = Total_branched_table[,1]
Total_branched_table = Total_branched_table[,2:31]

Total_spliced_table = read.delim('Combined_spliced_counts.txt', header = T)
row.names(Total_spliced_table) = Total_spliced_table[,1]
Total_spliced_table = Total_spliced_table[,2:31]

#normalization factor
spike_in_counts = colSums(Total_counts_table[307:311,])
spike_in_counts = spike_in_counts/colSums(Total_counts_table)
norm_fac  = c(spike_in_counts[1]/spike_in_counts[1:10], spike_in_counts[11]/spike_in_counts[11:20],spike_in_counts[21]/spike_in_counts[21:30])
time_points = c(0,92,124,160,210,304,392,543,725,990,0,90,120,150,209,301,390,544,720,984,0,91,126,159,217,298,392,544,728,993)


spike_1 = as.numeric(Total_counts_table[307,])
spike_1 = spike_1/colSums(Total_counts_table)
spike_1  = c(spike_1[1]/spike_1[1:10], spike_1[11]/spike_1[11:20], spike_1[21]/spike_1[21:30])


spike_2 = as.numeric(Total_counts_table[308,])
spike_2 = spike_2/colSums(Total_counts_table)
spike_2  = c(spike_2[1]/spike_2[1:10], spike_2[11]/spike_2[11:20], spike_2[21]/spike_2[21:30])


spike_3 = as.numeric(Total_counts_table[309,])
spike_3 = spike_3/colSums(Total_counts_table)
spike_3  = c(spike_3[1]/spike_3[1:10], spike_3[11]/spike_3[11:20], spike_3[21]/spike_3[21:30])


spike_4 = as.numeric(Total_counts_table[310,])
spike_4 = spike_4/colSums(Total_counts_table)
spike_4  = c(spike_4[1]/spike_4[1:10], spike_4[11]/spike_4[11:20], spike_4[21]/spike_4[21:30])


spike_5 = as.numeric(Total_counts_table[311,])
spike_5 = spike_5/colSums(Total_counts_table)
spike_5  = c(spike_5[1]/spike_5[1:10], spike_5[11]/spike_5[11:20], spike_5[21]/spike_5[21:30])

d = data.frame(time = time_points, spike_1 = spike_1, spike_2 = spike_2, spike_3 = spike_3, spike_4 = spike_4, spike_5 = spike_5)

#remove outliers

d_wooutliers = d
for (i in 1:10){
  out = boxplot.stats(as.numeric(c(d[i,2:6], d[10+i,2:6], d[2+i,2:6])))$out
  if (length(out)>0){
    for (i in 1:length(out)){
      d_wooutliers[which(d_wooutliers == out[i], arr.ind = T)] = NA
    }
  }  
}
#I'm considering these ones below outliers as well.
d_wooutliers[7,3] = NA
d_wooutliers[8,3] = NA
d_wooutliers[9,6] = NA
d_wooutliers[29,3] = NA
d_wooutliers[29,4] = NA
d_wooutliers[29,5] = NA

#compute normalization factor

norm_fac = vector()
for (i in 1:10){
  norm_fac = c(norm_fac, exp(mean(log(na.omit(as.numeric(c(d_wooutliers[i,2:6], d_wooutliers[10+i,2:6], d_wooutliers[20+i,2:6])))))))
}  
norm_fac = c(rep(norm_fac,3))


#normalize counts and establish introns to remove from further analysis


Total_unspliced_table_norm = mapply('*', Total_unspliced_table, norm_fac)
Total_unspliced_back = apply(Total_unspliced_table_norm, MARGIN = 1, function(x) mean(c(x[1],x[11],x[21])))
Total_unspliced_table_norm = cbind(Total_unspliced_table_norm[,1:10]-Total_unspliced_back,Total_unspliced_table_norm[,11:20]-Total_unspliced_back,Total_unspliced_table_norm[,21:30]-Total_unspliced_back)
Total_unspliced_table_norm[which(Total_unspliced_table_norm < 0)] = 0

Total_pre1st_step_table_norm = mapply('*', Total_pre1st_step_table, norm_fac)
Total_pre1st_back = apply(Total_pre1st_step_table_norm, MARGIN = 1, function(x) mean(c(x[1],x[11],x[21])))
Total_pre1st_step_table_norm= cbind(Total_pre1st_step_table_norm[,1:10]-Total_pre1st_back,Total_pre1st_step_table_norm[,11:20]-Total_pre1st_back,Total_pre1st_step_table_norm[,21:30]-Total_pre1st_back)
Total_pre1st_step_table_norm[which(Total_pre1st_step_table_norm < 0)] = 0

Total_branched_table_norm = mapply('*', Total_branched_table, norm_fac)
Total_branched_back = apply(Total_branched_table_norm, MARGIN = 1, function(x) mean(c(x[1],x[11],x[21])))
Total_branched_table_norm= cbind(Total_branched_table_norm[,1:10]-Total_branched_back,Total_branched_table_norm[,11:20]-Total_branched_back,Total_branched_table_norm[,21:30]-Total_branched_back)
Total_branched_table_norm[which(Total_branched_table_norm < 0)] = 0

Total_spliced_table_norm = mapply('*', Total_spliced_table, norm_fac)
Total_spliced_back = apply(Total_spliced_table_norm, MARGIN = 1, function(x) mean(c(x[1],x[11],x[21])))
Total_spliced_table_norm= cbind(Total_spliced_table_norm[,1:10]-Total_spliced_back,Total_spliced_table_norm[,11:20]-Total_spliced_back,Total_spliced_table_norm[,21:30]-Total_spliced_back)
Total_spliced_table_norm[which(Total_spliced_table_norm < 0)] = 0

#these are mostly annotated predicted introns with no evidence of splicing in our data. 
erroneous_genes = c("SNR17A;1","SNR17B;1",'SPBP8B7.06;1','SPBC800.04c;1','SPBC18E5.06;1','SPAC1805.13;1','SPAC15E1.03;1',
                    'YPR202W;1',
                    'YPL283C;1',
                    'YOR318C;1',
                    'YNL339C;1',
                    'YML133C;1',
                    'YLR464W;1',
                    'YLR202C;1',
                    'YLL067C;1',
                    'YLL066C;1',
                    'YJR112W-A;1',
                    'YJR079W;1',
                    'YHR218W;1',
                    'YGR296W;1',
                    'YDR424C;1',
                    'YEL076C-A;1',
                    'YDR535C;1',
                    'YBL111C;1',
                    'YJL225C;1',
                    'YIL177C;1',
                    'YHL050C;1',
                    'YBR219C;1',
                    'YLR054C;1'
)


#t_off calculation
#0 and 90 second time points are removed from her on and Several outlier samples are removed from here on.  

t = c(time_points[3:6], time_points[9:10], time_points[13:15], time_points[17:20], time_points[23:28], time_points[30])
x = vector()

for (i in 1:nrow(Total_pre1st_step_table)){
  gene = i
  if(row.names(Total_pre1st_step_table)[i] %in% erroneous_genes){
    next
  }
  y = c(Total_pre1st_step_table_norm[gene,3:6], Total_pre1st_step_table_norm[gene, 9:10], Total_pre1st_step_table_norm[gene,13:15], Total_pre1st_step_table_norm[gene, 17:20],Total_pre1st_step_table_norm[gene, 23:28], Total_pre1st_step_table_norm[gene, 30])
  m = NA
  tr <- tryCatch(
    m<-nlsLM(y ~ u/s * (1-(exp(-s*(t-x)))), start = list(u = 500, s = .0005, x = 90), lower = c(0,0,0), weights = 1/t^2), 
    error = function(e) e) 
  if(inherits(tr, "error")){
    x = c(x, NA)
  }
  else if(length(which(y == 0)) > 4){
    x = c(x, NA)
  }
  else{
    x = c(x, coef(m)[3])
  }
}
t_off = median(na.omit(x))



#combined rates models

setwd('//cbsumbgfs1.biohpc.cornell.edu/storage/MBG-LAB-Pleiss/mag456/4-thio U labeling/Time course 5.2019/data/Pipeline_V9/Updated modeling method/')
colnames = c("synthesis_coupled_model",
             "synthesis_coupled_model_90%CI_lower", 
             "synthesis_coupled_model_90%CI_upper", 
             "step1_rate_coupled_model",
             "step1_rate_coupled_model_90%CI_lower",
             "step1_rate_coupled_model_90%CI_upper",
             "step2_rate_coupled_model",
             "step2_rate_coupled_model_90%CI_lower",
             "step2_rate_coupled_model_90%CI_upper")
Rates_mat = matrix(nrow = nrow(Total_unspliced_table), ncol = length(colnames))
colnames(Rates_mat) <- colnames
row.names(Rates_mat) <- row.names(Total_unspliced_table)

t = c(time_points[3:6], time_points[9:10], time_points[13:15], time_points[17:20], time_points[23:28], time_points[30])

for (i in 1:nrow(Total_unspliced_table)){
  gene = i
  if(row.names(Rates_mat)[i] %in% erroneous_genes){
    Rates_mat[i,1:9] = NA
    next
  }
  y = c(Total_pre1st_step_table_norm[gene,3:6], Total_pre1st_step_table_norm[gene, 9:10], Total_pre1st_step_table_norm[gene,13:15], Total_pre1st_step_table_norm[gene, 17:20],Total_pre1st_step_table_norm[gene, 23:28], Total_pre1st_step_table_norm[gene, 30])
  x = c(Total_branched_table_norm[gene,3:6], Total_branched_table_norm[gene, 9:10], Total_branched_table_norm[gene,13:15], Total_branched_table_norm[gene, 17:20],Total_branched_table_norm[gene, 23:28], Total_branched_table_norm[gene, 30])
  d = NA
  p = NA
  d = data.frame(val = c(y,x), time = rep(t,2), ispre = c(rep(1,length(y)), rep(0,length(y))),isbra = c(rep(0,length(x)), rep(1,length(x))))
  tr <- tryCatch(
    p <- nlsLM(val ~ isbra*((u/(u2*(u2-u1))) * (u2*(1-exp(-u1*(time-t_off))) - u1*(1-exp(-u2*(t-t_off))))) + 
    ispre*(u/u1 * (1-(exp(-u1*(time-t_off))))), start = list(u = 1000, u1 = 0.05, u2 = .005), lower = c(0,0,0), 
    weights = 1/time^2, data = d), 
    error = function(e) e)
  if(inherits(tr, "error")){
    Rates_mat[i,1:9] = NA
  }
  else if(length(which(y == 0)) > 3){
    Rates_mat[i,1:9] = NA
  }
  else if(length(which(x == 0)) > 3){
    Rates_mat[i,1:9] = NA
  }
  else{
    Rates_mat[i,1] = coef(p)[1]
    Rates_mat[i,2] = confint2(p, level = 0.9)[1,1]
    Rates_mat[i,3] = confint2(p, level = 0.9)[1,2]
    Rates_mat[i,4] = log(2)/coef(p)[2]
    Rates_mat[i,5] = log(2)/confint2(p, level = 0.9)[2,1]
    Rates_mat[i,6] = log(2)/confint2(p, level = 0.9)[2,2]
    Rates_mat[i,7] = log(2)/coef(p)[3]
    Rates_mat[i,8] = log(2)/confint2(p, level = 0.9)[3,1]
    Rates_mat[i,9] = log(2)/confint2(p, level = 0.9)[3,2]
  }
}
write.csv(Rates_mat, "Combined_rates_coupled_model_WT.csv")
