
setwd('C:/Users/mike/Desktop/PWM_cere_splicesites')
library(ggplot2)
library(ggseqlogo)
intron_seq = read.delim('cere_intron_sequences.txt', header = F, stringsAsFactors = F)
BP_seq = read.delim('cere_BP_sequences.txt', header = F, stringsAsFactors = F)
RPGs = read.delim('RPGs.txt', header = T, stringsAsFactors = F)
Non_RPGs = read.delim('Non_RPGs.txt', header = T, stringsAsFactors = F)
intron_seq = intron_seq[order(intron_seq$V1),]
BP_seq = BP_seq[order(BP_seq$V1),]

#compute background base frequency
A = 0
Te = 0
C = 0
G = 0
l = 0
for (i in 1:nrow(intron_seq)){
  seq = intron_seq$V2[i]
  A = A+lengths(regmatches(seq, gregexpr("A", seq)))
  Te = Te+lengths(regmatches(seq, gregexpr("T", seq)))
  C = C+lengths(regmatches(seq, gregexpr("C", seq)))
  G = G+lengths(regmatches(seq, gregexpr("G", seq)))
  l = l+nchar(seq)
}
A = A/l
Te = Te/l
C = C/l
G = G/l

#five_PWM 5'SS +/- 3 bases
five_count_mat = matrix(nrow = 4, ncol = 12)
five_count_mat[] <- 0
row.names(five_count_mat) = c("A", "C", "U", "G")
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],1,12)
  for (j in 1:12){
    base = substr(seq,j,j)
    if (base == "A"){
      five_count_mat[1,j] = five_count_mat[1,j]+1
    }
    else if (base == "C"){
      five_count_mat[2,j] = five_count_mat[2,j]+1
    }
    else if (base == "T"){
      five_count_mat[3,j] = five_count_mat[3,j]+1
    }
    else if (base == "G"){
      five_count_mat[4,j] = five_count_mat[4,j]+1
    }
  }
}
five_PPM = five_count_mat[]/291 
five_PWM = five_PPM
five_PWM[1,] = log2(five_PWM[1,]/A)
five_PWM[2,] = log2(five_PWM[2,]/C)
five_PWM[3,] = log2(five_PWM[3,]/Te)
five_PWM[4,] = log2(five_PWM[4,]/G)


#3'SS 14bp upstream 3 bases downstream
three_count_mat = matrix(nrow = 4, ncol = 20)
three_count_mat[] <- 0
row.names(three_count_mat) = c("A", "C", "U", "G")
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],nchar(intron_seq$V2[i])-19,nchar(intron_seq$V2[i]))
  for (j in 1:20){
    base = substr(seq,j,j)
    if (base == "A"){
      three_count_mat[1,j] = three_count_mat[1,j]+1
    }
    else if (base == "C"){
      three_count_mat[2,j] = three_count_mat[2,j]+1
    }
    else if (base == "T"){
      three_count_mat[3,j] = three_count_mat[3,j]+1
    }
    else if (base == "G"){
      three_count_mat[4,j] = three_count_mat[4,j]+1
    }
  }
}
three_PPM = three_count_mat[]/291 
three_PWM = three_PPM
three_PWM[1,] = log2(three_PWM[1,]/A)
three_PWM[2,] = log2(three_PWM[2,]/C)
three_PWM[3,] = log2(three_PWM[3,]/Te)
three_PWM[4,] = log2(three_PWM[4,]/G)

#BP'SS +/-3 bp upstream/downstream
BP_count_mat = matrix(nrow = 4, ncol = 13)
BP_count_mat[] <- 0
row.names(BP_count_mat) = c("A", "C", "U", "G")
for (i in 1:nrow(BP_seq)){
  seq = BP_seq$V2[i]
  for (j in 1:13){
    base = substr(seq,j,j)
    if (base == "A"){
      BP_count_mat[1,j] = BP_count_mat[1,j]+1
    }
    else if (base == "C"){
      BP_count_mat[2,j] = BP_count_mat[2,j]+1
    }
    else if (base == "T"){
      BP_count_mat[3,j] = BP_count_mat[3,j]+1
    }
    else if (base == "G"){
      BP_count_mat[4,j] = BP_count_mat[4,j]+1
    }
  }
}
BP_PPM = BP_count_mat[]/291 
BP_PWM = BP_PPM
BP_PWM[1,] = log2(BP_PWM[1,]/A)
BP_PWM[2,] = log2(BP_PWM[2,]/C)
BP_PWM[3,] = log2(BP_PWM[3,]/Te)
BP_PWM[4,] = log2(BP_PWM[4,]/G)



#Splice site score calculation

scores = data.frame(intron = intron_seq$V1)

#5'SS +/- 3 bases
five = vector()
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],1,12)
  score = 0
  for (j in 1:nchar(seq)){
    base = substr(seq,j,j)
    if (base == "A"){
      score = score + five_PWM[1,j]
    }
    else if (base == "C"){
      score = score + five_PWM[2,j]
    }
    else if (base == "T"){
      score = score + five_PWM[3,j]
    }
    else if (base == "G"){
      score = score + five_PWM[4,j]
    }
  }
  five = c(five, score)
}

#3'SS 9 bases upstream of YAG 3 bases downstream (3'SS region score)
three_9up = vector()
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],nchar(intron_seq$V2[i])-14,nchar(intron_seq$V2[i]))
  score = 0
  for (j in 1:nchar(seq)){
    base = substr(seq,j,j)
    if (base == "A"){
      score = score + three_PWM[1,j+5]
    }
    else if (base == "C"){
      score = score + three_PWM[2,j+5]
    }
    else if (base == "T"){
      score = score + three_PWM[3,j+5]
    }
    else if (base == "G"){
      score = score + three_PWM[4,j+5]
    }
  }
  three_9up = c(three_9up, score)
}

#3'SS 2 bases upstream of YAG 3 bases downstream (3'SS score)
three_2up = vector()
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],nchar(intron_seq$V2[i])-7,nchar(intron_seq$V2[i]))
  score = 0
  for (j in 1:nchar(seq)){
    base = substr(seq,j,j)
    if (base == "A"){
      score = score + three_PWM[1,j+12]
    }
    else if (base == "C"){
      score = score + three_PWM[2,j+12]
    }
    else if (base == "T"){
      score = score + three_PWM[3,j+12]
    }
    else if (base == "G"){
      score = score + three_PWM[4,j+12]
    }
  }
  three_2up = c(three_2up, score)
}

#3'SS 8 bases upstream to 4 bases upstream (polyU tract)
U_tract_small = vector()
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],nchar(intron_seq$V2[i])-14,nchar(intron_seq$V2[i])-9)
  score = 0
  for (j in 1:nchar(seq)){
    base = substr(seq,j,j)
    if (base == "A"){
      score = score + three_PWM[1,j+5]
    }
    else if (base == "C"){
      score = score + three_PWM[2,j+5]
    }
    else if (base == "T"){
      score = score + three_PWM[3,j+5]
    }
    else if (base == "G"){
      score = score + three_PWM[4,j+5]
    }
  }
  U_tract_small = c(U_tract_small, score)
}

#Branch point +/- 3 bases
BP = vector()
for (i in 1:nrow(BP_seq)){
  seq = BP_seq$V2[i]
  score = 0
  for (j in 1:nchar(seq)){
    base = substr(seq,j,j)
    if (base == "A"){
      score = score + BP_PWM[1,j]
    }
    else if (base == "C"){
      score = score + BP_PWM[2,j]
    }
    else if (base == "T"){
      score = score + BP_PWM[3,j]
    }
    else if (base == "G"){
      score = score + BP_PWM[4,j]
    }
  }
  BP = c(BP, score)
}

scores$fiveSS = five
scores$threeSS9up = three_9up
scores$threeSS2up = three_2up
scores$BP = BP
scores$U_tract_small = U_tract_small
write.csv(scores, "scores.csv")



#RPG vs non-RPG at the 3'SS
#Non_RPGs

three_count_mat = matrix(nrow = 4, ncol = 20)
three_count_mat[] <- 0
row.names(three_count_mat) = c("A", "C", "U", "G")
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],nchar(intron_seq$V2[i])-19,nchar(intron_seq$V2[i]))
  if (intron_seq$V1[i] %in% Non_RPGs[,1]){
    for (j in 1:20){
      base = substr(seq,j,j)
      if (base == "A"){
        three_count_mat[1,j] = three_count_mat[1,j]+1
      }
      else if (base == "C"){
        three_count_mat[2,j] = three_count_mat[2,j]+1
      }
      else if (base == "T"){
        three_count_mat[3,j] = three_count_mat[3,j]+1
      }
      else if (base == "G"){
        three_count_mat[4,j] = three_count_mat[4,j]+1
      }
    }
  }
}
three_PPM_nonRPG = three_count_mat[]/203 
three_PWM_nonRPG = three_PPM_nonRPG
three_PWM_nonRPG[1,] = log2(three_PWM_nonRPG[1,]/A)
three_PWM_nonRPG[2,] = log2(three_PWM_nonRPG[2,]/C)
three_PWM_nonRPG[3,] = log2(three_PWM_nonRPG[3,]/Te)
three_PWM_nonRPG[4,] = log2(three_PWM_nonRPG[4,]/G)


#RPGs
three_count_mat = matrix(nrow = 4, ncol = 20)
three_count_mat[] <- 0
row.names(three_count_mat) = c("A", "C", "U", "G")
for (i in 1:nrow(intron_seq)){
  seq = substr(intron_seq$V2[i],nchar(intron_seq$V2[i])-19,nchar(intron_seq$V2[i]))
  if (intron_seq$V1[i] %in% RPGs[,1]){
    for (j in 1:20){
      base = substr(seq,j,j)
      if (base == "A"){
        three_count_mat[1,j] = three_count_mat[1,j]+1
      }
      else if (base == "C"){
        three_count_mat[2,j] = three_count_mat[2,j]+1
      }
      else if (base == "T"){
        three_count_mat[3,j] = three_count_mat[3,j]+1
      }
      else if (base == "G"){
        three_count_mat[4,j] = three_count_mat[4,j]+1
      }
    }
  }
}
three_PPM_RPG = three_count_mat[]/105 
three_PWM_RPG = three_PPM_RPG
three_PWM_RPG[1,] = log2(three_PWM_RPG[1,]/A)
three_PWM_RPG[2,] = log2(three_PWM_RPG[2,]/C)
three_PWM_RPG[3,] = log2(three_PWM_RPG[3,]/Te)
three_PWM_RPG[4,] = log2(three_PWM_RPG[4,]/G)

#Weblogo plots

ggseqlogo(five_PPM)
ggsave("five_ppm.pdf")
ggseqlogo(BP_PPM)
ggsave("BP_ppm.pdf")
ggseqlogo(three_PPM)
ggsave("three_ppm.pdf")
ggseqlogo(three_PPM_RPG)
ggsave("three_ppm_RPG.pdf")
ggseqlogo(three_PPM_nonRPG)
ggsave("three_ppm_nRPG.pdf")


# of differences from consensus 5SS and BP

BP_con = "TACTAAC"
five_SS_con = "GTATGT"

sub_dist = data.frame(gene = BP_seq[,1])
sub_dist$BP = apply(BP_seq, MARGIN = 1, function(z) mapply(function(x,y) sum(x!=y),strsplit(BP_con,""),strsplit(substr(z[2],4,10)
                                                                                                                ,"")))
sub_dist$BP_seq = apply(BP_seq, MARGIN = 1, function(z) substr(z[2],4,10))
sub_dist$BP_seq_extended = apply(BP_seq, MARGIN = 1, function(z) substr(z[2],2,12))

sub_dist$five_SS = apply(intron_seq, MARGIN = 1, function(z) mapply(function(x,y) sum(x!=y),strsplit(five_SS_con,""),strsplit(substr(z[2],4,9)
                                                                                                                              ,"")))
sub_dist$five_SS_seq = apply(intron_seq, MARGIN = 1, function(z) substr(z[2],4,9))
sub_dist$five_SS_seq_extended = apply(intron_seq, MARGIN = 1, function(z) substr(z[2],2,11))

write.csv(sub_dist, "SS_substitution_dist.csv")


