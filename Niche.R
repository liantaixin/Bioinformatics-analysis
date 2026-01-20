rm(list=ls())
library(MicroNiche)

library(reshape2)
#import
df_imported <- read.csv("D:/phD/16s分析模块/16s_sci_3/micronich/hem/hm/hemASV.csv")
sample_env_imported <- read.csv("D:/phD/16s分析模块/16s_sci_3/micronich/hem/hm/environment.csv")
#groups
sample_info_imported <- sample_env_imported$Sample
env_As_imported <- sample_env_imported$tAs
# Feinsinger(PS)
result_ps <- feinsingers.PS(
  df = df_imported,         
  R = 4,                    
  sampleInfo = sample_info_imported, 
  envInfo = env_As_imported, 
  q = 1.65                   
print(head(result_ps))

write.csv(result_ps, file = "Feinsinger's Proportional Similarity Index (PS).csv")

# Hurlbert(Bn)
result_hurlbert <- hurlberts.Bn(
  df = df_imported,
  R = 4,
  sampleInfo = sample_info_imported,
  envInfo = env_As_imported,
  q = 1.65
)
print(head(result_hurlbert))
write.csv(result_hurlbert, file = "Hurlbert's niche breadth index (Bn).csv")

# Levins(Bn)
result_levins_bn <- levins.Bn(
  df = df_imported,
  R = 4,
  sampleInfo = sample_info_imported,
  q = 1.65
)
print(head(result_levins_bn))
write.csv(result_levins_bn, file = "Levins' niche breadth index (Bn).csv")

# Levin Overlap index
overlap_levins <- levins.overlap(
  df = df_imported,
  q = 1.65
)
print(overlap_levins[1:6, 1:6])
write.csv(overlap_levins, file = "overlap_levins.csv")


# Proportional Overlap
overlap_prop <- proportional.overlap(
  df = df_imported,
  sampleInfo = sample_info_imported,
  envInfo = env_As_imported,
  q = 1.65
)
print(overlap_prop[1:6, 1:6])
write.csv(overlap_prop, file = "Proportional Overlap.csv")
