
  ## Collect each genotype's Cap Mags; convert the list of lists into a data frame whose length is determined by the longer of the two lists
Cap_Mags = list(WT_Cap_Mags = as.numeric(c(A1_11_8_Cap_Mags, A1_12_4_Cap_Mags, A3_12_4_Cap_Mags, B1_10_30_Cap_Mags, B1_11_8_Cap_Mags, B1_12_11_Cap_Mags, B2_10_30_Cap_Mags, B2_11_8_Cap_Mags, B2_12_11_Cap_Mags, B3.2_12_4_Cap_Mags, B3_10_30_Cap_Mags, B3_12_4_Cap_Mags, B4_10_30_Cap_Mags, B4_12_11_Cap_Mags, C3_12_11_Cap_Mags)),
                Mut_Cap_Mags = as.numeric(c(A1_11_6_Cap_Mags, A3_11_6_Cap_Mags, B4_11_6_Cap_Mags, A4_11_11_Cap_Mags, B1_11_11_Cap_Mags, B3_11_11_Cap_Mags, B1_12_6_Cap_Mags, B3_12_6_Cap_Mags, A1_12_13_Cap_Mags, A3_12_13_Cap_Mags, A4_12_13_Cap_Mags, B2_12_13_Cap_Mags, B4_12_13_Cap_Mags, C1late_12_13_Cap_Mags, C2_12_13_Cap_Mags)))

Cap_Mags = data.frame(lapply(Cap_Mags, "length<-", max(lengths(Cap_Mags))))

colnames(Cap_Mags) = c("WT Cap Mags", "Mut Cap Mags")
  ## Export the table to the server under the directory where all the scripts, raw data, and raw traces are stored
#write.csv(Cap_Mags, file = "M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/Cap Mags.csv", row.names = F, col.names = T)


AITC_Mags = list(WT_AITC_Mags = as.numeric(c(A4_10_30_AITC_Mags, B4_10_30_AITC_Mags, C1_10_30_AITC_Mags, A1_11_8_AITC_Mags, B1_11_8_AITC_Mags, B3_11_8_AITC_Mags, A4_11_15_AITC_Mags, B2_11_15_AITC_Mags, B4_11_15_AITC_Mags, B2_12_4_AITC_Mags, B4_12_4_AITC_Mags, C1_12_4_AITC_Mags, C2_12_4_AITC_Mags, A3_12_11_AITC_Mags, B3_12_11_AITC_Mags, C1_12_11_AITC_Mags, C2_12_11_AITC_Mags)), Mut_AITC_Mags = as.numeric(c(A1_11_6_AITC_Mags, A3_11_6_AITC_Mags, A1_11_20_AITC_Mags, A3_11_20_AITC_Mags, B1_11_20_AITC_Mags, B3_11_20_AITC_Mags, B4_11_20_AITC_Mags, C1_11_20_AITC_Mags, A3_12_6_AITC_Mags, B2_12_6_AITC_Mags, C1_12_6_AITC_Mags, C2_12_6_AITC_Mags, C3_12_6_AITC_Mags, A1_12_13_AITC_Mags, B1_12_13_AITC_Mags, B3_12_13_AITC_Mags, C3_12_13_AITC_Mags)))

AITC_Mags = data.frame(lapply(AITC_Mags, "length<-", max(lengths(AITC_Mags))))

colnames(AITC_Mags) = c("WT AITC Mags", "Mut AITC Mags")
#write.csv(AITC_Mags, file = "M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/AITC Mags.csv", row.names = F, col.names = T)


Beta_Mags = list(WT_Beta_Mags = as.numeric(c(A4_11_15_B_Mags, B4_11_15_B_Mags, B1_12_4_B_Mags, B2_12_4_B_Mags, C2_12_4_B_Mags, A1_12_11_B_Mags, B1_12_11_B_Mags, B3_12_11_B_Mags, C1_12_11_B_Mags, C2_12_11_B_Mags, C3_12_11_B_Mags)), Mut_Beta_Mags = as.numeric(c(A3_11_11_B_Mags, B3_11_11_B_Mags, A1_11_20_B_Mags, A3_11_20_B_Mags, B2_11_20_B_Mags, B3_11_20_B_Mags, B3_12_6_B_Mags, B4_12_6_B_Mags, C2_12_6_B_Mags, C3_12_6_B_Mags, A3_12_13_B_Mags, A4_12_13_B_Mags, B2_12_13_B_Mags, B3_12_13_B_Mags, B4_12_13_B_Mags, C3_12_13_B_Mags)))

Beta_Mags = data.frame(lapply(Beta_Mags, "length<-", max(lengths(Beta_Mags))))
colnames(Beta_Mags) = c("WT Beta Mags", "Mut Beta Mags")
#write.csv(Beta_Mags, file = "M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/Beta Mags.csv", row.names = F, col.names = T)


CYM_Mags = list(WT_CYM_Mags = as.numeric(c(A4_10_30_CYM_Mags, B1_10_30_CYM_Mags, B2_10_30_CYM_Mags, B3_10_30_CYM_Mags, C1_10_30_CYM_Mags, A3_11_8_CYM_Mags, A4_11_8_CYM_Mags, B3_11_8_CYM_Mags, B4_11_8_CYM_Mags, A1_11_15_CYM_Mags, B2_11_15_CYM_Mags, B3_11_15_CYM_Mags, A4_12_4_CYM_Mags, C1_12_4_CYM_Mags)), Mut_CYM_Mags = as.numeric(c(B1_11_6_CYM_Mags, B2_11_6_CYM_Mags, B3_11_6_CYM_Mags, B1_11_11_CYM_Mags, B2_11_11_CYM_Mags, B4_11_11_CYM_Mags, C1_11_20_CYM_Mags, A1_12_6_CYM_Mags, A3_12_6_CYM_Mags, A4_12_6_CYM_Mags, B1_12_6_CYM_Mags, C1_12_6_CYM_Mags)))

CYM_Mags = data.frame(lapply(CYM_Mags, "length<-", max(lengths(CYM_Mags))))
colnames(CYM_Mags) = c("WT CYM Mags", "Mut CYM Mags")
#write.csv(CYM_Mags, file = "M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/CYM Mags.csv", row.names = F, col.names = T)


LY_Mags = list(WT_LY_Mags = as.numeric(c(A4_10_30_LY_Mags, B2_10_30_LY_Mags, C1_10_30_LY_Mags, A3_11_8_LY_Mags, A4_11_8_LY_Mags, B2_11_8_LY_Mags, B4_11_8_LY_Mags, A1_11_15_LY_Mags, A3_12_4_LY_Mags, B1_12_4_LY_Mags, A1_12_11_LY_Mags, A3_12_11_LY_Mags, B2_12_11_LY_Mags, B4_12_11_LY_Mags)), Mut_LY_Mags = as.numeric(c(A4_11_6_LY_Mags, B1_11_6_LY_Mags, A1_11_11_LY_Mags, B2_11_11_LY_Mags, B4_11_11_LY_Mags, A4_11_20_LY_Mags, B4_11_20_LY_Mags, B2_12_6_LY_Mags, C2_12_13_LY_Mags)))

LY_Mags = data.frame(lapply(LY_Mags, "length<-", max(lengths(LY_Mags))))
colnames(LY_Mags) = c("WT LY Mags", "Mut LY Mags")
write.csv(LY_Mags, file = "M:/Erik/Data/Calcium Imaging Analysis/Fura2 Images/Paper Data Copies/LY Mags.csv", row.names = F, col.names = T)
