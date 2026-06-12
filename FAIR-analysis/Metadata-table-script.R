library("tidyverse")
# ============================================= # 
# GSE65099
# ============================================= # 

glimpse(GSE65099_meta) # To see case or control

GSE65099_sample_meta <- data.frame(
  sample_id = c("GSM1587300", "GSM1587301",
                "GSM1587302", "GSM1587303",
                "GSM1587304", "GSM1587305",
                "GSM1587306", "GSM1587307",
                "GSM1587308", "GSM1587309", 
                "GSM1587310", "GSM1587311", 
                "GSM1587312", "GSM1587313",
                "GSM1587314", "GSM1587316",
                "GSM1587317", "GSM1587318",
                "GSM1587319"),
  condition = c(rep(("(Repeated) Implantation Failure"), 10),
                rep("RPL", 9)), # "Infertile" "RPL"
  BMI = c(25, 25, 24, 21, 22, 24.5, 22, 33, 31, 18,
          26, 22, 23, 25, 24, 32, 28, 26, NA),
  AGE = c(38, 40, 31, 32, 44, 37, 40, 38, 31, 35,
          38, 38, 38, 40, 38, 39, 39, 31, 38),
  SAMPLE_TIMEFRAME = c("LH+8", "LH+8", "LH+7", "LH+8", "LH+7", "LH+9", "LH+6", "LH+7", "LH+7", "LH+9", 
                       "LH+9", "LH+10", "LH+8", "LH+10", "LH+6", "LH+7", "LH+6", "LH+8", "LH+8"),
  SAMPLE_TIMEFRAME_DAYS = c(21, 21, 20, 21, 20, 22, 19, 20, 20, 22,
                            22, 23, 21, 23, 19, 20, 19, 21, 21),
  SAMPLE_TIMEFRAME_WEEKS = c(rep(NA, 19)),
  tissue_type = "Endometrium",
  tissue_ontology_id = "UBERON:0001295",
  Pregnancy_Status = "Not Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C82475") %>% 
  mutate(
    condition_ontology_id = case_when(
      condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
      condition == "Infertile" ~ "HP:0008222", # Female Infertility
      condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
      condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
      TRUE                     ~ NA_character_
      )
    )

glimpse(GSE65099_sample_meta)
# ============================================= # 
# GSE113790
# ============================================= # 
glimpse(GSE113790_meta) # To see case or control


GSE113790_sample_meta <- data.frame(
  sample_id = c("GSM3119483", "GSM3119484", 
                "GSM3119485", "GSM3119486",
                "GSM3119487", "GSM3119488"),
  condition = c("RPL", "RPL", 
                "RPL", "Control",
                "Control", "Control"),
  BMI = c(rep(NA, 6)),
  AGE = c(25, 23,
          24, 24,
          24, 25),
  SAMPLE_TIMEFRAME = c(rep("7.083 ± 1.567 gestational weeks", 6)),
  SAMPLE_TIMEFRAME_DAYS = c(rep(NA, 6)),
  SAMPLE_TIMEFRAME_WEEKS = c(rep(7, 6)),
  tissue_type = "Decidua",
  tissue_ontology_id = "UBERON:0002450",
  Pregnancy_Status = "Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C124295") %>% 
  mutate(
    condition_ontology_id = case_when(
      condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
      condition == "Infertile" ~ "HP:0008222", # Female Infertility
      condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
      condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
      TRUE                     ~ NA_character_
    )
  )

# ============================================= # 
# GSE161969
# ============================================= # 
glimpse(GSE161969_meta) # To see case or control

GSE178535_sample_meta <- data.frame(
  sample_id = c("GSM4928875", "GSM4928876",
                "GSM4928877", "GSM4928878",
                "GSM4928879", "GSM4928880",
                "GSM4928881"),
  condition = c("Control", "Control",
                "Control", 
                "RPL", "RPL",
                "RPL","RPL"),
  BMI = c(rep(NA, 7)),
  AGE = c(rep(NA, 7)),
  SAMPLE_TIMEFRAME = c(rep("28-82 gestational days", 7)),
  SAMPLE_TIMEFRAME_DAYS = c(rep("28-82 gestational days", 7)),
  SAMPLE_TIMEFRAME_WEEKS = c(rep("4-11 gestational weeks", 7)),
  tissue_type = "Decidua",
  tissue_ontology_id = "UBERON:0002450",
  Pregnancy_Status = "Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C124295") %>% 
  mutate(
    condition_ontology_id = case_when(
    condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
    condition == "Infertile" ~ "HP:0008222", # Female Infertility
    condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
    condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
    TRUE                     ~ NA_character_   
  )
)

glimpse(GSE178535_sample_meta)
# ============================================= # 
# GSE178535
# ============================================= # 
glimpse(GSE178535_meta) # To see case or control
print(GSE178535_meta$characteristics_ch1.1)
GSE178535_sample_meta <- data.frame(
  sample_id = c("GSM5393469",
                "GSM5393470",
                "GSM5393471", 
                "GSM5393472",
                "GSM5393473", 
                "GSM5393474"),
  condition = c(rep("Control", 3),
                rep("RPL", 3)),
  BMI = c(rep(NA, 6)),
  AGE = c(30,
          31,
          34,
          33,
          33,
          33),
  SAMPLE_TIMEFRAME = c("sample ga/(d): 62",
                       "sample ga/(d): 57",
                       "sample ga/(d): 56",
                       "sample ga/(d): 55",
                       "sample ga/(d): 61",
                       "sample ga/(d): 54"),
  SAMPLE_TIMEFRAME_DAYS = c(62,
                            57,
                            56,
                            55,
                            61,
                            54),
  SAMPLE_TIMEFRAME_WEEKS = c(8.85,
                             8.14,
                             8,
                             7.86,
                             8.71,
                             7.71),
  tissue_type = "Decidua",
  tissue_ontology_id = "UBERON:0002450",
  Pregnancy_Status = "Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C124295") %>% 
  mutate(
    condition_ontology_id = case_when(
      condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
      condition == "Infertile" ~ "HP:0008222", # Female Infertility
      condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
      condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
      TRUE                     ~ NA_character_   
    )
  )

# ============================================= # 
# GSE183555
# ============================================= # 
glimpse(GSE183555_meta) # To see case or control
GSE183555_sample_meta <- data.frame(
  sample_id = c("GSM5591871", "GSM5591872",
                "GSM5591873", "GSM5591874",
                "GSM5591875", "GSM5591876",
                "GSM5591877", "GSM5591878",
                "GSM5591879", "GSM5591880"),
  condition = c("Control", "Control", 
                "Control", "Control",
                "Control",
                "RPL", "RPL", 
                "RPL", "RPL", 
                "RPL"),
  BMI = c(rep(NA, 10)),
  AGE = c(rep(NA, 10)),
  SAMPLE_TIMEFRAME = c(rep("window of implantation", 10)),
  SAMPLE_TIMEFRAME_DAYS = c(rep(14, 10)),
  SAMPLE_TIMEFRAME_WEEKS = c(rep(NA, 10)),
  tissue_type = "Endometrial Gland",
  tissue_ontology_id = "UBERON:0002451",
  Pregnancy_Status = "Not Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C82475") %>% 
  mutate(
    condition_ontology_id = case_when(
      condition == "RPL"       ~ "MESH:D000033", # Recurrent Pregnancy Loss
      condition == "Infertile" ~ "MESH:D007247", # Female Infertility
      condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy
      TRUE                     ~ NA_character_
    )
  )
glimpse(GSE183555_sample_meta)
# ============================================= # 
# GSE22490
# ============================================= # 
glimpse(GSE22490_meta_data)
print(GSE22490_meta_data$title)
GSE22490_sample_meta <- data.frame(
  sample_id = c("GSM558664", "GSM558665",
                "GSM558666", "GSM558667", 
                "GSM558668", "GSM558669",
                "GSM558670", "GSM558671",
                "GSM558672", "GSM558673"),
  condition = c(rep("Control", 6), 
                rep("RPL", 4)),
  BMI = c(rep(NA, 10)),
  AGE = c(rep(NA, 10)),
  SAMPLE_TIMEFRAME = c("8 weeks", "8 weeks", 
                       "8 weeks", "4 weeks",
                       "11 weeks", "13 weeks",
                       "14.6 weeks", "6.6 weeks",
                       "6.6 weeks", "8 weeks"),
  SAMPLE_TIMEFRAME_DAYS = c(rep(NA, 10)),
  SAMPLE_TIMEFRAME_WEEKS = c(8, 8,
                             8, 4,
                             11, 13,
                             14.6, 6.6,
                             6.6, 8),
  tissue_type = "Placenta",
  tissue_ontology_id = "UBERON:0001987",
  Pregnancy_Status = "Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C124295") %>% 
  mutate(
    condition_ontology_id = case_when(
      condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
      condition == "Infertile" ~ "HP:0008222", # Female Infertility
      condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
      condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
      TRUE                     ~ NA_character_   
    )
  )
glimpse(GSE22490_sample_meta)
print(GSE22490_meta_data$geo_accession)
# ============================================= # 
# GSE26787
# ============================================= # 
glimpse(GSE26787_meta_data)

GSE26787_sample_meta <- data.frame(
  sample_id = c("GSM659103","GSM659104", "GSM659105",
                "GSM659106", "GSM659107", "GSM659108",
                "GSM659109", "GSM659110", "GSM659111",
                "GSM659112", "GSM659113", "GSM659114",
                "GSM659115", "GSM659116", "GSM659117"),
  condition = c(rep("Control", 5),
                rep("(Repeated) Implantation Failure", 5),
                rep("RPL", 5)),
  BMI = c(rep(NA, 15)),
  AGE = c(rep(NA, 15)),
  SAMPLE_TIMEFRAME = c(rep("7-9 days after ovulation (middle luteal phase)", 15)),
  SAMPLE_TIMEFRAME_DAYS = c(rep(22, 15)),
  SAMPLE_TIMEFRAME_WEEKS = c(rep(NA, 15)),
  tissue_type = "Endometrial Gland",
  tissue_ontology_id = "UBERON:0002451",
  Pregnancy_Status = "Not Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C82475") %>% 
  mutate(
    condition_ontology_id = case_when(
      condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
      condition == "Infertile" ~ "HP:0008222", # Female Infertility
      condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
      condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
      TRUE                     ~ NA_character_   
    )
  )
glimpse(GSE26787_sample_meta)
# ============================================= # 
# GSE165004
# ============================================= #
glimpse(GSE165004_meta_data)

# print(GSE165004_meta_data$geo_accession)
GSE165004_rownames = c("GSM5024320", "GSM5024321", "GSM5024322",
 "GSM5024323", "GSM5024324", "GSM5024325",
 "GSM5024326", "GSM5024327", "GSM5024328",
 "GSM5024329", "GSM5024330", "GSM5024331",
 "GSM5024332", "GSM5024333", "GSM5024334",
 "GSM5024335", "GSM5024336", "GSM5024337",
 "GSM5024338", "GSM5024339", "GSM5024340",
 "GSM5024341", "GSM5024342", "GSM5024343",
 "GSM5024344", "GSM5024345", "GSM5024346",
 "GSM5024347", "GSM5024348", "GSM5024349",
 "GSM5024350", "GSM5024351", "GSM5024352",
 "GSM5024353", "GSM5024354", "GSM5024355",
 "GSM5024356", "GSM5024357", "GSM5024358",
 "GSM5024359", "GSM5024360", "GSM5024361",
 "GSM5024362", "GSM5024363", "GSM5024364",
 "GSM5024365", "GSM5024366", "GSM5024367",
 "GSM5024368", "GSM5024369", "GSM5024370",
 "GSM5024371", "GSM5024372", "GSM5024373",
 "GSM5024374", "GSM5024375", "GSM5024376",
 "GSM5024377", "GSM5024378", "GSM5024379",
 "GSM5024380", "GSM5024381", "GSM5024382",
 "GSM5024383", "GSM5024384", "GSM5024385",
 "GSM5024386", "GSM5024387", "GSM5024388",
 "GSM5024389", "GSM5024390", "GSM5024391")
all(GSE165004_rownames == GSE165004_meta_data$geo_accession)
rm(GSE165004_rownames)

GSE165004_sample_meta <- data.frame(
  sample_id = c("GSM5024320", "GSM5024321", "GSM5024322",
                "GSM5024323", "GSM5024324", "GSM5024325",
                "GSM5024326", "GSM5024327", "GSM5024328",
                "GSM5024329", "GSM5024330", "GSM5024331",
                "GSM5024332", "GSM5024333", "GSM5024334",
                "GSM5024335", "GSM5024336", "GSM5024337",
                "GSM5024338", "GSM5024339", "GSM5024340",
                "GSM5024341", "GSM5024342", "GSM5024343",
                "GSM5024344", "GSM5024345", "GSM5024346",
                "GSM5024347", "GSM5024348", "GSM5024349",
                "GSM5024350", "GSM5024351", "GSM5024352",
                "GSM5024353", "GSM5024354", "GSM5024355",
                "GSM5024356", "GSM5024357", "GSM5024358",
                "GSM5024359", "GSM5024360", "GSM5024361",
                "GSM5024362", "GSM5024363", "GSM5024364",
                "GSM5024365", "GSM5024366", "GSM5024367",
                "GSM5024368", "GSM5024369", "GSM5024370",
                "GSM5024371", "GSM5024372", "GSM5024373",
                "GSM5024374", "GSM5024375", "GSM5024376",
                "GSM5024377", "GSM5024378", "GSM5024379",
                "GSM5024380", "GSM5024381", "GSM5024382",
                "GSM5024383", "GSM5024384", "GSM5024385",
                "GSM5024386", "GSM5024387", "GSM5024388",
                "GSM5024389", "GSM5024390", "GSM5024391"),
  condition = c(rep("Control", 24),
                rep("RPL", 24),
                rep("(Repeated) Implantation Failure", 24)),
  BMI = c(rep(NA, 72)),
  AGE = c(rep(NA, 72)),
  SAMPLE_TIMEFRAME = c(rep("Mid-luteal phase", 72)),
  SAMPLE_TIMEFRAME_DAYS = c(rep(22, 72)),
  SAMPLE_TIMEFRAME_WEEKS = c(rep(NA, 72)),
  tissue_type = "Endometrial Gland",
  tissue_ontology_id = "UBERON:0002451",
  Pregnancy_Status = "Not Pregnant",
  Pregnancy_Status_ontology_id = "NCIT:C82475") %>% 
  mutate(
    condition_ontology_id = case_when(
    condition == "RPL"       ~ "HP:0200067", # Recurrent Pregnancy Loss
    condition == "Infertile" ~ "HP:0008222", # Female Infertility
    condition == "Control"   ~ "EFO:0001461",  # Experimental Control / Healthy,
    condition == "(Repeated) Implantation Failure" ~ "HP:0033712",
    TRUE                     ~ NA_character_)
    )

rownames(GSE165004_sample_meta)
# ============================================= # 
# Metadata_table

# write.table(NAME, file = " _SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)

#
Folder_of_your_choice <- rstudioapi::selectDirectory()
setwd(Folder_of_your_choice)
list.files()
# GSE65099_sample_meta
exists("GSE65099_sample_meta")
write.table(GSE65099_sample_meta, file = "GSE65099_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=TRUE)


# GSE113790_sample_meta
exists("GSE113790_sample_meta")
write.table(GSE113790_sample_meta, file = "GSE113790_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)


# GSE178535_sample_meta
exists("GSE178535_sample_meta")
write.table(GSE178535_sample_meta, file = "GSE178535_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)


# GSE183555_sample_meta
exists("GSE183555_sample_meta")
write.table(GSE183555_sample_meta, file = "GSE183555_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)


# GSE22490_sample_meta
exists("GSE22490_sample_meta")
write.table(GSE22490_sample_meta, file = "GSE22490_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)

# GSE26787_sample_meta
exists("GSE26787_sample_meta")
write.table(GSE26787_sample_meta, file = "GSE26787_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)


# GSE165004_sample_meta
exists("GSE165004_sample_meta")
write.table(GSE165004_sample_meta, file = "GSE165004_SAMPLE_Metadata_table.csv", sep = ";", dec = ".", row.names=FALSE)
