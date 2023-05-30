library(here)

source(here("base","src.R"))

#read host data
meta <- read_csv(here("data/raw_data/tank_metadata.csv"))
counts <- read_csv(here("data/raw_data/tank_counts.csv"))


#read ysi data
ysi_location <- list.files(path=here("data/raw_data/ysi_data"))
ysi <- readr::read_csv(here("data","raw_data","ysi_data",ysi_location), id = "filename", locale = locale(encoding = "latin1"))

#write ysi as singular csv
#write.csv(ysi, file = here("data/raw_data/ysi.csv"))

#declare object for experiment start date
start_date <- mdy("06/27/2022")

#fix up names for ysi data
ysi %<>% rename('temp_f' = '°F-21A104910', "do_per" = "DO %-21A106679", "do_conc" = "DO mg/L-21A106679",
                "ph" = "pH-21B101988", "ph_mv" = "pH mV-21B101988", "chlor_rfu" = "Chlorophyll RFU-21B101736",
                "chlor_conc" = "Chlorophyll ug/L-21B101736")
ysi %<>% mutate(bucket = as.factor(as.numeric(gsub("Tank ", "", Site))))
ysi %<>% mutate(temp_c = (temp_f-32)*(5/9))
ysi$Date <- mdy(ysi$Date)


#reclass variables in meta and counts data
meta$bucket <- as.factor(meta$bucket)
counts$bucket <- as.factor(counts$bucket)


#add metadata to count data
data <- left_join(counts, meta, by = "bucket")
data$community <- as.factor(data$community)
data$sample_date <- mdy(data$sample_date)




#zoop abundance data
abundance <- data %>% mutate(total = rowSums(across(juvenile_daphnia_uninf:cerio_inf), na.rm = T), 
                             daphnia = rowSums(across(juvenile_daphnia_uninf:adult_daphnia_inf), na.rm = T),
                             juve = juvenile_daphnia_inf + juvenile_daphnia_uninf,
                             adult = adult_daphnia_inf + adult_daphnia_uninf,
                             cerio = cerio_inf + cerio_uninf,
                             daph_ratio = daphnia/total,
                             adult_ratio = adult/daphnia,
                             juve_ratio = juve/daphnia,
                             treatment = paste(community, temp, sep="_"))

summ_abund <- abundance %>% group_by(treatment, sample_date, temp, community) %>% summarize(mean_juve = mean(juve),
                                                                                            juve_var = var(juve),
                                                                                            juve_se = sqrt(var(juve)/length(juve)),
                                                                                            mean_adult = mean(adult),
                                                                                            adult_var = var(adult),
                                                                                            adult_se = sqrt(var(adult)/length(adult)),
                                                                                            mean_cerio = mean(cerio),
                                                                                            cerio_var = var(cerio),
                                                                                            cerio_se = sqrt(var(cerio)/length(cerio)),
                                                                                            mean_daphnia = mean(daphnia),
                                                                                            daphnia_var = var(daphnia),
                                                                                            daphnia_se = sqrt(var(daphnia)/length(daphnia)),
                                                                                            mean_total = mean(total),
                                                                                            total_var = var(total),
                                                                                            total_se = sqrt(var(total)/length(total)))

#zoop infection data
infections <- abundance %>% mutate(daphnia_inc = (juvenile_daphnia_inf+adult_daphnia_inf),
                                   cerio_inc = (cerio_inf),
                                   total_inc = (juvenile_daphnia_inf+adult_daphnia_inf+cerio_inf),
                                   daphnia_prev = (juvenile_daphnia_inf+adult_daphnia_inf)/daphnia,
                                   cerio_prev = (cerio_inf)/cerio,
                                   total_prev = (juvenile_daphnia_inf+adult_daphnia_inf+cerio_inf)/total)

summ_inf <- infections %>% group_by(treatment, sample_date, temp, community) %>% 
  summarize(n = n(),
            mean_inc = mean(total_inc, na.rm = T),
            inc_se = sqrt(var(total_inc, na.rm = T)/length(total_inc)),
            mean_daph_inc = mean(daphnia_inc, na.rm = T),
            daph_inc_se = sqrt(var(daphnia_inc, na.rm = T)/length(daphnia_inc)),
            mean_cer_inc = mean(cerio_inc, na.rm = T),
            cer_inc_se = sqrt(var(cerio_inc, na.rm = T)/length(cerio_inc)),
            daph_prev = mean(daphnia_prev),
            cerio_prev = mean(cerio_prev),
            total_prev = mean(total_prev))



summ_ysi <- ysi %>% group_by(bucket, Date) %>% summarize(temp_f = mean(temp_f),
                                                         temp_c = mean(temp_c),
                                                         do_per = mean(do_per),
                                                         do_conc = mean(do_conc),
                                                         ph = mean(ph),
                                                         ph_mv = mean(ph_mv),
                                                         chlor_rfu = mean(chlor_rfu),
                                                         chlor_conc = mean(chlor_conc))
summ_ysi %<>% left_join(.,meta, by = "bucket")



saveRDS(abundance, file = here("data/processed_data","abundance.rds"))
saveRDS(infections, file = here("data/processed_data","infections.rds"))
saveRDS(summ_abund, file = here("data/processed_data","summ_abund.rds"))
saveRDS(summ_inf, file = here("data/processed_data","summ_inf.rds"))
saveRDS(summ_ysi, file = here("data/processed_data","summ_ysi.rds"))



