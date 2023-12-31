load('../data/aedat.RDS')

aedat <- aedat %>% 
  mutate(AETOXGR = factor(AETOXGR, 0:3, labels=c("None","Mild","Moderate","Severe"))) %>% 
  group_by(USUBJID) %>%
  # End of study for patients without a severe event
  mutate(TTE_SEVERE = case_when(
    STUDYID=="PROTA" ~ 2,
    STUDYID=="PROTB" ~ 6),
    # Time of severe event for those that had one
    TTE_SEVERE = ifelse(AETOXGR=="Severe", TTE, TTE_SEVERE)
  ) %>% 
  ungroup()

aedat <- aedat %>% 
  group_by(USUBJID) %>% 
  arrange(USUBJID, TTE_SEVERE) %>% slice(1) %>%
  group_by(PBO) %>%
  mutate(Quartile = ifelse(PBO == "PBO", "PBO",
                           paste0("Q", ntile(CAVGSS, n = 4)))) %>% 
  ungroup()

ref_bwt <- median(aedat$BWT) # will want reference value later, so save as variable!

dat_mod <- aedat %>% mutate(
  BWT_norm = BWT - ref_bwt,
  PTTYPE = factor(PTTYPE, levels = c("PT1", "PT2", "HV"))
)
