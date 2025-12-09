################## Funktioner over i DST #######################################

## ATC

### Glucose-lowering medication ---------------------------------------------

# DPP-4 inhibitors
dpp4i <- c(
  "A10BH", "A10BD07", "A10BD08", "A10BD10", "A10BD11",
  "A10BD13", "A10BD19", "A10BD21", "A10BD24", 
  "A10BD25", "A10BD27"
)

# SGLT-2 inhibitors
sglt2i <- c(
  "A10BK", "A10BD15", "A10BD16", "A10BD20", "A10BD23",
  "A10BD19", "A10BD21", "A10BD24", "A10BD25", "A10BD27"
)

# GLP-1 receptor agonists
glp1ra <- c("A10BJ", "A10AE56", "A10AE54")

# Metformin (expanded ranges)
metformin <- c(
  "A10BA", "A10BD02", "A10BD03", "A10BD05", "A10BD07", 
  "A10BD08", "A10BD10", "A10BD11",
  "A10BD13","A10BD14","A10BD15","A10BD16","A10BD17","A10BD18",
  "A10BD20", "A10BD22", "A10BD23", "A10BD25", "A10BD26", "A10BD27"
)

# Sulfonylureas
sulfonylureas <- c("A10BB", "A10BD01", "A10BD02", "A10BD04", "A10BD06")

# Insulin subclasses
fast_acting_insulin <- "A10AB"
intermediate_acting_insulin <- "A10AC"
premixed_insulin <- "A10AD"
long_acting_insulin <- "A10AE"

# Thiazolidinediones (expanded range A10BD03–A10BD06)
tzd <- c("A10BG", "A10BD03", "A10BD04", "A10BD05", "A10BD06")

# Acarbose
acarbose <- "A10BF"


### Lipid-lowering drugs ----------------------------------------------------

other_lipid_modifiers <- c("C10AB", "C10AC", "C10AD", "C10AX")


### Antiplatelet ------------------------------------------------------------

# Acetylsalicylic acid
asa <- c("A01AD05", "B01AC06", "N02BA01", "B01AC56")

# Other antiplatelets: all B01AC except B01AC06 & B01AC56
# Represented as exclusion:
other_antiplatelets <- c("B01AC")  # with exclusion applied during analysis


### Antihypertensive medication --------------------------------------------

acei_arb <- c("C09A", "C09B", "C09C", "C09D")
ccb <- c("C08C", "C08D")
thiazides <- "C03A"
beta_blockers <- "C07"
loop_diuretics <- "C03C"
aldosterone_antagonists <- "C03D"


### Other Drugs -------------------------------------------------------------

oral_anticoagulants <- c("B01AA", "B01AE", "B01AF")
nitrates <- "C01DA"
antiarrhythmics <- c("C01AA05", "C01BC04", "C01BD01")
weight_loss_drugs <- "A08A"
glucocorticoids <- "H02AB"
antidepressants <- "N06A"
antipsychotics <- "N05A"
anxiolytics <- "N05B"
benzodiazepines_sleep <- "N05CF"
inhalants <- "R03"
thyroid_therapy <- "H03"


# ICD

# Myocardial Infarction (MI)
mi_icd10 <- c("DI21", "DI23", "DI24")

# Non-MI Ischemic Heart Disease (Non-MI IHD)
non_mi_ihd_icd10 <- c("DI20", "DI25")

# Heart Failure (HF)
hf_icd10 <- c("DI110", "DI50", "DI130", "DI132")

# Hemorrhagic Stroke
hem_stroke_icd10 <- c("DI60", "DI61", "DI62", "DI691", "DI692")

# Ischemic Stroke
isch_stroke_icd10 <- c("DI63", "DI64", "DI693", "DI694", "DI698")

# Transient Cerebral Ischemia
tci_icd10 <- c("DG45")

# Atrial Fibrillation
af_icd10 <- c("DI48")

# Macrovascular Atherosclerotic Disease
# I70 excluding I702  → include "I70" but NOT "I702"
macro_athe_icd10 <- c("DI70", "DI71")

# Peripheral Arteriosclerotic Disease (PAD)
pad_icd10 <- c(
  "DI739", 
  "DI702",   # included because explicitly listed
  "DE105", "DE115", "DE125", "DE135", "DE145",
  "DI749"
)

# Lower Extremity Amputations (LEA)
major_lea_icd10 <- c("DZ896", "DZ896A", "DZ897")
medium_lea_icd10 <- c("DZ895", "DZ895A", "DZ895B")
minor_lea_icd10 <- c("DZ894", "DZ894A")

# Chronic Kidney Disease (CKD)
ckd_icd10 <- c("DE102", "DE112", "DE132", "DE142", "DN19")

# End-Stage Kidney Disease
eskd_icd10 <- c("DN185", "DZ992")

# Severe CKD
severe_ckd_icd10 <- c("DN184")

# Moderate CKD
moderate_ckd_icd10 <- c("DN183", "DN189")

# Kidney Transplant
kidney_tx_icd10 <- c("DZ940")

# Albuminuria  — ICD-10 has no codes (NPU only), so empty
albuminuria_icd10 <- character(0)

# Hypertensive Disease
htn_icd10 <- c("DI10", "DI11", "DI12", "DI13", "DI15")

# Retinopathy
retinopathy_icd10 <- c(
  "DH334", "DH350D", "DH352",
  "DH360B", "DH360D", "DH360E", "DH360F", "DH360G", "DH360H", "DH360I", "DH360J", "DH360K",
  "DH368D", "DH420", "DH431",
  "DE103", "DE113", "DE123", "DE133", "DE143"
)

# Neuropathy
neuropathy_icd10 <- c(
  "DE104", "DE114", "DE124", "DE134", "DE144",
  "DG590", "DG632", "DG990", "DG629"
)

# Hospitalization for Hypoglycemia
hypoglycemia_icd10 <- c(
  "DE100", "DE110", "DE120", "DE130", "DE140",
  "DE160", "DE161", "DE161B", "DE162",
  "DT383"
)

# Diabetic Ketoacidosis
dka_icd10 <- c("DE101", "DE111", "DE121", "DE131", "DE141")

### LAB

# HbA1c
hba1c_codes <- c("NPU03835", "NPU27300")

# Total cholesterol
total_cholesterol_codes <- c("NPU18412", "NPU01566")

# LDL cholesterol
ldl_codes <- c("NPU10171", "NPU01568", "DNK35308")

# HDL cholesterol
hdl_codes <- c("NPU10157", "NPU01567")

# Creatinine
creatinine_codes <- c("NPU04998", "NPU18016", "NPU17559", "NPU09101", "NPU01807")

# eGFR
egfr_codes <- c("DNK35131", "DNK35301", "DNK35302", "DNK35303")

# U-albumin/creatinine ratio
ualb_crea_codes <- c("NPU19661", "NPU28842", "NPU03918")

# ALAT
alat_codes <- c("NPU19651", "NPU53495")

# Kreatinkinase (CK)
ck_codes <- c("NPU19656", "NPU19750")


