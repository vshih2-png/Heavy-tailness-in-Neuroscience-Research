
library(readxl)
library(DescTools)
library(moments)
library(ggplot2)
library(ctrialsgov)
library(tidyverse)

## C Ratios ##

#study 1
data1 <- read_excel("/Users/vanita/Downloads/StatStudy/APEx_FinalData_20201021 (1).xlsx")

#control group
control1 <- data1$amyloid_Mean6_CB[data1$apx_group == "NoSupport"]

q_control1 <- quantile(control1, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control1 <- q_control1["95%"] / q_control1["75%"]

C95_control1

#treatment group
treat1 <- data1$amyloid_Mean6_CB[data1$apx_group == "AExSupport"]

q_treat1 <- quantile(treat1, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat1 <- q_treat1["95%"] / q_treat1["75%"]

C95_treat1

#study 2
data2 <-  read_excel("/Users/vanita/Downloads/StatStudy/Cucos et al., MYD88_mRNAlevels_Mice.xlsx")
#control group
control2 <- data2$`MYD88_T0 (baseline)`[data2$Case == "WT"]

q_control2 <- quantile(control2, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control2 <- q_control2["95%"] / q_control2["75%"]

C95_control2

#treatment group
treat2 <- data2$`MYD88_T0 (baseline)`[data2$Case == "AT_transgenic"]

q_treat2 <- quantile(treat2, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat2 <- q_treat2["95%"] / q_treat2["75%"]

C95_treat2

#study 3
data3 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/DATA _ SUBMISSION 10.27.  25 - .xlsx",
  sheet = "y MAZE "
)
#control group
control3 <- data3$'D/A + Veh'

q_control3 <- quantile(control3, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control3 <- q_control3["95%"] / q_control3["75%"]

C95_control3

#treatment group
treat3 <- data3$`D/A +INS`

q_treat3 <- quantile(treat3, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat3 <- q_treat3["95%"] / q_treat3["75%"]

C95_treat3

#study 4
data4 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Fig4D-F.xlsx",
  sheet = "Fig.4F"
)
#control group
control4 <- data4$`AD + VEH`

q_control4 <- quantile(control4, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control4 <- q_control4["95%"] / q_control4["75%"]

C95_control4

#treatment group
treat4 <- data4$`AD+DDL`

q_treat4 <- quantile(treat4, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat4 <- q_treat4["95%"] / q_treat4["75%"]

C95_treat4

#study 5
data5 <- read_excel(
  "PST-001_GMR-rough_eye_data-1.xlsx",
  sheet = "Sheet2(own sheet I added)"
)
#control group
control5 <- data5$Control

q_control5 <- quantile(control5, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control5 <- q_control5["95%"] / q_control5["75%"]

C95_control5

#treatment group
treat5 <- data5$Tau

q_treat5 <- quantile(treat5, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat5 <- q_treat5["95%"] / q_treat5["75%"]

C95_treat5

#study 6
data6 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Davidowitz_et_al_2023_Fig_1-8_data.xlsx",
  sheet = "Fig. 8"
)
#control group
control6 <- data6$Veh

q_control6 <- quantile(control6, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control6 <- q_control6["95%"] / q_control6["75%"]

C95_control6

#treatment group
treat6 <-  data6$`40`

q_treat6 <- quantile(treat6, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat6 <- q_treat6["95%"] / q_treat6["75%"]

C95_treat6

#study 7
data7 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Crenshaw+et+al+data+file.xlsx",
  sheet = "Figure 4B"
)
#control group
control7 <- data7$wt

q_control7 <- quantile(control7, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control7 <- q_control7["95%"] / q_control7["75%"]

C95_control7

#treatment group
treat7 <- data7$mut

q_treat7 <- quantile(treat7, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat7 <- q_treat7["95%"] / q_treat7["75%"]

C95_treat7

#study 8
data8 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/SourceData_Fig_2.xlsx",
  sheet = "Fig2C_polarization"
)
#control group
control8 <- data8$'5XFAD'

q_control8 <- quantile(control8, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control8 <- q_control8["95%"] / q_control8["75%"]

C95_control8

#treatment group
treat8 <- data8$'5XFAD;cKO'

q_treat8 <- quantile(treat8, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat8 <- q_treat8["95%"] / q_treat8["75%"]

C95_treat8

#study 9
data9 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Davidowitz_et_al._2025_Data_for_DRYAD_for_JNC_revised_07-22-2025.xlsx",
  sheet = "Fig. 2"
)
#control group
control9 <- data9$Veh

q_control9 <- quantile(control9, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control9 <- q_control9["95%"] / q_control9["75%"]

C95_control9

#treatment group
treat9 <- data9$'80 mg/kg'

q_treat9 <- quantile(treat9, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat9 <- q_treat9["95%"] / q_treat9["75%"]

C95_treat9


#study 10
data10 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Data_PONE-D-10-14787.xlsx",
  sheet = "Vanita's own sheet"
)
#control group
control10 <- data10$'TBI+Vehicle'

q_control10 <- quantile(control10, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control10 <- q_control10["95%"] / q_control10["75%"]

C95_control10

#treatment group
treat10 <- data10$'TBI+ AVL-3288'

q_treat10 <- quantile(treat10, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat10 <- q_treat10["95%"] / q_treat10["75%"]

C95_treat10

#study 11
data11 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure+4.xlsx",
  sheet = "Fig 4a"
)
#control group
control11 <- data11$control

q_control11 <- quantile(control11, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control11 <- q_control11["95%"] / q_control11["75%"]

C95_control11

#treatment group
treat11 <- data11$'2.5mM metformin'

q_treat11 <- quantile(treat11, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat11 <- q_treat11["95%"] / q_treat11["75%"]

C95_treat11

#study 12
data12 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure4_FSCV_Data_Analysis_edited2.xlsx",
  sheet = "Sheet 2"
)
#control group
control12 <- data12$'Vehicle/PFF'

q_control12 <- quantile(control12, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control12 <- q_control12["95%"] / q_control12["75%"]

C95_control12

#treatment group
treat12 <- data12$'Dieldrin/PFF'

q_treat12 <- quantile(treat12, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat12 <- q_treat12["95%"] / q_treat12["75%"]

C95_treat12

#study 13
data13 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figures_1_and_2_Dataset_.xlsx",
  sheet = "Figure 1B"
)
#control group
control13 <- data13$WT

q_control13 <- quantile(control13, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control13 <- q_control13["95%"] / q_control13["75%"]

C95_control13

#treatment group
treat13 <- data13$'OGT-1;EEL-1'

q_treat13 <- quantile(treat13, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat13 <- q_treat13["95%"] / q_treat13["75%"]

C95_treat13

#study 14
data14 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/elife-57155-Figure 2-source data 1.xlsx",
  sheet = "Figure 2E"
)
#control group
control14 <-  data14$'Ctrl (% Freezing )'

q_control14 <- quantile(control14, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control14 <- q_control14["95%"] / q_control14["75%"]

C95_control14

#treatment group
treat14 <- data14$'Stim (% Freezing )'

q_treat14 <- quantile(treat14, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat14 <- q_treat14["95%"] / q_treat14["75%"]

C95_treat14

#study 15
data15 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Paternal_atrazine_exposure_on_F1_zebrafish_behavior.xlsx",
  sheet = "Vanita's Sheet"
)
#control group
control15 <- data15$Control2

q_control15 <- quantile(control15, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control15 <- q_control15["95%"] / q_control15["75%"]

C95_control15

#treatment group
treat15 <- data15$AZT0.32

q_treat15 <- quantile(treat15, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat15 <- q_treat15["95%"] / q_treat15["75%"]

C95_treat15

#study 16
data16 <- read.csv(
  "/Users/vanita/Downloads/StatStudy/Behavior_Data.csv"
)

#control group
control16 <- data16$DistanceEPM[data16$Genotype == "C57"]

q_control16 <- quantile(control16, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control16 <- q_control16["95%"] / q_control16["75%"]

C95_control16

#treatment group
treat16 <- data16$DistanceEPM[data16$Genotype == "IL4KO"]

q_treat16 <- quantile(treat16, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat16 <- q_treat16["95%"] / q_treat16["75%"]

C95_treat16

#study 17
data17 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/SourceData_Fig_2 copy (study17).xlsx",
  sheet = "Fig2D_EPM"
)
#control group
control17 <-  data17$Control

q_control17 <- quantile(control17, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control17 <- q_control17["95%"] / q_control17["75%"]

C95_control17

#treatment group
treat17 <- data17$'PERK cKO'

q_treat17 <- quantile(treat17, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat17 <- q_treat17["95%"] / q_treat17["75%"]

C95_treat17

#study 18
data18 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/NCOMMS-20-01106_data_source_file.xlsx",
  sheet = "Fig. 3"
)
#control group
control18 <-  data18$wt

q_control18 <- quantile(control18, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control18 <- q_control18["95%"] / q_control18["75%"]

C95_control18

#treatment group
treat18 <- data18$'mut'

q_treat18 <- quantile(treat18, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat18 <- q_treat18["95%"] / q_treat18["75%"]

C95_treat18

#study 19
data19 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure4B.xlsx",
  sheet = "Nile Red Relative to Control"
)
#control group
control19 <-  data19$Control

q_control19 <- quantile(control19, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control19 <- q_control19["95%"] / q_control19["75%"]

C95_control19

#treatment group
treat19 <- data19$Upregulation

q_treat19 <- quantile(treat19, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat19 <- q_treat19["95%"] / q_treat19["75%"]

C95_treat19

#study 20
data20 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure_7_raw_data copy.xlsx",
  sheet = "Vanita's Own Sheet"
)
#control group
control20 <-  data20$Control

q_control20 <- quantile(control20, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control20 <- q_control20["95%"] / q_control20["75%"]

C95_control20

#treatment group
treat20 <- data20$'B. pertussis'

q_treat20 <- quantile(treat20, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat20 <- q_treat20["95%"] / q_treat20["75%"]

C95_treat20

#study 21
data21 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure5_Spontaneous Repetitive Behaviors.xlsx",
  sheet = "Vanita's Own Sheet"
)
#control group
control21 <-  data21$WT

q_control21 <- quantile(control21, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control21 <- q_control21["95%"] / q_control21["75%"]

C95_control21

#treatment group
treat21 <- data21$'Mut'

q_treat21 <- quantile(treat21, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat21 <- q_treat21["95%"] / q_treat21["75%"]

C95_treat21

#study 22
data22 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/sucrose_preference_context_group.xlsx",
  sheet = "Sheet2"
)
#control group
control22 <-  data22$'uninjured sham, not shocked'

q_control22 <- quantile(control22, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control22 <- q_control22["95%"] / q_control22["75%"]

C95_control22

#treatment group
treat22 <- data22$'rcTBI, not shocked'

q_treat22 <- quantile(treat22, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat22 <- q_treat22["95%"] / q_treat22["75%"]

C95_treat22

#study 23
data23 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Raw+data+for+dyrad.xlsx",
  sheet = "Vanita's Own Sheet"
)
#control group
control23 <-  data23$'WT mice'

q_control23 <- quantile(control23, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control23 <- q_control23["95%"] / q_control23["75%"]

C95_control23

#treatment group
treat23 <- data23$'miR-137 Tg mice'

q_treat23 <- quantile(treat23, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat23 <- q_treat23["95%"] / q_treat23["75%"]

C95_treat23

#study 24
data24 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Bix+supplemental+data.xlsx",
  sheet = "Vanita's Own Sheet"
)
#control group
control24 <-  data24$Saline

q_control24 <- quantile(control24, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control24 <- q_control24["95%"] / q_control24["75%"]

C95_control24

#treatment group
treat24 <- data24$'BIX injection'

q_treat24 <- quantile(treat24, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat24 <- q_treat24["95%"] / q_treat24["75%"]

C95_treat24

#study 25
data25 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Guillot_2025_MCHdata.xlsx"
)
#control group
control25 <-  data25$Baseline

q_control25 <- quantile(control25, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control25 <- q_control25["95%"] / q_control25["75%"]

C95_control25

#treatment group
treat25 <- data25$MCH

q_treat25 <- quantile(treat25, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat25 <- q_treat25["95%"] / q_treat25["75%"]

C95_treat25

#study 26
data26 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Different_intensity_exercise_and_verbal_spatial_working_memory_Data.xlsx"
)
#control group
control26 <-  data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 4]

q_control26 <- quantile(control26, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control26 <- q_control26["95%"] / q_control26["75%"]

C95_control26

#treatment group
treat26 <- data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 2]

q_treat26 <- quantile(treat26, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat26 <- q_treat26["95%"] / q_treat26["75%"]

C95_treat26

#study 27
data27 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/EXP4_RM.xlsx"
)
#control group
control27 <- data27$'POST_iTBS+sham-γtACS'

q_control27 <- quantile(control27, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control27 <- q_control27["95%"] / q_control27["75%"]

C95_control27

#treatment group
treat27 <- data27$'POST_iTBS+γtACS'

q_treat27 <- quantile(treat27, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat27 <- q_treat27["95%"] / q_treat27["75%"]

C95_treat27

#study 28
data28 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Wellmind_ms_data.xlsx"
)
#control group
control28 <- data28$'Post- selfcompassion'[data28$'Group' == 'Control']

q_control28 <- quantile(control28, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control28 <- q_control28["95%"] / q_control28["75%"]

C95_control28

#treatment group
treat28 <- data28$'Post- selfcompassion'[data28$'Group' == 'WellMind']

q_treat28 <- quantile(treat28, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat28 <- q_treat28["95%"] / q_treat28["75%"]

C95_treat28

#study 29
data29 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/MSTP_ms_Data_Dryad.xlsx"
)
#control group
control29 <-  na.omit(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'B'])

q_control29 <- quantile(control29, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control29 <- q_control29["95%"] / q_control29["75%"]

C95_control29

#treatment group
treat29 <-   na.omit(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'A'])

q_treat29 <- quantile(treat29, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat29 <- q_treat29["95%"] / q_treat29["75%"]

C95_treat29

#study 30
data30 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/characteristics_behavioral_data.xlsx"
)
#control group
control30 <-  data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'waiting group']

q_control30 <- quantile(control30, probs = c(0.75, 0.95), na.rm = TRUE)

C95_control30 <- q_control30["95%"] / q_control30["75%"]

C95_control30

#treatment group
treat30 <-    data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'diet group']

q_treat30 <- quantile(treat30, probs = c(0.75, 0.95), na.rm = TRUE)

C95_treat30 <- q_treat30["95%"] / q_treat30["75%"]

C95_treat30

#combine into vectors
C_control <- c(
  C95_control1, C95_control2, C95_control3, C95_control4, C95_control5,
  C95_control6, C95_control7, C95_control8, C95_control9, C95_control10,
  C95_control11, C95_control12, C95_control13, C95_control14, C95_control15,
  C95_control16, C95_control17, C95_control18, C95_control19, C95_control20,
  C95_control21, C95_control22, C95_control23, C95_control24, C95_control25,
  C95_control26, C95_control27, C95_control28, C95_control29, C95_control30
)

C_treat <- c(
  C95_treat1, C95_treat2, C95_treat3, C95_treat4, C95_treat5,
  C95_treat6, C95_treat7, C95_treat8, C95_treat9, C95_treat10,
  C95_treat11, C95_treat12, C95_treat13, C95_treat14, C95_treat15,
  C95_treat16, C95_treat17, C95_treat18, C95_treat19, C95_treat20,
  C95_treat21, C95_treat22, C95_treat23, C95_treat24, C95_treat25,
  C95_treat26, C95_treat27, C95_treat28, C95_treat29, C95_treat30
)

#compute summary statistics
summary_C <- tibble(
  Group = c("Control", "Treatment"),
  Mean = c(mean(C_control, na.rm = TRUE),
           mean(C_treat, na.rm = TRUE)),
  Median = c(median(C_control, na.rm = TRUE),
             median(C_treat, na.rm = TRUE)),
  Q1 = c(quantile(C_control, 0.25, na.rm = TRUE),
         quantile(C_treat, 0.25, na.rm = TRUE)),
  Q3 = c(quantile(C_control, 0.75, na.rm = TRUE),
         quantile(C_treat, 0.75, na.rm = TRUE))
)
summary_C

## KURTOSIS ##

#STUDY 1
data1 <- read_excel("/Users/vanita/Downloads/StatStudy/APEx_FinalData_20201021 (1).xlsx")

#skew/kurtosis - CONTROL
k_control1 <- Kurt(
  data1$amyloid_Mean6_CB[data1$apx_group == "NoSupport"],
  conf.level = 0.95,
  na.rm = TRUE
)
k_control1
#skew/kurtosis - Treatment
k_treat1 <- Kurt(
  data1$amyloid_Mean6_CB[data1$apx_group == "AExSupport"],
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat1

# number in each group
n_control1 <- length(na.omit(data1$amyloid_Mean6_CB[data1$apx_group == "NoSupport"]))
n_treat1   <- length(na.omit(data1$amyloid_Mean6_CB[data1$apx_group == "AExSupport"]))
n_control1
n_treat1

#variance ratio
var_control1 <- var(data1$amyloid_Mean6_CB[data1$apx_group == "NoSupport"], na.rm=TRUE)
var_treat1   <- var(data1$amyloid_Mean6_CB[data1$apx_group == "AExSupport"], na.rm=TRUE)

var_ratio1 <- var_treat1 / var_control1
var_ratio1

study1_summary <- data.frame(
  study = "data1",
  kurt_control = k_control1[1],
  lwr_kurt_control = k_control1[2],
  upr_kurt_control = k_control1[3],
  kurt_treat = k_treat1[1],
  lwr_kurt_treat = k_treat1[2],
  upr_kurt_treat = k_treat1[3],
  var_ratio = var_ratio1,
  n_control = n_control1,
  n_treat = n_treat1
)
study1_summary



#STUDY 2
data2 <-  read_excel("/Users/vanita/Downloads/StatStudy/Cucos et al., MYD88_mRNAlevels_Mice.xlsx")
names(data2)
#skew/kurtosis - CONTROL
k_control2 <- Kurt(
  data2$`MYD88_T0 (baseline)`[data2$Case == "WT"],
  conf.level = 0.95,
  na.rm = TRUE
)
k_control2
#skew/kurtosis - Treatment
k_treat2 <- Kurt(
  data2$`MYD88_T0 (baseline)`[data2$Case == "AT_transgenic"],
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat2

# number in each group
n_control2 <- length(na.omit(data2$`MYD88_T0 (baseline)`[data2$Case == "WT"]))
n_treat2   <- length(na.omit(data2$`MYD88_T0 (baseline)`[data2$Case == "AT_transgenic"]))
n_control2
n_treat2

#variance ratio
var_control2 <- var(data2$`MYD88_T0 (baseline)`[data2$Case == "WT"], na.rm=TRUE)
var_treat2  <- var(data2$`MYD88_T0 (baseline)`[data2$Case == "AT_transgenic"], na.rm=TRUE)

var_ratio2 <- var_treat2 / var_control2
var_ratio2

study2_summary <- data.frame(
  study = "data2",
  kurt_control = k_control2[1],
  lwr_kurt_control = k_control2[2],
  upr_kurt_control = k_control2[3],
  kurt_treat = k_treat2[1],
  lwr_kurt_treat = k_treat2[2],
  upr_kurt_treat = k_treat2[3],
  var_ratio = var_ratio2,
  n_control = n_control2,
  n_treat = n_treat2
)
study2_summary


#STUDY 3
data3 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/DATA _ SUBMISSION 10.27.  25 - .xlsx",
  sheet = "y MAZE "
)

names(data3)

#skew/kurtosis - CONTROL
k_control3 <- Kurt(
  data3$'D/A + Veh',
  conf.level = 0.95,
  na.rm = TRUE
)
k_control3
#skew/kurtosis - Treatment
k_treat3 <- Kurt(
  data3$'D/A +INS',
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat3

# number in each group
n_control3 <- length(na.omit(data3$Ctrl))
n_treat3  <- length(na.omit(data3$`D/A +INS`))
n_control3
n_treat3

#variance ratio
var_control3 <- var(data3$Ctrl, na.rm=TRUE)
var_treat3 <- var(data3$`D/A +INS`, na.rm=TRUE)

var_ratio3 <- var_treat3 / var_control3
var_ratio3

study3_summary <- data.frame(
  study = "data3",
  kurt_control = k_control3[1],
  lwr_kurt_control = k_control3[2],
  upr_kurt_control = k_control3[3],
  kurt_treat = k_treat3[1],
  lwr_kurt_treat = k_treat3[2],
  upr_kurt_treat = k_treat3[3],
  var_ratio = var_ratio3,
  n_control = n_control3,
  n_treat = n_treat3
)
study3_summary

#STUDY 4
data4 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Fig4D-F.xlsx",
  sheet = "Fig.4F"
)

names(data4)

#skew/kurtosis - CONTROL
k_control4 <- Kurt(
  data4$`AD + VEH`,
  conf.level = 0.95,
  na.rm = TRUE
)
k_control4
#skew/kurtosis - Treatment
k_treat4 <- Kurt(
  data4$`AD+DDL`,
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat4

# number in each group
n_control4 <- length(na.omit( data4$`AD + VEH`))
n_treat4  <- length(na.omit(data4$`AD+DDL`))
n_control4
n_treat4

#variance ratio
var_control4 <- var(data4$`AD + VEH`, na.rm=TRUE)
var_treat4 <- var(data4$`AD+DDL`, na.rm=TRUE)

var_ratio4 <- var_treat4 / var_control4
var_ratio4

study4_summary <- data.frame(
  study = "data4",
  kurt_control = k_control4[1],
  lwr_kurt_control = k_control4[2],
  upr_kurt_control = k_control4[3],
  kurt_treat = k_treat4[1],
  lwr_kurt_treat = k_treat4[2],
  upr_kurt_treat = k_treat4[3],
  var_ratio = var_ratio4,
  n_control = n_control4,
  n_treat = n_treat4
)
study4_summary

#STUDY 5
data5 <- read_excel(
  "PST-001_GMR-rough_eye_data-1.xlsx",
  sheet = "Sheet2(own sheet I added)"
)

names(data5)

#skew/kurtosis - CONTROL
k_control5 <- Kurt(
  data5$Control,
  conf.level = 0.95,
  na.rm = TRUE
)
k_control5
#skew/kurtosis - Treatment
k_treat5 <- Kurt(
  data5$Tau,
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat5

# number in each group
n_control5 <- length(na.omit( data5$Control))
n_treat5  <- length(na.omit(data5$Tau))
n_control5
n_treat5

#variance ratio
var_control5 <- var(data5$Control, na.rm=TRUE)
var_treat5 <- var(data5$Tau, na.rm=TRUE)

var_ratio5 <- var_treat5 / var_control5
var_ratio5

study5_summary <- data.frame(
  study = "data5",
  kurt_control = k_control5[1],
  lwr_kurt_control = k_control5[2],
  upr_kurt_control = k_control5[3],
  kurt_treat = k_treat5[1],
  lwr_kurt_treat = k_treat5[2],
  upr_kurt_treat = k_treat5[3],
  var_ratio = var_ratio5,
  n_control = n_control5,
  n_treat = n_treat5
)
study5_summary

#STUDY 6
data6 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Davidowitz_et_al_2023_Fig_1-8_data.xlsx",
  sheet = "Fig. 8"
)
names(data6)

#skew/kurtosis - CONTROL
k_control6 <- Kurt(
  data6$Veh,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control6
#skew/kurtosis - Treatment
k_treat6 <- Kurt(
  data6$`40`,      # backticks around the column name
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat6
# number in each group
n_control6 <- length(na.omit(data6$Veh))
n_treat6  <- length(na.omit( data6$`40`))
n_control6
n_treat6

#variance ratio
var_control6 <- var(data6$Veh, na.rm=TRUE)
var_treat6 <- var(data6$`40`, na.rm=TRUE)

var_ratio6 <- var_treat6 / var_control6
var_ratio6

study6_summary <- data.frame(
  study = "data6",
  kurt_control = k_control6[1],
  lwr_kurt_control = k_control6[2],
  upr_kurt_control = k_control6[3],
  kurt_treat = k_treat6[1],
  lwr_kurt_treat = k_treat6[2],
  upr_kurt_treat = k_treat6[3],
  var_ratio = var_ratio6,
  n_control = n_control6,
  n_treat = n_treat6
)
study6_summary

#STUDY 7
data7 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Crenshaw+et+al+data+file.xlsx",
  sheet = "Figure 4B"
)

# Check column names
colnames(data7)
head(data7)  # check the first few rows

#skew/kurtosis - CONTROL
k_control7 <- Kurt(
  data7$wt,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control7
#skew/kurtosis - Treatment
k_treat7 <- Kurt(
  data7$mut,      # backticks around the column name
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat7
# number in each group
n_control7 <- length(na.omit(data7$wt))
n_treat7  <- length(na.omit(data7$mut))
n_control7
n_treat7

#variance ratio
var_control7 <- var(data7$wt, na.rm=TRUE)
var_treat7 <- var(data7$mut, na.rm=TRUE)

var_ratio7 <- var_treat7 / var_control7
var_ratio7

study7_summary <- data.frame(
  study = "data7",
  kurt_control = k_control7[1],
  lwr_kurt_control = k_control7[2],
  upr_kurt_control = k_control7[3],
  kurt_treat = k_treat7[1],
  lwr_kurt_treat = k_treat7[2],
  upr_kurt_treat = k_treat7[3],
  var_ratio = var_ratio7,
  n_control = n_control7,
  n_treat = n_treat7
)
study7_summary

#STUDY 8
data8 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/SourceData_Fig_2.xlsx",
  sheet = "Fig2C_polarization"
)
names(data8)

#skew/kurtosis - CONTROL
k_control8 <- Kurt(
  data8$'5XFAD',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control8
#skew/kurtosis - Treatment
k_treat8 <- Kurt(
  data8$'5XFAD;cKO',
  conf.level = 0.95,
  na.rm = TRUE
)

# number in each group
n_control8 <- length(na.omit( data8$'5XFAD'))
n_treat8  <- length(na.omit(  data8$'5XFAD;cKO'))
n_control8
n_treat8

#variance ratio
var_control8 <- var(data8$'5XFAD', na.rm=TRUE)
var_treat8 <- var(data8$'5XFAD;cKO', na.rm=TRUE)

var_ratio8 <- var_treat8 / var_control8
var_ratio8

study8_summary <- data.frame(
  study = "data8",
  kurt_control = k_control8[1],
  lwr_kurt_control = k_control8[2],
  upr_kurt_control = k_control8[3],
  kurt_treat = k_treat8[1],
  lwr_kurt_treat = k_treat8[2],
  upr_kurt_treat = k_treat8[3],
  var_ratio = var_ratio8,
  n_control = n_control8,
  n_treat = n_treat8
)
study8_summary

#STUDY 9
data9 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Davidowitz_et_al._2025_Data_for_DRYAD_for_JNC_revised_07-22-2025.xlsx",
  sheet = "Fig. 2"
)
names(data9)

#skew/kurtosis - CONTROL
k_control9 <- Kurt(
  data9$Veh,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control9
#skew/kurtosis - Treatment
k_treat9 <- Kurt(
  data9$'80 mg/kg',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat9
# number in each group
n_control9 <- length(na.omit( data9$Veh))
n_treat9  <- length(na.omit( data9$'80 mg/kg'))
n_control9
n_treat9

#variance ratio
var_control9 <- var(data9$Veh, na.rm=TRUE)
var_treat9 <- var(data9$'80 mg/kg', na.rm=TRUE)

var_ratio9 <- var_treat9 / var_control9
var_ratio9

study9_summary <- data.frame(
  study = "data9",
  kurt_control = k_control9[1],
  lwr_kurt_control = k_control9[2],
  upr_kurt_control = k_control9[3],
  kurt_treat = k_treat9[1],
  lwr_kurt_treat = k_treat9[2],
  upr_kurt_treat = k_treat9[3],
  var_ratio = var_ratio9,
  n_control = n_control9,
  n_treat = n_treat9
)
study9_summary

#STUDY 10
data10 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Data_PONE-D-10-14787.xlsx",
  sheet = "Vanita's own sheet"
)
names(data10)

#skew/kurtosis - CONTROL
k_control10 <- Kurt(
  data10$'TBI+Vehicle',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control10
#skew/kurtosis - Treatment
k_treat10 <- Kurt(
  data10$'TBI+ AVL-3288',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat10
# number in each group
n_control10 <- length(na.omit(data10$'TBI+Vehicle'))
n_treat10  <- length(na.omit(data10$'TBI+ AVL-3288'))
n_control10
n_treat10

#variance ratio
var_control10 <- var(data10$'TBI+Vehicle', na.rm=TRUE)
var_treat10 <- var(data10$'TBI+ AVL-3288', na.rm=TRUE)

var_ratio10 <- var_treat10 / var_control10
var_ratio10

study10_summary <- data.frame(
  study = "data10",
  kurt_control = k_control10[1],
  lwr_kurt_control = k_control10[2],
  upr_kurt_control = k_control10[3],
  kurt_treat = k_treat10[1],
  lwr_kurt_treat = k_treat10[2],
  upr_kurt_treat = k_treat10[3],
  var_ratio = var_ratio10,
  n_control = n_control10,
  n_treat = n_treat10
)
study10_summary

#STUDY 11
data11 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure+4.xlsx",
  sheet = "Fig 4a"
)
names(data11)

#skew/kurtosis - CONTROL
k_control11 <- Kurt(
  data11$control,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control11
#skew/kurtosis - Treatment
k_treat11 <- Kurt(
  data11$'2.5mM metformin',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat11
# number in each group
n_control11 <- length(na.omit(  data11$control))
n_treat11  <- length(na.omit(data11$'2.5mM metformin'))
n_control11
n_treat11

#variance ratio
var_control11 <- var(data11$control, na.rm=TRUE)
var_treat11 <- var(data11$'2.5mM metformin', na.rm=TRUE)

var_ratio11 <- var_treat11 / var_control11
var_ratio11

study11_summary <- data.frame(
  study = "data11",
  kurt_control = k_control11[1],
  lwr_kurt_control = k_control11[2],
  upr_kurt_control = k_control11[3],
  kurt_treat = k_treat11[1],
  lwr_kurt_treat = k_treat11[2],
  upr_kurt_treat = k_treat11[3],
  var_ratio = var_ratio11,
  n_control = n_control11,
  n_treat = n_treat11
)
study11_summary

#STUDY 12
data12 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure4_FSCV_Data_Analysis_edited2.xlsx",
  sheet = "Sheet 2"
)
names(data12)

#skew/kurtosis - CONTROL
k_control12 <- Kurt(
  data12$'Vehicle/PFF',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control12
#skew/kurtosis - Treatment
k_treat12 <- Kurt(
  data12$'Dieldrin/PFF',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat12
# number in each group
n_control12 <- length(na.omit(  data12$'Vehicle/PFF'))
n_treat12  <- length(na.omit(data12$'Dieldrin/PFF'))
n_control12
n_treat12

#variance ratio
var_control12 <- var(data12$'Vehicle/PFF', na.rm=TRUE)
var_treat12 <- var(data12$'Dieldrin/PFF', na.rm=TRUE)

var_ratio12 <- var_treat12 / var_control12
var_ratio12

study12_summary <- data.frame(
  study = "data12",
  kurt_control = k_control12[1],
  lwr_kurt_control = k_control12[2],
  upr_kurt_control = k_control12[3],
  kurt_treat = k_treat12[1],
  lwr_kurt_treat = k_treat12[2],
  upr_kurt_treat = k_treat12[3],
  var_ratio = var_ratio12,
  n_control = n_control12,
  n_treat = n_treat12
)
study12_summary

#STUDY 13
data13 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figures_1_and_2_Dataset_.xlsx",
  sheet = "Figure 1B"
)
names(data13)

#skew/kurtosis - CONTROL
k_control13 <- Kurt(
  data13$WT,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control13
#skew/kurtosis - Treatment
k_treat13 <- Kurt(
  data13$'OGT-1;EEL-1',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat13
# number in each group
n_control13 <- length(na.omit(data13$WT))
n_treat13  <- length(na.omit(data13$'OGT-1;EEL-1'))
n_control13
n_treat13

#variance ratio
var_control13 <- var(data13$WT, na.rm=TRUE)
var_treat13 <- var(data13$'OGT-1;EEL-1', na.rm=TRUE)

var_ratio13 <- var_treat13 / var_control13
var_ratio13

study13_summary <- data.frame(
  study = "data13",
  kurt_control = k_control13[1],
  lwr_kurt_control = k_control13[2],
  upr_kurt_control = k_control13[3],
  kurt_treat = k_treat13[1],
  lwr_kurt_treat = k_treat13[2],
  upr_kurt_treat = k_treat13[3],
  var_ratio = var_ratio13,
  n_control = n_control13,
  n_treat = n_treat13
)
study13_summary

#STUDY 14
data14 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/elife-57155-Figure 2-source data 1.xlsx",
  sheet = "Figure 2E"
)
names(data14)

#skew/kurtosis - CONTROL
k_control14 <- Kurt(
  data14$'Ctrl (% Freezing )',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control14
#skew/kurtosis - Treatment
k_treat14 <- Kurt(
  data14$'Stim (% Freezing )',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat14
# number in each group
n_control14 <- length(na.omit(data14$'Ctrl (% Freezing )'))
n_treat14  <- length(na.omit(data14$'Stim (% Freezing )'))
n_control14
n_treat14

#variance ratio
var_control14 <- var(data14$'Ctrl (% Freezing )', na.rm=TRUE)
var_treat14 <- var(data14$'Stim (% Freezing )', na.rm=TRUE)

var_ratio14 <- var_treat14 / var_control14
var_ratio14

study14_summary <- data.frame(
  study = "data14",
  kurt_control = k_control14[1],
  lwr_kurt_control = k_control14[2],
  upr_kurt_control = k_control14[3],
  kurt_treat = k_treat14[1],
  lwr_kurt_treat = k_treat14[2],
  upr_kurt_treat = k_treat14[3],
  var_ratio = var_ratio14,
  n_control = n_control14,
  n_treat = n_treat14
)
study14_summary

#STUDY 15
data15 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Paternal_atrazine_exposure_on_F1_zebrafish_behavior.xlsx",
  sheet = "Vanita's Sheet"
)
names(data15)

#skew/kurtosis - CONTROL
k_control15 <- Kurt(
  data15$Control2,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control15
#skew/kurtosis - Treatment
k_treat15 <- Kurt(
  data15$AZT0.32,
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat15
# number in each group
n_control15 <- length(na.omit(data15$Control2))
n_treat15  <- length(na.omit( data15$AZT0.32))
n_control15
n_treat15

#variance ratio
var_control15 <- var(data15$Control2, na.rm=TRUE)
var_treat15<- var(data15$AZT0.32, na.rm=TRUE)

var_ratio15 <- var_treat15 / var_control15
var_ratio15

study15_summary <- data.frame(
  study = "data15",
  kurt_control = k_control15[1],
  lwr_kurt_control = k_control15[2],
  upr_kurt_control = k_control15[3],
  kurt_treat = k_treat15[1],
  lwr_kurt_treat = k_treat15[2],
  upr_kurt_treat = k_treat15[3],
  var_ratio = var_ratio15,
  n_control = n_control15,
  n_treat = n_treat15
)
study15_summary

#STUDY 16
data16 <- read.csv(
  "/Users/vanita/Downloads/StatStudy/Behavior_Data.csv"
)

names(data16)

control16 <- data16$DistanceEPM[data16$Genotype == "C57"]
treat16   <- data16$DistanceEPM[data16$Genotype == "IL4KO"]

# CONTROL (C57)
k_control16 <- Kurt(
  control16,
  conf.level = 0.95,
  na.rm = TRUE
)
k_control16

# TREATMENT (IL4KO)
k_treat16 <- Kurt(
  treat16,
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat16

# number in each group
n_control16 <- length(na.omit(control16))
n_treat16   <- length(na.omit(treat16))
n_control16
n_treat16

#variance ratio
var_control16 <- var(control16, na.rm=TRUE)
var_treat16 <- var(treat16, na.rm=TRUE)

var_ratio16 <- var_treat16 / var_control16
var_ratio16

study16_summary <- data.frame(
  study = "data16",
  kurt_control = k_control16[1],
  lwr_kurt_control = k_control16[2],
  upr_kurt_control = k_control16[3],
  kurt_treat = k_treat16[1],
  lwr_kurt_treat = k_treat16[2],
  upr_kurt_treat = k_treat16[3],
  var_ratio = var_ratio16,
  n_control = n_control16,
  n_treat = n_treat16
)
study16_summary

#STUDY 17
data17 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/SourceData_Fig_2 copy (study17).xlsx",
  sheet = "Fig2D_EPM"
)
names(data17)

#skew/kurtosis - CONTROL
k_control17 <- Kurt(
  data17$Control,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control17
#skew/kurtosis - Treatment
k_treat17 <- Kurt(
  data17$'PERK cKO',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat17
# number in each group
n_control17 <- length(na.omit(data17$Control))
n_treat17  <- length(na.omit(data17$'PERK cKO'))
n_control17
n_treat17

#variance ratio
var_control17 <- var(data17$Control, na.rm=TRUE)
var_treat17 <- var(data17$'PERK cKO', na.rm=TRUE)

var_ratio17 <- var_treat17 / var_control17
var_ratio17

study17_summary <- data.frame(
  study = "data17",
  kurt_control = k_control17[1],
  lwr_kurt_control = k_control17[2],
  upr_kurt_control = k_control17[3],
  kurt_treat = k_treat17[1],
  lwr_kurt_treat = k_treat17[2],
  upr_kurt_treat = k_treat17[3],
  var_ratio = var_ratio17,
  n_control = n_control17,
  n_treat = n_treat17
)
study17_summary

#STUDY 18
data18 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/NCOMMS-20-01106_data_source_file.xlsx",
  sheet = "Fig. 3"
)
names(data18)

#skew/kurtosis - CONTROL
k_control18 <- Kurt(
  data18$wt,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control18
#skew/kurtosis - Treatment
k_treat18 <- Kurt(
  data18$'mut',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat18
# number in each group
n_control18 <- length(na.omit(data18$wt))
n_treat18  <- length(na.omit(data18$'mut'))
n_control18
n_treat18

#variance ratio
var_control18 <- var(data18$wt, na.rm=TRUE)
var_treat18 <- var(data18$'mut', na.rm=TRUE)

var_ratio18 <- var_treat18 / var_control18
var_ratio18

study18_summary <- data.frame(
  study = "data18",
  kurt_control = k_control18[1],
  lwr_kurt_control = k_control18[2],
  upr_kurt_control = k_control18[3],
  kurt_treat = k_treat18[1],
  lwr_kurt_treat = k_treat18[2],
  upr_kurt_treat = k_treat18[3],
  var_ratio = var_ratio18,
  n_control = n_control18,
  n_treat = n_treat18
)
study18_summary

#STUDY 19
data19 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure4B.xlsx",
  sheet = "Nile Red Relative to Control"
)
names(data19)

#skew/kurtosis - CONTROL
k_control19 <- Kurt(
  data19$Control,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control19
#skew/kurtosis - Treatment
k_treat19 <- Kurt(
  data19$Upregulation,
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat19
# number in each group
n_control19 <- length(na.omit(data19$Control))
n_treat19  <- length(na.omit(data19$Upregulation))
n_control19
n_treat19

#variance ratio
var_control19 <- var(data19$Control, na.rm=TRUE)
var_treat19 <- var(data19$Upregulation, na.rm=TRUE)

var_ratio19 <- var_treat19 / var_control19
var_ratio19

study19_summary <- data.frame(
  study = "data19",
  kurt_control = k_control19[1],
  lwr_kurt_control = k_control19[2],
  upr_kurt_control = k_control19[3],
  kurt_treat = k_treat19[1],
  lwr_kurt_treat = k_treat19[2],
  upr_kurt_treat = k_treat19[3],
  var_ratio = var_ratio19,
  n_control = n_control19,
  n_treat = n_treat19
)
study19_summary

#STUDY 20
data20 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure_7_raw_data copy.xlsx",
  sheet = "Vanita's Own Sheet"
)
names(data20)

#skew/kurtosis - CONTROL
k_control20 <- Kurt(
  data20$Control,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control20
#skew/kurtosis - Treatment
k_treat20 <- Kurt(
  data20$'B. pertussis',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat20
# number in each group
n_control20 <- length(na.omit( data20$Control))
n_treat20 <- length(na.omit(data20$'B. pertussis'))
n_control20
n_treat20

#variance ratio
var_control20 <- var(data20$Control, na.rm=TRUE)
var_treat20 <- var(data20$'B. pertussis', na.rm=TRUE)

var_ratio20 <- var_treat20 / var_control20
var_ratio20

study20_summary <- data.frame(
  study = "data20",
  kurt_control = k_control20[1],
  lwr_kurt_control = k_control20[2],
  upr_kurt_control = k_control20[3],
  kurt_treat = k_treat20[1],
  lwr_kurt_treat = k_treat20[2],
  upr_kurt_treat = k_treat20[3],
  var_ratio = var_ratio20,
  n_control = n_control20,
  n_treat = n_treat20
)
study20_summary

#STUDY 21
data21 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Figure5_Spontaneous Repetitive Behaviors.xlsx",
  sheet = "Vanita's Own Sheet"
)
names(data21)

#skew/kurtosis - CONTROL
k_control21 <- Kurt(
  data21$WT,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control21
#skew/kurtosis - Treatment
k_treat21 <- Kurt(
  data21$'Mut',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat21
# number in each group
n_control21 <- length(na.omit( data21$WT))
n_treat21 <- length(na.omit(data21$'Mut'))
n_control21
n_treat21

#variance ratio
var_control21 <- var(data21$WT, na.rm=TRUE)
var_treat21 <- var(data21$'Mut', na.rm=TRUE)

var_ratio21 <- var_treat21 / var_control21
var_ratio21

study21_summary <- data.frame(
  study = "data21",
  kurt_control = k_control21[1],
  lwr_kurt_control = k_control21[2],
  upr_kurt_control = k_control21[3],
  kurt_treat = k_treat21[1],
  lwr_kurt_treat = k_treat21[2],
  upr_kurt_treat = k_treat21[3],
  var_ratio = var_ratio21,
  n_control = n_control21,
  n_treat = n_treat21
)
study21_summary

#STUDY 22
data22 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/sucrose_preference_context_group.xlsx",
  sheet = "Sheet2"
)
names(data22)

#skew/kurtosis - CONTROL
k_control22 <- Kurt(
  data22$'uninjured sham, not shocked',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control22
#skew/kurtosis - Treatment
k_treat22 <- Kurt(
  data22$'rcTBI, not shocked',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat22
# number in each group
n_control22 <- length(na.omit(data22$'uninjured sham, not shocked'))
n_treat22 <- length(na.omit(data22$'rcTBI, not shocked'))
n_control22
n_treat22

#variance ratio
var_control22 <- var(data22$'uninjured sham, not shocked', na.rm=TRUE)
var_treat22 <- var(data22$'rcTBI, not shocked', na.rm=TRUE)

var_ratio22 <- var_treat22 / var_control22
var_ratio22

study22_summary <- data.frame(
  study = "data22",
  kurt_control = k_control22[1],
  lwr_kurt_control = k_control22[2],
  upr_kurt_control = k_control22[3],
  kurt_treat = k_treat22[1],
  lwr_kurt_treat = k_treat22[2],
  upr_kurt_treat = k_treat22[3],
  var_ratio = var_ratio22,
  n_control = n_control22,
  n_treat = n_treat22
)
study22_summary

#STUDY 23
data23 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Raw+data+for+dyrad.xlsx",
  sheet = "Vanita's Own Sheet"
)
names(data23)

#skew/kurtosis - CONTROL
k_control23 <- Kurt(
  data23$'WT mice',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control23
#skew/kurtosis - Treatment
k_treat23 <- Kurt(
  data23$'miR-137 Tg mice',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat23
# number in each group
n_control23 <- length(na.omit( data23$'WT mice'))
n_treat23 <- length(na.omit( data23$'miR-137 Tg mice'))
n_control23
n_treat23

#variance ratio
var_control23 <- var(data23$'WT mice', na.rm=TRUE)
var_treat23 <- var(data23$'miR-137 Tg mice', na.rm=TRUE)

var_ratio23 <- var_treat23 / var_control23
var_ratio23

study23_summary <- data.frame(
  study = "data23",
  kurt_control = k_control23[1],
  lwr_kurt_control = k_control23[2],
  upr_kurt_control = k_control23[3],
  kurt_treat = k_treat23[1],
  lwr_kurt_treat = k_treat23[2],
  upr_kurt_treat = k_treat23[3],
  var_ratio = var_ratio23,
  n_control = n_control23,
  n_treat = n_treat23
)
study23_summary

#STUDY 24
data24 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Bix+supplemental+data.xlsx",
  sheet = "Vanita's Own Sheet"
)
names(data24)

#skew/kurtosis - CONTROL
k_control24 <- Kurt(
  data24$Saline,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control24
#skew/kurtosis - Treatment
k_treat24 <- Kurt(
  data24$'BIX injection',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat24
# number in each group
n_control24 <- length(na.omit(data24$Saline))
n_treat24 <- length(na.omit(data24$'BIX injection'))
n_control24
n_treat24

#variance ratio
var_control24 <- var(data24$Saline, na.rm=TRUE)
var_treat24 <- var(data24$'BIX injection', na.rm=TRUE)

var_ratio24 <- var_treat24 / var_control24
var_ratio24

study24_summary <- data.frame(
  study = "data24",
  kurt_control = k_control24[1],
  lwr_kurt_control = k_control24[2],
  upr_kurt_control = k_control24[3],
  kurt_treat = k_treat24[1],
  lwr_kurt_treat = k_treat24[2],
  upr_kurt_treat = k_treat24[3],
  var_ratio = var_ratio24,
  n_control = n_control24,
  n_treat = n_treat24
)
study24_summary

#STUDY 25
data25 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Guillot_2025_MCHdata.xlsx"
)
names(data25)

#skew/kurtosis - CONTROL
k_control25 <- Kurt(
  data25$Baseline,
  conf.level = 0.95,
  na.rm = TRUE
)

k_control25
#skew/kurtosis - Treatment
k_treat25 <- Kurt(
  data25$MCH,
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat25
# number in each group
n_control25 <- length(na.omit(data25$Baseline))
n_treat25 <- length(na.omit(data25$MCH))
n_control25
n_treat25

#variance ratio
var_control25 <- var(data25$Baseline, na.rm=TRUE)
var_treat25 <- var(data25$MCH, na.rm=TRUE)

var_ratio25 <- var_treat25 / var_control25
var_ratio25

study25_summary <- data.frame(
  study = "data25",
  kurt_control = k_control25[1],
  lwr_kurt_control = k_control25[2],
  upr_kurt_control = k_control25[3],
  kurt_treat = k_treat25[1],
  lwr_kurt_treat = k_treat25[2],
  upr_kurt_treat = k_treat25[3],
  var_ratio = var_ratio25,
  n_control = n_control25,
  n_treat = n_treat25
)
study25_summary

#STUDY 26
data26 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Different_intensity_exercise_and_verbal_spatial_working_memory_Data.xlsx"
)
names(data26)

#skew/kurtosis - CONTROL
k_control26 <- Kurt(
  data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 4],
  conf.level = 0.95,
  na.rm = TRUE
)
k_control26
#skew/kurtosis - Treatment
k_treat26 <- Kurt(
  data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 2],
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat26

# number in each group
n_control26 <- length(na.omit(data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 4]))
n_treat26   <- length(na.omit(data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 2]))
n_control26
n_treat26

#variance ratio
var_control26 <- var(data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 4], na.rm=TRUE)
var_treat26   <- var(data26$'AC-VM'[data26$'Group(1=LG,2=MG,3=HG,4=CG)' == 2], na.rm=TRUE)

var_ratio26 <- var_treat26 / var_control26
var_ratio26

study26_summary <- data.frame(
  study = "data26",
  kurt_control = k_control26[1],
  lwr_kurt_control = k_control26[2],
  upr_kurt_control = k_control26[3],
  kurt_treat = k_treat26[1],
  lwr_kurt_treat = k_treat26[2],
  upr_kurt_treat = k_treat26[3],
  var_ratio = var_ratio26,
  n_control = n_control26,
  n_treat = n_treat26
)
study26_summary

#STUDY 27
data27 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/EXP4_RM.xlsx"
)
names(data27)

#skew/kurtosis - CONTROL
k_control27 <- Kurt(
  data27$'POST_iTBS+sham-γtACS',
  conf.level = 0.95,
  na.rm = TRUE
)

k_control27
#skew/kurtosis - Treatment
k_treat27 <- Kurt(
  data27$'POST_iTBS+γtACS',
  conf.level = 0.95,
  na.rm = TRUE
)

k_treat27
# number in each group
n_control27 <- length(na.omit(data27$'POST_iTBS+sham-γtACS'))
n_treat27 <- length(na.omit(data27$'POST_iTBS+γtACS'))
n_control27
n_treat27

#variance ratio
var_control27 <- var(data27$'POST_iTBS+sham-γtACS', na.rm=TRUE)
var_treat27 <- var(data27$'POST_iTBS+γtACS', na.rm=TRUE)

var_ratio27 <- var_treat27 / var_control27
var_ratio27

study27_summary <- data.frame(
  study = "data27",
  kurt_control = k_control27[1],
  lwr_kurt_control = k_control27[2],
  upr_kurt_control = k_control27[3],
  kurt_treat = k_treat27[1],
  lwr_kurt_treat = k_treat27[2],
  upr_kurt_treat = k_treat27[3],
  var_ratio = var_ratio27,
  n_control = n_control27,
  n_treat = n_treat27
)
study27_summary

#STUDY 28
data28 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/Wellmind_ms_data.xlsx"
)
names(data28)

#skew/kurtosis - CONTROL
k_control28 <- Kurt(
  data28$'Post- selfcompassion'[data28$'Group' == 'Control'],
  conf.level = 0.95,
  na.rm = TRUE
)
k_control28
#skew/kurtosis - Treatment
k_treat28 <- Kurt(
  data28$'Post- selfcompassion'[data28$'Group' == 'WellMind'],
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat28

# number in each group
n_control28 <- length(na.omit(data28$'Post- selfcompassion'[data28$'Group' == 'Control']))
n_treat28   <- length(na.omit(data28$'Post- selfcompassion'[data28$'Group' == 'WellMind']))
n_control28
n_treat28

#variance ratio
var_control28 <- var(data28$'Post- selfcompassion'[data28$'Group' == 'Control'], na.rm=TRUE)
var_treat28  <- var(data28$'Post- selfcompassion'[data28$'Group' == 'WellMind'], na.rm=TRUE)

var_ratio28 <- var_treat28 / var_control28
var_ratio28

study28_summary <- data.frame(
  study = "data28",
  kurt_control = k_control28[1],
  lwr_kurt_control = k_control28[2],
  upr_kurt_control = k_control28[3],
  kurt_treat = k_treat28[1],
  lwr_kurt_treat = k_treat28[2],
  upr_kurt_treat = k_treat28[3],
  var_ratio = var_ratio28,
  n_control = n_control28,
  n_treat = n_treat28
)
study28_summary

#STUDY 29
data29 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/MSTP_ms_Data_Dryad.xlsx"
)
names(data29)

#skew/kurtosis - CONTROL
k_control29 <- Kurt(
  na.omit(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'B']),
  conf.level = 0.95,
  na.rm = TRUE
)
k_control29
#skew/kurtosis - Treatment
k_treat29 <- Kurt(
  na.omit(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'A']),
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat29

# number in each group
n_control29 <- length(na.omit(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'B']))
n_treat29   <- length(na.omit(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'A']))
n_control29
n_treat29

#variance ratio
var_control29 <- var(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'B'], na.rm=TRUE)
var_treat29  <- var(data29$'Post_SA_Speed'[data29$'GroupIDs' == 'A'], na.rm=TRUE)

var_ratio29 <- var_treat29 / var_control29
var_ratio29

study29_summary <- data.frame(
  study = "data29",
  kurt_control = k_control29[1],
  lwr_kurt_control = k_control29[2],
  upr_kurt_control = k_control29[3],
  kurt_treat = k_treat29[1],
  lwr_kurt_treat = k_treat29[2],
  upr_kurt_treat = k_treat29[3],
  var_ratio = var_ratio29,
  n_control = n_control29,
  n_treat = n_treat29
)
study29_summary

#STUDY 30
data30 <- read_excel(
  "/Users/vanita/Downloads/StatStudy/characteristics_behavioral_data.xlsx"
)
names(data30)

#skew/kurtosis - CONTROL
k_control30 <- Kurt(
  data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'waiting group'],
  conf.level = 0.95,
  na.rm = TRUE
)
k_control30
#skew/kurtosis - Treatment
k_treat30 <- Kurt(
  data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'diet group'],
  conf.level = 0.95,
  na.rm = TRUE
)
k_treat30

# number in each group
n_control30 <- length(na.omit(data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'waiting group']))
n_treat30   <- length(na.omit(data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'diet group']))
n_control30
n_treat30

#variance ratio
var_control30 <- var(data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'waiting group'], na.rm=TRUE)
var_treat30  <- var(data30$'weight_change_t0_minus_t1_kg'[data30$'group' == 'diet group'], na.rm=TRUE)

var_ratio30 <- var_treat30 / var_control30
var_ratio30

study30_summary <- data.frame(
  study = "data30",
  kurt_control = k_control30[1],
  lwr_kurt_control = k_control30[2],
  upr_kurt_control = k_control30[3],
  kurt_treat = k_treat30[1],
  lwr_kurt_treat = k_treat30[2],
  upr_kurt_treat = k_treat30[3],
  var_ratio = var_ratio30,
  n_control = n_control30,
  n_treat = n_treat30
)
study30_summary

#combining studies into a MEGA data frame
all_studies <- rbind(study1_summary, study2_summary,study3_summary,study4_summary,
                     study5_summary,study6_summary, study7_summary,study8_summary,
                     study9_summary, study10_summary, study11_summary, study12_summary,
                     study13_summary, study14_summary, study15_summary, study16_summary,
                     study17_summary, study18_summary, study19_summary, study20_summary,
                     study21_summary,study22_summary, study23_summary, study24_summary,
                     study25_summary, study26_summary, study27_summary, study28_summary,
                     study29_summary, study30_summary

)
all_studies

library(dplyr)
library(tidyr)

forest_df <- all_studies %>%
  pivot_longer(
    cols = c(kurt_control, kurt_treat,
             lwr_kurt_control, lwr_kurt_treat,
             upr_kurt_control, upr_kurt_treat),
    names_to = c(".value", "group"),
    names_pattern = "(kurt|lwr_kurt|upr_kurt)_(control|treat)"
  )
forest_df


## FOREST PLOT BY DATA NUMBER ##
forest_df <- forest_df %>%
  mutate(
    group = factor(group, levels = c("control", "treat"))
  )
forest_df <- forest_df %>%
  mutate(
    study_group = interaction(study, group, sep = " – ", lex.order = TRUE)
  )

ggplot(
  forest_df,
  aes(
    x = kurt,
    y = study_group,
    xmin = lwr_kurt,
    xmax = upr_kurt,
    color = group
  )
) +
  geom_point(size = 2) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Kurtosis (scaled, normal = 0)",
    y = "Study / Group",
    title = "Forest Plot of Kurtosis (with 95% CI) by Treatment Group",
    color = "Group"
  ) +
  theme_minimal()

#summary table
all_studies <- all_studies %>%
  mutate(
    total_n = n_control + n_treat
  )
sample_size_summary <- all_studies %>%
  summarise(
    mean_n   = mean(total_n),
    median_n = median(total_n),
    min_n    = min(total_n),
    max_n    = max(total_n)
  )
sample_size_summary

vr_summary <- all_studies %>%
  summarise(
    mean_vr   = mean(var_ratio),
    median_vr = median(var_ratio),
    min_vr    = min(var_ratio),
    max_vr    = max(var_ratio),
    iqr_vr    = IQR(var_ratio)
  )
vr_summary

#control kurtosis
kurt_control_summary <- all_studies %>%
  summarise(
    mean_kurt   = mean(kurt_control),
    median_kurt = median(kurt_control),
    min_kurt    = min(kurt_control),
    max_kurt    = max(kurt_control)
  )
kurt_control_summary

#treatment kurtosis
kurt_treat_summary <- all_studies %>%
  summarise(
    mean_kurt   = mean(kurt_treat),
    median_kurt = median(kurt_treat),
    min_kurt    = min(kurt_treat),
    max_kurt    = max(kurt_treat)
  )
kurt_treat_summary

summary_table <- tibble(
  Metric = c(
    "Total sample size",
    "Variance ratio (VR)",
    "Kurtosis (Control)",
    "Kurtosis (Treatment)",
    "C ratio (Control)",
    "C ratio (Treatment)"
  ),
  Mean = c(
    mean(all_studies$total_n),
    mean(all_studies$var_ratio),
    mean(all_studies$kurt_control),
    mean(all_studies$kurt_treat),
    mean(C_control, na.rm = TRUE),
    mean(C_treat, na.rm = TRUE)

  ),
  Median = c(
    median(all_studies$total_n),
    median(all_studies$var_ratio),
    median(all_studies$kurt_control),
    median(all_studies$kurt_treat),
    median(C_control, na.rm = TRUE),
    median(C_treat, na.rm = TRUE)
  ),
  Q1 = c(
    quantile(all_studies$total_n, 0.25, na.rm = TRUE),
    quantile(all_studies$var_ratio, 0.25, na.rm = TRUE),
    quantile(all_studies$kurt_control, 0.25, na.rm = TRUE),
    quantile(all_studies$kurt_treat, 0.25, na.rm = TRUE),
    quantile(C_control, 0.25, na.rm = TRUE),
    quantile(C_treat, 0.25, na.rm = TRUE)
  ),
  Q3 = c(
    quantile(all_studies$total_n, 0.75, na.rm = TRUE),
    quantile(all_studies$var_ratio, 0.75, na.rm = TRUE),
    quantile(all_studies$kurt_control, 0.75, na.rm = TRUE),
    quantile(all_studies$kurt_treat, 0.75, na.rm = TRUE),
    quantile(C_control, 0.75, na.rm = TRUE),
    quantile(C_treat, 0.75, na.rm = TRUE)
  )
)

summary_table

library(kableExtra)
library(dplyr)

summary_table %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  kbl(
    caption = "<span style='color:#333333; font-weight:bold;'>
    Summary of Sample Size, Variance Ratio, Kurtosis, and C Ratios Across Studies
    </span>",
    align = c("l", rep("c", 6)),
    booktabs = TRUE,
    linesep = ""
  ) %>%
  kable_styling(
    full_width = TRUE,
    position = "center",
    font_size = 14,
    bootstrap_options = c("striped", "hover", "condensed")
  ) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, bold = TRUE)


## META DATA ##
kurt_se <- function(n_val) {
  se <- sqrt(24/n_val)
  return(se)
}

make_study_meta <- function(study_name,
                            k_control, k_treat,
                            n_control, n_treat,
                            var_control, var_treat) {
  
  data.frame(
    study = study_name,
    group = c("control", "treat"),
    kurt = c(k_control[1], k_treat[1]),
    lwr = c(k_control[2], k_treat[2]),
    upr = c(k_control[3], k_treat[3]),
    SE = c(
      kurt_se(n_control),
      kurt_se(n_treat) 
    ),
    n = c(n_control, n_treat),
    var = c(var_control, var_treat)
  )
}

study_list <- list(
  make_study_meta("data1",  k_control1,  k_treat1,  n_control1,  n_treat1,  var_control1,  var_treat1),
  make_study_meta("data2",  k_control2,  k_treat2,  n_control2,  n_treat2,  var_control2,  var_treat2),
  make_study_meta("data3",  k_control3,  k_treat3,  n_control3,  n_treat3,  var_control3,  var_treat3),
  make_study_meta("data4",  k_control4,  k_treat4,  n_control4,  n_treat4,  var_control4,  var_treat4),
  make_study_meta("data5",  k_control5,  k_treat5,  n_control5,  n_treat5,  var_control5,  var_treat5),
  make_study_meta("data6",  k_control6,  k_treat6,  n_control6,  n_treat6,  var_control6,  var_treat6),
  make_study_meta("data7",  k_control7,  k_treat7,  n_control7,  n_treat7,  var_control7,  var_treat7),
  make_study_meta("data8",  k_control8,  k_treat8,  n_control8,  n_treat8,  var_control8,  var_treat8),
  make_study_meta("data9",  k_control9,  k_treat9,  n_control9,  n_treat9,  var_control9,  var_treat9),
  make_study_meta("data10", k_control10, k_treat10, n_control10, n_treat10, var_control10, var_treat10),
  
  make_study_meta("data11", k_control11, k_treat11, n_control11, n_treat11, var_control11, var_treat11),
  make_study_meta("data12", k_control12, k_treat12, n_control12, n_treat12, var_control12, var_treat12),
  make_study_meta("data13", k_control13, k_treat13, n_control13, n_treat13, var_control13, var_treat13),
  make_study_meta("data14", k_control14, k_treat14, n_control14, n_treat14, var_control14, var_treat14),
  make_study_meta("data15", k_control15, k_treat15, n_control15, n_treat15, var_control15, var_treat15),
  make_study_meta("data16", k_control16, k_treat16, n_control16, n_treat16, var_control16, var_treat16),
  make_study_meta("data17", k_control17, k_treat17, n_control17, n_treat17, var_control17, var_treat17),
  make_study_meta("data18", k_control18, k_treat18, n_control18, n_treat18, var_control18, var_treat18),
  make_study_meta("data19", k_control19, k_treat19, n_control19, n_treat19, var_control19, var_treat19),
  make_study_meta("data20", k_control20, k_treat20, n_control20, n_treat20, var_control20, var_treat20),
  
  make_study_meta("data21", k_control21, k_treat21, n_control21, n_treat21, var_control21, var_treat21),
  make_study_meta("data22", k_control22, k_treat22, n_control22, n_treat22, var_control22, var_treat22),
  make_study_meta("data23", k_control23, k_treat23, n_control23, n_treat23, var_control23, var_treat23),
  make_study_meta("data24", k_control24, k_treat24, n_control24, n_treat24, var_control24, var_treat24),
  make_study_meta("data25", k_control25, k_treat25, n_control25, n_treat25, var_control25, var_treat25),
  make_study_meta("data26", k_control26, k_treat26, n_control26, n_treat26, var_control26, var_treat26),
  make_study_meta("data27", k_control27, k_treat27, n_control27, n_treat27, var_control27, var_treat27),
  make_study_meta("data28", k_control28, k_treat28, n_control28, n_treat28, var_control28, var_treat28),
  make_study_meta("data29", k_control29, k_treat29, n_control29, n_treat29, var_control29, var_treat29),
  make_study_meta("data30", k_control30, k_treat30, n_control30, n_treat30, var_control30, var_treat30)
)

all_meta <- bind_rows(study_list)
all_meta 

library(metafor)

res <- rma(yi = kurt, sei = SE, mods = ~ group, data = all_meta)

summary(res)
#grouptreat estimate being 0.4569 shows that treatment tends to increase kurtosis but the pvalue of 0.3831 
#shows that it is not significant. I^2 = 76.11% means 76% of variability in effect sizes is 
#NOT due to sampling error, confirmed by how the test for residual heterogeneity is significant. 
#Our moderator (treatment vs control) accounts for 0.52% of heterogeneity -> Whether something is “treatment” or “control” barely explains differences in kurtosis

#Overall there is no evidence treatment increase kurtosis, VERY high variability across studies, moderator is a weak explanation of kurtosis differences

## BOOTSTRAPPING FOR VR ##
list_vr <- list(
  var_ratio1,
  var_ratio2,
  var_ratio3,
  var_ratio4,
  var_ratio5,
  var_ratio6,
  var_ratio7,
  var_ratio8,
  var_ratio9,
  var_ratio10,
  var_ratio11,
  var_ratio12,
  var_ratio13,
  var_ratio14,
  var_ratio15,
  var_ratio16,
  var_ratio17,
  var_ratio18,
  var_ratio19,
  var_ratio20,
  var_ratio21,
  var_ratio22,
  var_ratio23,
  var_ratio24,
  var_ratio25,
  var_ratio26,
  var_ratio27,
  var_ratio28,
  var_ratio29,
  var_ratio30
)
list_vr

df_vr <- data.frame(
  study = 1:30,
  vr = unlist(list_vr)
)
df_vr

library(boot)

boot_fun <- function(data, indices) {
  d <- data[indices, ]
  mean(d$vr)
}

#bootstrap mean sample VR 
set.seed(123)

b <- boot(
  data = df_vr,
  statistic = boot_fun,
  R = 5000
)
b
#get bootstrap percentile interval (perc) and adjusted bootstrap percentile (BCa) interval
boot.ci(b, conf = 0.95, type = c("perc", "bca"))

#the estimate shows 1.88 suggesting treatment variance is higher, however VR is inherently right skewed 
#which could be misleading. Our CI's do not contain 1, which is consistent

#bootstrap log mean sample VR
library(boot)

boot_fun_log <- function(data, indices) {
  d <- data[indices, ]
  mean(log(d$vr))
}

set.seed(123)

b_log <- boot(
  data = df_vr,
  statistic = boot_fun_log,
  R = 5000
)
b_log

#get bootstrap percentile interval (perc) and adjusted bootstrap percentile (BCa) interval
boot.ci(b_log, conf = 0.95, type = c("perc", "bca"))


#when we logged the mean sample VR, we get 0.0766 as our estimate. Back transforming via 
#exp(0.0766) we get 1.07961 which is close to 1, suggesting not much of a difference in treatment 
#variance compared to control variance. Our CI's contain 0 now, showing that logging the VR values 
#affects the significance of the VR we bootstrapped 
