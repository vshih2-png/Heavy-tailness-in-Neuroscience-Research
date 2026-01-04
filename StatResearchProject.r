
library(readxl)
library(DescTools)
library(moments)
library(ggplot2)

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
  data3$Ctrl,
  conf.level = 0.95,
  na.rm = TRUE
)
k_control3
#skew/kurtosis - Treatment
k_treat3 <- Kurt(
  data3$`D/A +INS`,
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


#combining studies into a MEGA data frame 
all_studies <- rbind(study1_summary, study2_summary,study3_summary
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

ggplot(forest_df,
       aes(x = kurt,
           y = interaction(study, group),
           xmin = lwr_kurt,
           xmax = upr_kurt)) +
  geom_point(size = 2) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Kurtosis (scaled, normal = 0)",
    y = "Study / Group",
    title = "Forest Plot of Kurtosis Across Studies",
    subtitle = "Assessing heavy-tailedness"
  ) +
  theme_minimal()


