```R
options(here.quiet = TRUE)
library(tidyverse)
library(survival)
library(survminer)
library(here)
```

    ── [1mAttaching core tidyverse packages[22m ──── tidyverse 2.0.0 ──
    [32m✔[39m [34mdplyr    [39m 1.1.3     [32m✔[39m [34mreadr    [39m 2.1.4
    [32m✔[39m [34mforcats  [39m 1.0.0     [32m✔[39m [34mstringr  [39m 1.5.0
    [32m✔[39m [34mggplot2  [39m 3.4.3     [32m✔[39m [34mtibble   [39m 3.2.1
    [32m✔[39m [34mlubridate[39m 1.9.2     [32m✔[39m [34mtidyr    [39m 1.3.0
    [32m✔[39m [34mpurrr    [39m 1.0.2     
    ── [1mConflicts[22m ────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m    masks [34mstats[39m::lag()
    [36mℹ[39m Use the conflicted package ([3m[34m<http://conflicted.r-lib.org/>[39m[23m) to force all conflicts to become errors
    Loading required package: ggpubr
    
    
    Attaching package: ‘survminer’
    
    
    The following object is masked from ‘package:survival’:
    
        myeloma
    
    
    here() starts at /diskmnt/Projects/Users/ysong/project/jupyter_lab/MMRF_analysis/IA2024Paper1/Fig1/scripts
    



```R
df <-read.delim("../data/baseline_263_visit_clinical.txt")
# Format the categorical variable to appropriate factor levels.
df <- df %>%
  mutate(
    Sex = d_pt_sex,
    Age = d_dx_amm_age,
    BMI = d_dx_amm_bmi,
    Transplant = d_amm_tx_asct_1st,
    PFS = censpfs,
    PFS.time = ttcpfs,
    OS = censos,
    OS.time = ttcos,
    Treatment = case_when(
      d_tx_induction_cat == "chemo_imid_pi_steroid" ~ "imid_pi_steroid",
      d_tx_induction_cat == "chemo_pi_steroid" ~ "pi_steroid",
      TRUE ~ d_tx_induction_cat
    ),
    risk = case_when(
      davies_based_risk %in% c("no_risk_data", "not_calculable") ~ NA_character_,
      TRUE ~ as.character(davies_based_risk)
    ),
    ECOG = case_when(
      d_dx_amm_ecog %in% c("3", "4") ~ ">=3",
      TRUE ~ as.character(d_dx_amm_ecog)
    ),
    Race = case_when(
      d_pt_race_1 %in% c("asian_nos", "unknown") ~ "others/unknown",
      d_pt_race_1 %in% c("black_african_american") ~ "black",
      TRUE ~ as.character(d_pt_race_1)
    ),
    ISS_stage = case_when(
      d_dx_amm_iss_stage %in% c("1") ~ "stage I",
      d_dx_amm_iss_stage %in% c("2") ~ "stage II",
      d_dx_amm_iss_stage %in% c("3") ~ "stage III",
      TRUE ~ as.character(d_dx_amm_iss_stage)
    )
  ) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("pi_steroid", "imid_pi_steroid", "imid_steroid")),
    risk = factor(risk, levels = c("standard_risk", "high_risk")),
    ECOG = factor(ECOG, levels = c("0", "1", "2", ">=3")),
    Race = factor(Race, levels = c("white", "black", "others/unknown")),
    ISS_stage = factor(ISS_stage, levels = c("stage I", "stage II", "stage III"))
  )

# Display the first few rows of the transformed data frame
df %>%
  select(Sex, BMI, Transplant, Treatment, risk, ECOG, Race, ISS_stage)%>% head()
```


<table class="dataframe">
<caption>A data.frame: 6 × 8</caption>
<thead>
	<tr><th></th><th scope=col>Sex</th><th scope=col>BMI</th><th scope=col>Transplant</th><th scope=col>Treatment</th><th scope=col>risk</th><th scope=col>ECOG</th><th scope=col>Race</th><th scope=col>ISS_stage</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>male  </td><td>31.20</td><td>no </td><td>pi_steroid     </td><td>high_risk    </td><td>0</td><td>white</td><td>stage III</td></tr>
	<tr><th scope=row>2</th><td>male  </td><td>24.26</td><td>no </td><td>imid_pi_steroid</td><td>NA           </td><td>1</td><td>white</td><td>stage I  </td></tr>
	<tr><th scope=row>3</th><td>male  </td><td>25.21</td><td>yes</td><td>imid_pi_steroid</td><td>high_risk    </td><td>1</td><td>white</td><td>stage II </td></tr>
	<tr><th scope=row>4</th><td>female</td><td>33.12</td><td>no </td><td>pi_steroid     </td><td>NA           </td><td>2</td><td>white</td><td>NA       </td></tr>
	<tr><th scope=row>5</th><td>male  </td><td>42.78</td><td>no </td><td>imid_pi_steroid</td><td>high_risk    </td><td>1</td><td>white</td><td>stage II </td></tr>
	<tr><th scope=row>6</th><td>male  </td><td>25.73</td><td>yes</td><td>imid_pi_steroid</td><td>standard_risk</td><td>0</td><td>white</td><td>stage II </td></tr>
</tbody>
</table>




```R
# Figure 1B Hazard Ratio of clinical features associated with PFS
options(repr.plot.width=11, repr.plot.height=12)

formula <- as.formula(paste("Surv(PFS.time, PFS) ~", "Age +Sex+ Race+ BMI +ISS_stage + Transplant +Treatment+risk+ECOG"))
model <- coxph(formula, data = as.data.frame(df))
p <- ggforest(model, data = as.data.frame(df),
              main = 'Hazard ratio of progression free survival', 
              cpositions = c(0.0, 0.25, 0.45),
              fontsize = 1.3,
              refLabel = '1.00',
              noDigits = 2)

p
```


    
![png](output_2_0.png)
    

