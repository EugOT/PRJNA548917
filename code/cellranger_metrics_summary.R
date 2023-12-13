library(tidyverse)
metrics_summary_SRR9303096 <- read_csv("/data/PRJNA548917/cellranger/SRR9303096/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303096")
metrics_summary_SRR9303097 <- read_csv("/data/PRJNA548917/cellranger/SRR9303097/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303097")
metrics_summary_SRR9303098 <- read_csv("/data/PRJNA548917/cellranger/SRR9303098/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303098")
metrics_summary_SRR9303099 <- read_csv("/data/PRJNA548917/cellranger/SRR9303099/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303099")
metrics_summary_SRR9303100 <- read_csv("/data/PRJNA548917/cellranger/SRR9303100/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303100")
metrics_summary_SRR9303101 <- read_csv("/data/PRJNA548917/cellranger/SRR9303101/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303101")
metrics_summary_SRR9303102 <- read_csv("/data/PRJNA548917/cellranger/SRR9303102/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303102")
metrics_summary_SRR9303103 <- read_csv("/data/PRJNA548917/cellranger/SRR9303103/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303103")
metrics_summary_SRR9303104 <- read_csv("/data/PRJNA548917/cellranger/SRR9303104/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303104")
metrics_summary_SRR9303105 <- read_csv("/data/PRJNA548917/cellranger/SRR9303105/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303105")
metrics_summary_SRR9303106 <- read_csv("/data/PRJNA548917/cellranger/SRR9303106/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303106")
metrics_summary_SRR9303107 <- read_csv("/data/PRJNA548917/cellranger/SRR9303107/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303107")
metrics_summary_SRR9303108 <- read_csv("/data/PRJNA548917/cellranger/SRR9303108/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303108")
metrics_summary_SRR9303109 <- read_csv("/data/PRJNA548917/cellranger/SRR9303109/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303109")
metrics_summary_SRR9303110 <- read_csv("/data/PRJNA548917/cellranger/SRR9303110/outs/metrics_summary.csv") %>% mutate(Run = "SRR9303110")
metrics_summary <-
    bind_rows(
        metrics_summary_SRR9303096,
        metrics_summary_SRR9303097,
        metrics_summary_SRR9303098,
        metrics_summary_SRR9303099,
        metrics_summary_SRR9303100,
        metrics_summary_SRR9303101,
        metrics_summary_SRR9303102,
        metrics_summary_SRR9303103,
        metrics_summary_SRR9303104,
        metrics_summary_SRR9303105,
        metrics_summary_SRR9303106,
        metrics_summary_SRR9303107,
        metrics_summary_SRR9303108,
        metrics_summary_SRR9303109,
        metrics_summary_SRR9303110)

metrics_summary |>
    select("Estimated Number of Cells", "Run")

write_tsv(metrics_summary, "/data/PRJNA548917/metrics_summary.tsv")

