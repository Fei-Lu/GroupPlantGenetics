
library(tidyverse)
args <- commandArgs(trailingOnly = T)


# fdResDir <- "003_twoWay_proportion0.1_genetration100/004_fd/003_fdRes/"
# individualLocalAncestryDir <- "003_twoWay_proportion0.1_genetration100/004_fd/004_individualLocalAncestry/"
fdResDir <- args[1]
individualLocalAncestryDir <- args[2]

df <- dir(fdResDir,full.names = T) %>% 
    map(~read_csv(.x) %>% 
            mutate(AdmixtureProportion=str_remove(str_split(basename(.x),"_")[[1]][[2]],"proportion"),AdmixtureGenetration=str_remove(str_split(basename(.x),"_")[[1]][[3]],"genetration"),DivergenceT=str_remove(str_split(basename(.x),"_")[[1]][[4]],"divergence"),
                   DivergenceTime=str_remove(DivergenceT,"Ne"),
                   P2=str_c("tsk_",str_remove(str_split(basename(.x),"_")[[1]][[6]],".csv")))) %>% 
    bind_rows()

df_fd_global <- df %>% 
    group_by(P2,AdmixtureProportion,AdmixtureGenetration,DivergenceTime) %>% 
    mutate(fdModify=ifelse(is.na(D) | D <= 0, 0, ifelse(fd > 1, 0, ifelse(fd < 0, 0, fd)))) %>%
    summarise(fd_mean=mean(fdModify))

df_fd <- df %>% 
    left_join(y = df_fd_global, by = c("P2","AdmixtureProportion","AdmixtureGenetration","DivergenceTime")) %>% 
    mutate(fdModify=ifelse(is.na(D) | D <= 0, 0, ifelse(fd > 1, 0, ifelse(fd < 0, 0, fd)))) %>%
    mutate(IfIntrogressionByfd=ifelse(fdModify>fd_mean,1,0)) %>% 
    select(-fdModify) %>% 
    select(scaffold:sitesUsed, P2,AdmixtureProportion,AdmixtureGenetration,DivergenceTime,IfIntrogressionByfd)

df_fd %>% 
    mutate(Donor=ifelse(IfIntrogressionByfd==1, 0, 4)) %>% 
    mutate(P2=str_c("Proportion",AdmixtureProportion,"_genetration",AdmixtureGenetration,"_divergence",DivergenceTime,"Ne_",P2)) %>% 
    select(-c(AdmixtureProportion,AdmixtureGenetration,DivergenceTime)) %>% 
    group_by(P2) %>% 
    group_walk(~write_tsv(.x, file.path(individualLocalAncestryDir,str_c("SimulatedLocalAncestry_",.y,"simulated.txt.gz"))))


