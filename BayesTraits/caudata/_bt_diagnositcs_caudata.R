##############################################################################################
##############################################################################################
##############################################################################################
### BayesTraits - diagnostics and exploration

# This script was used to evaluate the MCMC outputs from BayesTraits, for the Caudata
# dataset of Liedtke et al. 2022.
##############################################################################################
##############################################################################################
##############################################################################################

##############################################################################################
# set working directory and load required libraries

setwd("~/Documents/amphibian_diversity_project/2021/_SI_scripts/BayesTraits/caudata/")

library(tidyverse)
library(coda)
library(qgraph)

##############################################################################################
# load BayesTraits output

# list files
log_files<-list.files(pattern = "Log.txt")
log_files

# extract run names
model_names<-gsub(log_files, pattern=".*\\_r(\\d+).*.*", replacement = paste0("run ","\\1"))
model_names


# load all results
bt_out<-list()
for(i in 1:length(log_files)){
  bt_out[[i]]<-read_tsv(log_files[i], skip=grep(pattern = "Lh", read_lines(log_files[i], n_max=100))-1, col_names = TRUE) #%>%
}


### add run names to list
names(bt_out)<-model_names

# stack as single dataframe
bt_df<-bind_rows(bt_out, .id="model")


##############################################################################################
# Plot chains

# plot likelihood traces [overlapping]

ggplot(bt_df, aes(x=Iteration, y=Lh, color=model)) +
  geom_line(alpha=0.5) +
  #geom_smooth(se = FALSE) +
  theme_bw()


###trim excess burnin if necessary

if(FALSE){
  
  for(i in 1:length(bt_out)) {
    bt_out[[i]]<-bt_out[[i]] %>%
      slice_tail(prop=0.9)
  }
  
  bt_df<-bind_rows(bt_out, .id="model")
  
  ggplot(bt_df, aes(x=Iteration, y=Lh)) +
    geom_line() +
    geom_smooth(se = FALSE) +
    facet_wrap(~model, ncol=2,scales = "free") + 
    theme_bw()
}


# plot likelihood densities 
ggplot(bt_df, aes(x=Lh, fill=model)) +
  geom_density(alpha=0.5, color=NA) +
  theme_bw() 

##############################################################################################
#  get MCMC stats

# converts log file into mcmc object to use for coda package (note: cannot take non-numerics)

MCMCdata<-list()
for(i in 1:length(bt_out)) {
  MCMCdata[[i]]<-coda::mcmc(bt_df[bt_df$model==names(bt_out)[i],grep(colnames(bt_df), pattern="^L|^q|^R")])
}  
names(MCMCdata)<-names(bt_out)

# calculates effective sample sizes
lapply(MCMCdata, FUN =  effectiveSize) %>%
  bind_rows(.id = "model") %>%
  t()

## check whether sampling rate is high enough to avoid autocorrelation
acfplot(MCMCdata[[1]])


################################################################################################
# Summarise best/first chain

#### summarize rates
chain1<-bt_out[[1]]

### summarize all important variables
rates_sum<-chain1 %>%
  select(starts_with("q")) %>%
  pivot_longer(everything(), names_to = "transition", values_to = "rate") %>%
  group_by(transition) %>%
  summarise(mean=mean(rate),
            median=median(rate),
            perc_non_zero=100-(sum(rate==0)/length(rate))*100,
            ESS=effectiveSize(coda::mcmc(rate)),
            HPD=HPDinterval(coda::mcmc(rate)) %>%
              as_tibble() %>%
              mutate_all(round, 3) %>%
              unite("x",c(lower, upper),sep=",") %>% pull(x)
  ) %>%
  arrange(desc(perc_non_zero))

rates_sum
write.csv(rates_sum, "chain1_summary.csv")

#### plot model string rankings
table(chain1$`Model string`) %>% sort()

chain1 %>%
  count(`Model string`) %>%
  arrange(desc(n)) %>%
  slice_max(n, n=25,with_ties = FALSE) %>%
  mutate("Model Rank"=factor(1:n())) %>%
  ggplot(aes(x=`Model Rank`, y=n)) +
  geom_bar(stat="identity") +
  ylab("number of iterations") +
  theme(legend.position = "none") +
  theme_bw()
  
  
#### Plot distributions of rates [with mean and medians]

chain1 %>%
  select(starts_with("q")) %>%
  pivot_longer(everything(), names_to = "transition", values_to = "rate") %>%
  ggplot(aes(x=rate, fill=transition, color=transition)) +
  #geom_histogram(alpha=0.5, color=NA, bins = 15) +
  geom_density(alpha=0.75, fill="grey50",color=NA) +
  geom_vline(data=rates_sum, aes(xintercept=mean)) +
  geom_vline(data=rates_sum, aes(xintercept=median), linetype="dashed") +
  theme_bw() +
  facet_wrap(~transition, scales = "free") +
  theme(legend.position = "none")

### Plot transition network for 0% not zero 

non_zero_rates<-rates_sum %>%
  select(perc_non_zero, transition) %>%
  mutate(transition=str_remove(transition, "q*")) %>%
  separate(transition, into=c("from","to"), sep = 1) %>%
  pivot_wider(names_from = to, values_from = perc_non_zero) %>%
  data.frame(row.names = "from") %>%
  select(rownames(.))

rep.cols<-c("#24b9e9","#55FF7E","#BF9C00","#009E73","#f6776f","gold")
names(rep.cols)<-c("A","S","T","D","V","P")

qgraph(non_zero_rates,
       #title=,
       #nodes:
       labels=rownames(non_zero_rates),
       shape="rectangle",
       vsize=10,
       node.height=0.5,
       label.cex=1.25,
       label.color=rep.cols[rownames(non_zero_rates)],
       # edges:
       edge.labels=T,
       edge.label.cex=2,
       minimum=50, # threshold for edges to include
       edge.width=1.5,
       posCol="black",
       fade=T,
       colFactor=0.9,
       # themes and layout:
       #theme="TeamFortress",
       layout="circular",
       directed=TRUE)




