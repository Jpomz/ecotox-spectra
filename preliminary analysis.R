# sort data 
# Justin Pomeranz
# jfpomeranz@gmail.com

library(plyr)
library(tidyverse)
breaks<-2^seq(-25,15)

# function to bin individual weights, and center bins at x = 0
bin.center <- function(data, var, breaks, ...){
  binned_hist = hist(data[[var]], breaks = breaks,
                     include.lowest = TRUE,
                     plot = FALSE)
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  break_width = breaks_offset - breaks_orig
  count = binned_hist$counts 
  dataout = data.frame(
    # count = binned_hist$counts,
    # mids = binned_hist$mids,
    # break_width = break_width,
    # count_corrected = count_corrected,
    log_count_corrected = log10(count / break_width),
    log_mids = log10(binned_hist$mids)
  )
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  dataout$log_mids_center = scale(dataout$log_mids,
                                  center = TRUE,
                                  scale = FALSE)
  dataout
}


# data from Scott Morton Travis Schmidt email 19 Jan 2018
# Modified to csv, changed colnames "estimated DM weight per individual (mg)" --> "ind_mg"; "estimated DM weight for all individuals of specified size class (mg)" --> "tot_mg"; deleted "a", "b" and "base" columns from original

# row 1329 in original did not have a value in "a" column. added value of "0.0070" found in "Scott_working_taxa_list.xlsx"
# added "a" and "b" values to row 1685

# read in modified data
dat <- read_csv("axl007.csv")[,c(2:3, 9:10)]
# filter out rows where ind_mg / tot_mg ==0
dat <- dat %>% filter(ind_mg != 0)

# what are the A1 A2... ?
# are the treatments just day 30?
# there are "controls" in the lab_smp_id, and in the treatment columns. When lab_smp_id == Alphanumeric && treatment == control, what is that? another control sampled at day 30???

# control rows 1:615
ctrl <- dat[1:615,]
ctrl <- ctrl[rep(seq_len(nrow(ctrl)),
                 times = ctrl$num_indiv), ]
ctrl <- ctrl %>%
  separate(lab_smp_id, c("rep", "day"), sep = " Day ")

ctrl$rep <- gsub("Control ", "", fixed = TRUE, x = ctrl$rep)


ctrl.bin.center <- ddply(ctrl,
                         .(rep, day),
                         bin.center, var = "ind_mg", breaks = breaks)

ggplot(ctrl.bin.center, aes(x = log_mids_center,
                     y = log_count_corrected,
                     color = day)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm", alpha = 0.1) +
  theme_classic()

ctrl.lm <- lm(log_count_corrected ~ log_mids_center * rep*day, data = ctrl.bin.center)
summary(ctrl.lm)


# separate treatment into chem and conc(entrtion) 
dat <- dat %>%
   separate(treatment, c("chem", "conc"), fill = "right", remove = FALSE, convert = TRUE) %>% replace_na(replace = list(conc = 0))


# rep rows by num_ind
dat <- dat[rep(seq_len(nrow(dat)), times = dat$num_indiv), ]





# ddply solution ####
bin.df <- ddply(dat, .(chem, conc), bin, var = "ind_mg", breaks = breaks)

ggplot(bin.df, aes(x = log_break,
                   y = log_count_corrected,
                   color = chem)) +
  geom_point()

# color = chem, size = conc
ggplot(bin.df, aes(x = log_break,
                   y = log_count_corrected,
                   color = chem,
                   size = conc/100)) +
  geom_point() +
  theme_classic()

#













# dat %>% group_by(treatment) %>%
#   do(data.frame(bin(., var = "ind_mg", breaks = breaks)))

# splitting into list ####
# bin.list <- dat %>% 
#   split(.$treatment) %>% 
#   map(bin, var = "ind_mg", breaks = breaks)
# 
# bin.df <- ldply(bin.list)
# 
# ggplot(bin.df, aes(x = log_break,
#                    y = log_count_corrected,
#                    color = .id)) +
#   geom_point()
# 
# # color = chem, size = conc
# ggplot(bin.df, aes(x = log_break,
#                    y = log_count_corrected,
#                    color = chem,
#                    size = conc)) +
#   geom_point()

