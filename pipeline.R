library(tidyverse)
library(data.table)
library(utile.tools)

##

a = fread('review/HMP2/Metatranscriptomes/ecs.tsv')
a[1:15,1:5]

a2 = a %>%
  mutate(EC = str_remove(`# Gene Family`,":.*"),
         Name = str_remove(`# Gene Family`,".*: ")) %>%
  select(EC, Name, everything(),-`# Gene Family`) %>%
  group_by(EC) %>%
  slice(1) %>%
  ungroup()
a2
a_stat = a2 %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(EC,Name) %>%
  summarise(n = n(),
            n_more0 = length(count[count > 0]),
            min = min(count),
            min_more0 = min(count[count > 0]),
            max = max(count),
            max_more0 = max(count[count > 0]),
            mean = mean(count),
            mean_more0 = mean(count[count > 0]),
            median = median(count),
            median_more0 = median(count[count > 0]),
            sd = sd(count),
            sd_more0 = sd(count[count > 0]),) %>%
  rename(Name_EC = Name) %>%
  mutate(n_per_enz = round((n_more0/n)*100,2)) %>%
  select(EC, Name_EC, n_per_enz)

head(a_stat)

# all toxins all enzymes 67 - tox_enz2 ####

tox_name = as_tibble(fread('tables/tox_list.tab'))

tox_name = tox_name %>%
  select(Name,KEGG_ID) %>%
  distinct()
tox_name

tox_enz = fread('apr_2021/tox_in_react2.tab')
head(tox_enz)

tox_enz2 = as_tibble(tox_enz) %>%
  select(TOXIN, ENTRY, ENZYME) %>%
  filter(ENZYME != "") %>%
  distinct() %>%
  select(TOXIN, ENZYME) %>%
  distinct() %>%
  separate_rows(ENZYME, convert = TRUE,sep='; ') %>%
  distinct()
tox_enz2

tox_enz2 = tox_enz2  %>%
  left_join(tox_name,by=c('TOXIN' = 'KEGG_ID')) %>%
  rename(TOXIN_Name = Name) %>%
  distinct()
tox_enz2

length(unique(tox_enz2$TOXIN)) #67
length(unique(tox_enz2$TOXIN_Name)) #67
length(unique(tox_enz2$ENZYME)) #609

enz_n_tox = tox_enz2 %>%
  group_by(ENZYME) %>%
  summarise(n_tox_for_enz = length(unique(TOXIN)))
enz_n_tox

length(unique(enz_n_tox$ENZYME)) # 609


#### ecs (tox) - 2082(enz) x 735(pat) - ecs2 #####

ecs = fread('review/HMP2/Metatranscriptomes/ecs.tsv')

ecs2 = ecs %>%
  mutate(EC = str_remove(`# Gene Family`,":.*"),
         Name_enz = str_remove(`# Gene Family`,".*: ")) %>%
  select(EC, Name_enz, everything(),-`# Gene Family`) %>%
  group_by(EC) %>%
  slice(1) %>%
  ungroup()
head(ecs2)
length(unique(ecs2$EC)) #2082

### proteins - 910(prot) x 450 (enz) - prot2#####

prot = fread('review/HMP2/Proteomics/HMP2_proteomics_ecs.tsv')
head(prot)

prot2 = as_tibble(prot) %>%
  filter(Gene != 'UNGROUPED') %>%
  mutate(EC = str_remove(Gene,":.*"),
         Name = str_remove(Gene,".*: ")) %>%
  select(EC, Name, everything(),-Gene) %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(EC, Name) %>%
  summarise(n = n(),
            n_more0 = length(count[count > 0]),
            n_per = round((n_more0/n)*100,2),
            prot = paste(n_more0,n_per,sep='_')) %>%
  ungroup() %>%
  select(EC,prot)
prot2

### merging - tox_group - from toxines #####

tox_enz_ecs = ecs2 %>% 
  inner_join(tox_enz2,by=c('EC'='ENZYME')) %>%
  select(EC, Name_enz, TOXIN, TOXIN_Name,everything())
tox_enz_ecs

length(unique(tox_enz_ecs$EC)) #189
length(unique(tox_enz_ecs$TOXIN)) #47

tox_stat = tox_enz_ecs %>% 
  group_by(TOXIN, TOXIN_Name) %>%
  summarise(n_enz = n(),
            across(where(is.numeric),sum)) %>%
  ungroup() %>%
  mutate(n_enz = as.character(n_enz)) %>%
  #select(-n_enz) %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(TOXIN, TOXIN_Name, n_enz) %>%
  summarise(n = n(),
            n_more0 = length(count[count > 0]),
            n_per = round((n_more0/n)*100,2),
            min = min(count),
            min_more0 = min(count[count > 0]),
            max = max(count),
            max_more0 = max(count[count > 0]),
            mean = mean(count),
            mean_more0 = mean(count[count > 0]),
            median = median(count),
            median_more0 = median(count[count > 0]),
            sd = sd(count),
            sd_more0 = sd(count[count > 0]),)
tox_stat
  
quart = tox_enz_ecs %>% 
  select(-Name_enz, -TOXIN, -TOXIN_Name) %>%
  distinct() %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  distinct(count) %>%
  summarise(quantile(count)) %>%
  pull()
quart


enz_quart_prot = tox_enz_ecs %>% 
  group_by(EC,Name_enz) %>%
  distinct() %>%
  ungroup() %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(EC, Name_enz) %>%
  summarize(median_more0_tox = median(count[count > 0])) %>%
  ungroup() %>%
  mutate(group_enz_tox = case_when(
    between(median_more0_tox, quart[1], quart[2]) ~ "1",
    between(median_more0_tox, quart[2], quart[3]) ~ "2",
    between(median_more0_tox, quart[3], quart[4]) ~ "3",
    between(median_more0_tox, quart[4], quart[5]) ~ "4")) %>%
  left_join(prot2,by = c('EC'))
enz_quart_prot

tox_group = tox_enz_ecs %>% 
  select (EC, Name_enz, TOXIN, TOXIN_Name) %>%
  inner_join(enz_quart_prot,by=c('EC','Name_enz')) %>%
  distinct() %>%
  left_join(a_stat, by = 'EC') %>%
  mutate(prot2 = str_remove(prot,'.*_'),
         prot2 = paste(prot2,'%',sep=' '),
         prot2 = str_replace(prot2,'NA %','NA'),
         n_per_enz = paste(n_per_enz,'%',sep='')) %>%
  group_by(TOXIN, TOXIN_Name) %>%
  mutate(enz_group = paste0('EC: ',EC,'; Enz_Name: ',Name_enz,'; median_value: ',round(median_more0_tox,2),'; quartile: ',group_enz_tox, '; prot: ',prot,collapse = ";"),
         enz_group2 = paste0(Name_enz,' (mRNA: ',n_per_enz,', Q',group_enz_tox,', prot: ',prot2,')',collapse = '; ')) %>%
  ungroup() %>%
  select(TOXIN, TOXIN_Name, enz_group2) %>%
  distinct() %>%
  inner_join(tox_stat,by=c('TOXIN','TOXIN_Name')) %>%
  mutate(n_per2 = paste(n_per,'%',sep=' ')) %>%
  select(TOXIN_Name, n_enz, n_per2, enz_group2, everything())
tox_group


write.table(tox_group,file='review/tables/from_tox.tab',quote=T,row.names=F,col.names=T,sep='\t')
write.table(tox_group,file='/mnt/public_html/suvorova/toxins/review/from_tox.tab',quote=T,row.names=F,col.names=T,sep='\t')

write.table(tox_group,file='review/tables/from_tox_topap.tab',quote=F,row.names=F,col.names=T,sep='\t')
write.table(tox_group,file='/mnt/public_html/suvorova/toxins/review/from_tox_topap.tab',quote=F,row.names=F,col.names=T,sep='\t')


### 16S - 284 bac - bac_16S #####
a=fread('review/HMP2/16S/taxonomic_profiles.tsv',header=T)
head(a)

at=as_tibble(a)
head(at)

at2 = at %>%
  mutate(w1 = word(taxonomy,1),
         w2 = word(taxonomy,2),
         w3 = word(taxonomy,3),
         w4 = word(taxonomy,4),
         w5 = word(taxonomy,5),
         w6 = word(taxonomy,6)) 

at2 = at2 %>%
  mutate(w1 = str_replace(w1,';',''),
         w2 = str_replace(w2,';',''),
         w3 = str_replace(w3,';',''),
         w4 = str_replace(w4,';',''),
         w5 = str_replace(w5,';',''),
         w6 = str_replace(w6,';','')) 


at2 = at2 %>%
  mutate(w1 = str_replace(w1,'^__',''),
         w2 = str_replace(w2,'^__',''),
         w3 = str_replace(w3,'^__',''),
         w4 = str_replace(w4,'^__',''),
         w5 = str_replace(w5,'^__',''),
         w6 = str_replace(w6,'^__',''))

at2 = at2 %>%
  mutate(w1 = str_replace(w1,'^_',''),
         w2 = str_replace(w2,'^_',''),
         w3 = str_replace(w3,'^_',''),
         w4 = str_replace(w4,'^_',''),
         w5 = str_replace(w5,'^_',''),
         w6 = str_replace(w6,'^_',''))

at2 = at2 %>%
  select(w1:w6,everything(),-(`#OTU ID`))


bac_filt = as_tibble(fread('review/HMP2/16S/taxonomic_profiles_filtr.tab'))
bac_filt = bac_filt %>%
  filter(group == 'yes')
bac_filt

quart_bac = at2 %>%
  left_join(bac_filt, by=c('w5','w6')) %>%
  filter(group == 'yes') %>%
  select(w5,w7,everything(),-w1,-w2,-w3,-w4,-w6) %>%
  distinct() %>%
  group_by(w5,w7) %>%
  summarise(across(where(is.numeric),max)) %>%
  ungroup() %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  distinct(count) %>%
  summarise(quantile(count)) %>%
  pull()
quart_bac

bac_16S = at2 %>%
  left_join(bac_filt, by=c('w5','w6')) %>%
  filter(group == 'yes') %>%
  select(w5,w7,everything(),-w1,-w2,-w3,-w4,-w6) %>%
  distinct() %>%
  group_by(w5,w7) %>%
  summarise(across(where(is.numeric),max)) %>%
  ungroup() %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(w5,w7) %>%
  summarize(n_16S = n(),
            n_more0_16S = length(count[count > 0]),
            n_more0_16S_per = round((n_more0_16S/n_16S)*100,2),
            mean_more0_16S = mean(count[count > 0])) %>%
  ungroup() %>%
  mutate(Bac = paste(w5,w7,sep='_')) %>%
  mutate(group_16S = case_when(
    between(mean_more0_16S, quart_bac[1], quart_bac[2]) ~ "1",
    between(mean_more0_16S, quart_bac[2], quart_bac[3]) ~ "2",
    between(mean_more0_16S, quart_bac[3], quart_bac[4]) ~ "3",
    between(mean_more0_16S, quart_bac[4], quart_bac[5]) ~ "4")) %>%
  mutate(Bac = str_to_title(Bac)) %>%
  select(Bac, n_16S, n_more0_16S, n_more0_16S_per, mean_more0_16S, group_16S)
bac_16S

##### ecs (bac) #####

ecs = fread('review/HMP2/Metatranscriptomes/ecs.tsv')

ecs_bac = as_tibble(ecs) %>%
  mutate(EC = str_remove(`# Gene Family`,":.*")) %>%
  group_by(EC) %>%
  slice(-1) %>%
  ungroup() %>%
  mutate(Bac = str_remove(`# Gene Family`,".*\\|")) %>%
  select(-`# Gene Family`) %>%
  distinct() %>%
  mutate(test = str_count(Bac,'s__')) %>%
  filter(test == 1) %>%
  mutate(Bacteria = str_remove(Bac,'.*s__')) %>%
  select(Bacteria,EC,everything(),-Bac,-test)
ecs_bac

length(unique(ecs_bac$Bacteria))
length(unique(ecs_bac$EC))

quart_bac_enz = ecs_bac %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  distinct(count) %>%
  summarise(quantile(count)) %>%
  pull()
quart_bac_enz

bac_enz = ecs_bac %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(Bacteria,EC) %>%
  summarize(n = n(),
            n_more0 = length(count[count > 0]),
            n_per = round((n_more0/n)*100,2),
            median = median(count),
            median_more0 = median(count[count > 0])) %>%
  mutate(group_enz_tox = case_when(
    between(median_more0, quart_bac_enz[1], quart_bac_enz[2]) ~ "1",
    between(median_more0, quart_bac_enz[2], quart_bac_enz[3]) ~ "2",
    between(median_more0, quart_bac_enz[3],quart_bac_enz[4]) ~ "3",
    between(median_more0, quart_bac_enz[4], quart_bac_enz[5]) ~ "4")) %>%
  left_join(prot2,by = c('EC'))
bac_enz

bac_byallenz = ecs_bac %>%
  group_by(Bacteria) %>%
  summarise(n_enz = n(),
            across(where(is.numeric),sum)) %>%
  ungroup() %>%
  mutate(n_enz = as.character(n_enz)) %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(Bacteria, n_enz) %>%
  summarise(n = n(),
            n_more0 = length(count[count > 0]),
            n_per = round((n_more0/n)*100,2),
            min = min(count),
            min_more0 = min(count[count > 0]),
            max = max(count),
            max_more0 = max(count[count > 0]),
            mean = mean(count),
            mean_more0 = mean(count[count > 0]),
            median = median(count),
            median_more0 = median(count[count > 0]),
            sd = sd(count),
            sd_more0 = sd(count[count > 0]),)

bac_byallenz  

bac_enz_tox = ecs_bac %>%
  filter(EC %in% unique(tox_enz2$ENZYME)) %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(Bacteria,EC) %>%
  summarize(n = n(),
            n_more0 = length(count[count > 0]),
            n_per = round((n_more0/n)*100,2),
            median = median(count),
            median_more0 = median(count[count > 0])) %>%
  mutate(group_enz_tox = case_when(
    between(median_more0, quart_bac_enz[1], quart_bac_enz[2]) ~ "1",
    between(median_more0, quart_bac_enz[2], quart_bac_enz[3]) ~ "2",
    between(median_more0, quart_bac_enz[3],quart_bac_enz[4]) ~ "3",
    between(median_more0, quart_bac_enz[4], quart_bac_enz[5]) ~ "4")) %>%
  left_join(prot2,by = c('EC'))
bac_enz_tox

bac_byuroenz = ecs_bac %>%
  filter(EC %in% unique(tox_enz2$ENZYME)) %>%
  group_by(Bacteria) %>%
  summarise(n_enz_uro = n(),
            across(where(is.numeric),sum)) %>%
  ungroup() %>%
  mutate(n_enz_uro = as.character(n_enz_uro)) %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(Bacteria, n_enz_uro) %>%
  summarise(n_uro = n(),
            n_more0_uro = length(count[count > 0]),
            n_per_uro = round((n_more0_uro/n_uro)*100,2),
            min_uro = min(count),
            min_more0_uro = min(count[count > 0]),
            max_uro = max(count),
            max_more0_uro = max(count[count > 0]),
            mean_uro = mean(count),
            mean_more0_uro = mean(count[count > 0]),
            median_uro = median(count),
            median_more0_uro = median(count[count > 0]),
            sd_uro = sd(count),
            sd_more0_uro = sd(count[count > 0]),)

bac_byuroenz  

bac_tox = ecs_bac %>%
  select(Bacteria, EC) %>%
  inner_join(tox_enz2, by=c('EC' = 'ENZYME')) %>%
  select(-EC) %>%
  distinct() %>%
  group_by(Bacteria) %>%
  mutate(tox = paste0(TOXIN, '_',TOXIN_Name,sep='', collapse = '; '),
         tox2 = paste0(TOXIN_Name,sep='', collapse = '; ')) %>%
  ungroup() %>%
  select(Bacteria, tox, tox2) %>%
  mutate(n_tox = str_count(tox,';')+1) %>%
  distinct()
bac_tox 


bac_byalluro = bac_byallenz %>%
  left_join(bac_byuroenz, by='Bacteria') %>%
  left_join(bac_tox, by='Bacteria')
bac_byalluro

write.table(bac_byalluro,file='review/tables/from_bac.tab',quote=T,row.names=F,col.names=T,sep='\t')
write.table(bac_byalluro,file='/mnt/public_html/suvorova/toxins/review/from_bac.tab',quote=T,row.names=F,col.names=T,sep='\t')

#####

bac_byalluro2 = bac_byalluro %>%
  mutate(n_per2 = paste(n_per,"%",sep=''),
         Bacteria = str_replace_all(Bacteria,'_',' ')) %>%
  select(Bacteria, n_per2, n_tox, n_enz_uro, tox2)

write.table(bac_byalluro2,file='review/tables/from_bac_topap.tab',quote=F,row.names=F,col.names=T,sep='\t')
write.table(bac_byalluro2,file='/mnt/public_html/suvorova/toxins/review/from_bac_topap.tab',quote=F,row.names=F,col.names=T,sep='\t')

  



#### by_tox

EC_name = ecs2 %>%
  select(EC, Name_enz) %>%
  distinct()
EC_name
  
from_enz = ecs_bac %>%
  filter(EC %in% unique(tox_enz2$ENZYME)) %>%
  pivot_longer(where(is.numeric), names_to = "income", values_to = "count") %>%
  group_by(EC) %>%
  summarize(n_bac = length(unique(Bacteria)),
            max = max(count[count > 0]),
            median_more0 = median(count[count > 0])) %>%
  mutate(group_enz_tox = case_when(
    between(median_more0, quart_bac_enz[1], quart_bac_enz[2]) ~ "1",
    between(median_more0, quart_bac_enz[2], quart_bac_enz[3]) ~ "2",
    between(median_more0, quart_bac_enz[3],quart_bac_enz[4]) ~ "3",
    between(median_more0, quart_bac_enz[4], quart_bac_enz[5]) ~ "4")) %>%
  left_join(prot2,by = c('EC')) %>%
  left_join(tox_enz2,by=c('EC' = 'ENZYME')) %>%
  left_join(EC_name, by='EC') %>%
  group_by(EC, Name_enz, n_bac, max, median_more0, group_enz_tox, prot) %>%
  mutate(tox = paste0(TOXIN, '_',TOXIN_Name,sep='', collapse = '; '),
         tox2 = paste0(TOXIN_Name,sep='', collapse = '; ')) %>%
  select(-TOXIN, -TOXIN_Name) %>%
  distinct() %>%
  ungroup()

table(from_enz$group_enz_tox)

write.table(from_enz,file='review/tables/from_enz.tab',quote=T,row.names=F,col.names=T,sep='\t')
write.table(from_enz,file='/mnt/public_html/suvorova/toxins/review/from_enz.tab',quote=T,row.names=F,col.names=T,sep='\t')

from_enz2 = from_enz %>%
  left_join(a_stat, by = c('EC')) %>%
  mutate(prot2 = str_remove(prot,'.*_')) %>%
  mutate(n_per_enz = paste(n_per_enz,'%',sep=' '),
         prot2 = paste(prot2,'%',sep=' '),
         prot2 = str_replace(prot2,'NA %','NA')) %>%
  select(Name_EC, n_bac, group_enz_tox, n_per_enz,prot2,tox2)
from_enz2

write.table(from_enz2,file='review/tables/from_enz_topap.tab',quote=F,row.names=F,col.names=T,sep='\t')
write.table(from_enz2,file='/mnt/public_html/suvorova/toxins/review/from_enz_topap.tab',quote=F,row.names=F,col.names=T,sep='\t')


#### pictures

library(tidyverse)
library(data.table)

bac = fread('review/tables/from_bac.tab')
head(bac)

ggplot(bac) +
  aes(x = log(n_enz_uro), y = n_per_uro) +
  geom_point()






