quantiles <- projec %>% select(time, sim, starts_with('cases')) %>% select(-CasesAll) %>%
  mutate(ReportedAll = rowSums(.[grep("cases", names(.))], na.rm = TRUE)) %>%
  select(time, sim,  ReportedAll) %>%
  as_tibble() %>% 
  mutate(isdata = sim == "data") %>%
  gather(variable, value, -time, -sim, -isdata) %>% 
  group_by(time, isdata, variable) %>% 
  summarise( q05 = quantile(value, 0.025, na.rm = T),
             mean = mean(value, na.rm = T),
             q50 = quantile(value, 0.5, na.rm = T),
             q95 = quantile(value, 0.975, na.rm = T)) %>%
  ungroup() %>%
  mutate(isdata = ifelse(isdata, "data", "simulation"),
         date = yearsToDateTime(time)) %>% 
  #filter(date <= yearsToDate(2019.5)) %>%
  mutate(date = as.Date(round_date(date)))

p.sim <- ggplot(data = quantiles  %>% filter(isdata == 'simulation'),
                aes(x = date)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.1) +
  geom_line(aes(y = q50)) + 
  geom_point(data= quantiles  %>% filter(isdata == 'data'), aes(x = date, y = q50))
print(p.sim)


quantiles  %>% filter(isdata == 'simulation') %>% select(date, q05, q50,  q95) %>% write.csv('test.csv')