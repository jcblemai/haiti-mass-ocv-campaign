args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  # default departement
  args[1] = "Artibonite"
}
departement <- args[1]

load(paste0("sirb_cholera_pomped_", departement, ".rda")
     )


coef(sirb_cholera)["mu_B"] <- 365/5
coef(sirb_cholera)["betaB"] <- 0.1
#coef(sirb_cholera)["Rtot_0"] <-100000000

p <- simulate(sirb_cholera, nsim = 10, as.data.frame = T) %>% 
  gather(variable, value, -time, -rain, -sim)  %>% 
  ggplot(aes(x = time, y = value, color = sim)) + 
  #  geom_line(aes(y = cases), color = datacol, lwd = 0.2) 
  geom_line() + 
  facet_wrap(~variable, scales = "free_y") 


datacol <- "#ED0000"
print(p)