test <- sims_stoch %>% filter(.id %in% id_sims, var == "CC") %>%
  mutate(loc = factor(loc, levels = c("Lübeck", "Bogota"), labels = c("Lübeck", "Bogota"))) %>%
  filter(loc == "Lübeck")
test_det <- sims_det %>% 
  mutate(loc = factor(loc, levels = c("Lübeck", "Bogota"), labels = c("Lübeck", "Bogota"))) %>%
  filter(loc == "Lübeck")

pl_up1 <- ggplot(data = test, 
                mapping = aes(x = week, y = 1e2 * value / N_val, group = .id)) + 
  geom_line(color = "grey", alpha = 0.5) + 
  geom_line(data = test_det, mapping = aes(x = week, y = 1e2 * CC / N_val), color = "black") +
  facet_wrap(~ loc, scales = "free_y", nrow = 1) + 
  theme_classic() + 
  scale_x_continuous(expand = c(0,0)) + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank()) + 
  labs(x = "Time (weeks)", y = "Weekly cases (per 100)")
print(pl_up1)

test_main <- pars_main %>%
  mutate(loc = factor(loc, levels = c("Lübeck", "Bogota"), labels = c("Lübeck", "Bogota")))%>%
  filter(loc == "Lübeck")
pl_low1 <- test_main %>%
  group_by(loc, par) %>%
  mutate(m_mle = median(mle)) %>%
  ggplot(mapping = aes(x = mle)) + 
  geom_vline(xintercept = pars_true$true[pars_true$par == "e_Te"], linetype = "dotted") + 
  geom_density(alpha = 0.1, adjust = 2, color = "#f17522", fill = "#f17522") + 
  geom_vline(aes(xintercept = m_mle), color = "#f17522") + 
  geom_rug(alpha = 0.5, color = "#f17522") + 
  scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0,0)) + 
  facet_wrap(par ~. , scales = "free_y", ncol = 1) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1), colour = "black"),
        legend.position = "none") + 
  labs(x = "Parameter estimate", y = "Density", color = "", fill = "")
print(pl_low1)

pl_up1/pl_low1


test <- sims_stoch %>% filter(.id %in% id_sims, var == "CC") %>%
  mutate(loc = factor(loc, levels = c("Lübeck", "Bogota"), labels = c("Lübeck", "Bogota"))) %>%
  filter(loc == "Bogota")
test_det <- sims_det %>% 
  mutate(loc = factor(loc, levels = c("Lübeck", "Bogota"), labels = c("Lübeck", "Bogota"))) %>%
  filter(loc == "Bogota")

pl_up2 <- ggplot(data = test, 
                 mapping = aes(x = week, y = 1e2 * value / N_val, group = .id)) + 
  geom_line(color = "grey", alpha = 0.5) + 
  geom_line(data = test_det, mapping = aes(x = week, y = 1e2 * CC / N_val), color = "black") +
  facet_wrap(~ loc, scales = "free_y", nrow = 1) + 
  scale_x_continuous(expand = c(0,0)) + 
  theme_classic() + 
  theme(strip.background = element_blank(), 
        strip.text = element_blank()) + 
  labs(x = "Time (weeks)", y = "Weekly cases (per 100)")
print(pl_up2)

test_main <- pars_main %>%
  mutate(loc = factor(loc, levels = c("Lübeck", "Bogota"), labels = c("Lübeck", "Bogota")))%>%
  filter(loc == "Bogota")
pl_low2 <- test_main %>%
  group_by(loc, par) %>%
  mutate(m_mle = median(mle)) %>%
  ggplot(mapping = aes(x = mle)) + 
  geom_vline(xintercept = pars_true$true[pars_true$par == "e_Te"], linetype = "dotted") + 
  geom_density(alpha = 0.1, adjust = 2, color = "#7dd2c2", fill = "#7dd2c2") + 
  geom_vline(aes(xintercept = m_mle), color = "#7dd2c2") + 
  geom_rug(alpha = 0.5, color = "#7dd2c2") + 
  scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0,0)) + 
  facet_wrap(par ~. , scales = "free_y", ncol = 1) + 
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = rel(1), colour = "black"),
        legend.position = "none") + 
  labs(x = "Parameter estimate", y = "Density", color = "", fill = "")
print(pl_low2)

pl_up2/pl_low2

(pl_up1/pl_low1) | (pl_up2/pl_low2)





col = c("#0caa8e", "black")

fun_main_pl <- function(loc = "Pasto", legend = T) {
  pl1a <- sim_long %>%
    filter(state_var == "CC") %>%
    filter(R0 == 1.25) %>% 
    filter(loc_tidy_nm == loc) %>%
    ggplot(aes(x = week_no, y = beta_seas, color = Te_effect, group = Te_effect)) + 
    geom_line() + 
    scale_x_continuous("Time (weeks)") + 
    scale_y_continuous("Renormalized \n transmission rate") + 
    scale_color_manual("Effect of temperature", values = col,
                       labels = c(expression(paste("Indirect effect (", delta[Te], " = 0, ", delta[RH], " = -0.2)")),
                                  expression(paste("Total effect (", delta[Te], " = -0.2, ", delta[RH], " = -0.2)")))) +
    facet_wrap(.~loc_tidy_nm, scales = "free", ncol = 2) + 
    theme(legend.position = "top",
          legend.direction = "vertical",
          legend.justification = "left",
          strip.background = element_blank())
  
  if(legend == F) {pl1a <- pl1a + theme(legend.position = "none")}
  
  pl1b <- sim_long %>%
    filter(state_var == "CC") %>%
    filter(R0 == 1.25) %>% 
    mutate(Te_effect = case_when(str_detect(Te_effect, "Indirect") ~ "Indirect effect",
                                 str_detect(Te_effect, "Total") ~ "Total effect")) %>%
    filter(loc_tidy_nm == loc) %>%
    ggplot(aes(x = Te, y = beta_seas, color = Te_effect, group = Te_effect)) + 
    geom_point() + 
    scale_x_continuous("Temperature (°C)") + 
    scale_y_continuous("Renormalized \n transmission rate") + 
    facet_wrap(Te_effect~., scales = "free", ncol = 4) + 
    scale_color_manual("Effect of temperature", values = col) + 
    theme(legend.position = "none",
          strip.background = element_blank())
  
  pl1 <- pl1a / pl1b
  return(pl1) 
}

pl1 <- fun_main_pl(loc = "Lübeck") | fun_main_pl(loc = "Pasto", legend = F)




