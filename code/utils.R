# read in individual year data from https://www.nber.org/research/data/mortality-data-vital-statistics-nchs-multiple-cause-death-data
read_nchs_data <- function(folder, COD, icd8_codes, icd9_codes, icd10_codes) {
  
  pattern8 <- paste0("^(", paste(icd8_codes, collapse = "|"), ")")
  pattern9 <- paste0("^(", paste(icd9_codes, collapse = "|"), ")")
  pattern10 <- paste0("^(", paste(icd10_codes, collapse = "|"), ")")
  
  csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  
  # list to hold years individually
  dt_list <- vector("list", length(csv_files))
  
  for (i in seq_along(csv_files)) {
    file <- csv_files[i]
    file_name <- tools::file_path_sans_ext(basename(file))
    yr <- as.numeric(str_extract(file_name, "\\d{4}"))
    
    dt <- suppressWarnings(fread(file, select = c("age", "monthdth", "ucod", paste0("record_", 1:20))))
    
    dt[, year := yr]
    dt[, date := as.Date(paste0(year, "-", monthdth, "-01"))]
    dt <- dt[date >= as.Date("1968-09-01") & date <= as.Date("2021-08-01") & !grepl("999$", age)]
    
    if (COD == "flu MCOD") { # all 20 CODs, not just UCOD
      # unite `record_` columns into single "cod" column to search through. definitely not the most efficient process
      record_cols <- grep("^record_", names(dt), value = TRUE)
      dt[, cod := do.call(paste, c(.SD, sep = ", ")), .SDcols = record_cols]
      dt[, (record_cols) := NULL]  # remove original record_ columns
      
      # create death variable based on ICD patterns
      dt[, death := fifelse(year %in% 1968:1978 & (grepl(pattern8, ucod) | grepl(pattern8, cod)), 1,
                            fifelse(year %in% 1979:1998 & (grepl(pattern9, ucod) | grepl(pattern9, cod)), 1,
                                    fifelse(year %in% 1999:2021 & (grepl(pattern10, ucod) | grepl(pattern10, cod)), 1, 0)))]
      
    } else if (COD == "all cause") {
      dt[, death := 1]
    } else {
      dt[, death := fifelse(year %in% 1968:1978 & grepl(pattern8, ucod), 1,
                            fifelse(year %in% 1979:1998 & grepl(pattern9, ucod), 1,
                                    fifelse(year %in% 1999:2022 & grepl(pattern10, ucod), 1, 0)))]
    }
    
    dt <- dt[death == 1]
    
    dt[, age_unit := substr(age, 1, 1)] # age unit to help calculate birth year
    dt[, age_value := suppressWarnings(as.numeric(substr(age, 2, 3)))] # pre-2003
    dt[, age_value_03 := suppressWarnings(as.numeric(substr(age, 3, 4)))] # 4-digit format starting in 2003
    
    dt[, age_num := suppressWarnings(as.numeric(age))]
    
    dt[, age_recode := fifelse(
      year >= 2003 & age_unit == "1", age_num - 1000,
      fifelse(year >= 2003 & age_unit %in% c("2", "4", "5", "6"), 0,
              fifelse(year < 2003 & age_num < 200, age_num,
                      fifelse(year < 2003 & age_unit %in% c("2", "3", "4", "5", "6"), 0, NA_real_))))]
    
    # birth year
    dt[, birth_year := fifelse(age_recode > 0, year - age_recode, # older than 1 year
                               fifelse(age_unit %in% c("3", "4", "5", "6"), year, # less than a month old
                                       fifelse(year >= 2003 & age_unit == "2" & age_value_03 > monthdth, year - 1, # < 1 yr, born year prior
                                               fifelse(year >= 2003 & age_unit == "2" & age_value_03 <= monthdth, year, # < 1 yr, born this year
                                                       fifelse(year < 2003 & age_unit == "2" & age_value > monthdth, year - 1, # < 1 yr, born year prior
                                                               fifelse(year < 2003 & age_unit == "2" & age_value <= monthdth, year, # < 1 yr, born this year
                                                                       NA_integer_))))))]
    
    dt[, season := fifelse(year == 1968, "1968-1969",
                           paste0(ifelse(monthdth >= 9, year, year - 1), "-", ifelse(monthdth >= 9, year + 1, year)))]
    
    dt[, weight := fifelse(year == 1972, 2, 1)] # weight for 1972 50% sample
    
    dt_list[[file_name]] <- dt[, .(year, season, age_recode, birth_year, weight)]
    
    print(yr) # check progress
  }
  
  return(dt_list)
}


# turn long dataframe of exposures or mortality into a matrix for mortality models. years need not be specified unless splitting
# the data up
matrix_fun <- function(df, col, years, ages) {
  if(missing(years)) {
    years <- 1968:2020
  }
  
  if(missing(ages)) {
    ages <- 0:108
  }
  
  X <- matrix(NA, nrow = length(ages), ncol = length(years), dimnames = list(ages, years))
  
  # match ages and years
  age_index <- match(df$age_recode, ages)
  year_index <- match(df$year_recode, years)
  
  # fill values
  X[cbind(age_index, year_index)] <- df[[col]]
  
  return(X)
}




# return fitted mortality rates from fitted StMoMo 
fitted_values_long <- function(fitted_model, type, col_name) {
  if(missing(type)) {
    type = "rates"
  }
  if(missing(col_name)) {
    col_name <- "fitted_mx"
  }
  
  as.data.frame(fitted(fitted_model, type = type)) %>% 
    tibble::rownames_to_column("age_recode") %>% 
    pivot_longer(cols = -age_recode,
                 names_to = "year_recode",
                 values_to = col_name) %>% 
    mutate_if(is.character, as.numeric)
}



residuals_long <- function(fitted_model, col_name) {
  if(missing(col_name)) {
    col_name <- "residuals"
  }
  as.data.frame(residuals(fitted_model)[["residuals"]]) %>%
    tibble::rownames_to_column("age_recode") %>% 
    pivot_longer(cols = -age_recode,
                 names_to = "year_recode",
                 values_to = col_name) %>%
    mutate_if(is.character, as.numeric)
}


# objects to fit a Plat model. from StMoMo documentation
f2 <- function(x, ages) mean(ages) - x

constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages) { 
  nYears <- dim(wxt)[2]
  x <- ages 
  t <- 1:nYears 
  c <- (1 - tail(ages, 1)):(nYears - ages[1]) 
  xbar <- mean(x)
  phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit) 
  phi <- coef(phiReg) 
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c ^ 2 
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t 
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t ^ 2 - 2 * xbar * t) 
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x ^ 2 
  ci <- rowMeans(kt, na.rm = TRUE) 
  ax <- ax + ci[1] + ci[2] * (xbar - x) 
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2] 
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc) 
}

PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE, 
               periodAgeFun = c("1", f2), cohortAgeFun = "1", constFun = constPlat)



# incidence rate ratio table
rr_tab <- function(mod) {
  cov.mod <- vcovHC(mod, type="HC0")
  
  std.err <- sqrt(diag(cov.mod))
  
  r.est <- na.omit(cbind(estimate = coef(mod, complete = FALSE), 
                         "robust SE" = std.err,
                         "Pr(>|z|)" = 2 * pnorm(abs(coef(mod, complete = FALSE)/std.err), lower.tail = FALSE),
                         LL = coef(mod, complete = FALSE) - 1.96 * std.err,
                         UL = coef(mod, complete = FALSE) + 1.96 * std.err))
  
  deltas <- lapply(1:length(coef(mod, complete = FALSE)), function(i) as.formula(paste0("~ exp(x", i, ")")))
  
  s <- deltamethod(deltas, coef(mod, complete = FALSE), cov.mod)
  
  # exponentiate estimates and confidence intervals, dropping p values
  rexp.est <- exp(r.est[, -3])
  
  # Replace SEs with delta method estimates
  rexp.est[, "robust SE"] <- s
  
  return(round(rexp.est, 3))
}




# plotting functions

plot_residuals <- function(df, resids, limits) {
  if(missing(limits)) {
    limits <- c(-5.5, 5.5)
  }
  ggplot(df, aes(x = year_recode, y = age_recode, fill = resids)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "red", midpoint = 0, 
                       na.value = "transparent", name = "Excess \nmortality", limits = limits) +
  labs(x = "Influenza season", y = "Age", tag = "Lower \nmortality") +
  theme_minimal() + 
  theme(plot.tag.position = c(0.8875, 0.375),
        plot.tag = element_text(hjust = 0, size = 15),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15, family = "Helvetica"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) + 
  scale_x_continuous(expand = c(0, 0), breaks = seq(1970, 2020, by = 5), 
                     labels = c("1970", "1975", "1980", "1985", "1990", "1995", "2000",
                                "2005",  "2010", "2015", "2020")) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, by = 10))
}





#grid of all seasons for supplementary text
arrange_plots <- function(df, years, sigma, breaks, ylim, file) {
  df_filtered <- df %>% filter(year_recode %in% years)
  
  plot_list <- lapply(seq_along(years), function(idx) {
    i <- years[idx]
    df_year <- df_filtered %>% filter(year_recode == i)
    
    y_label <- ifelse(idx %% 3 == 1, "Influenza deaths per 100,000", "")  # only on the left side
    
    ggplot(df_year, na.rm = TRUE) +
      geom_ribbon(aes(x = birth_year, ymin = lower_flu_mx_p100k, ymax = upper_flu_mx_p100k), 
                  fill = "grey70") +
      geom_line(aes(x = birth_year, y = fitted_flu_mx_p100k), na.rm = TRUE) +
      geom_point(aes(x = birth_year, y = obs_flu_mx_p100k, color = imprinted_strain), 
                 alpha = 0.6, size = 2.25, na.rm = TRUE) +
      scale_color_manual(name = "Imprinted strain", 
                         values = c(pH1N1 = "#ff6698", H1N1_alpha = "#b70000", H1N1_beta = "#ff8900", 
                                    H1N1_gamma = "#ECC905", H2N2 = "#95C90F", H3N2 = "#3683C3", 
                                    "mixed/naive" = "#989898", unknown = "#555"),
                         labels = c("mixed/naive" = "Mixed/Naïve", unknown = "Unknown",
                                    pH1N1 = bquote(pH1N1["\u03b1"]), H1N1_alpha = bquote(H1N1["\u03b1"]),
                                    H1N1_beta = bquote(H1N1["\u03b2"]), H1N1_gamma = bquote(H1N1["\u03b3"])),
                         guide = guide_legend(ncol = 1)) +
      labs(title = ifelse(i != 2008, paste0("‎ ", df_year$season[1], ", ", df_year$circulating_strain[1]), 
                          "‎ 2008-2009, H1N1*"),
           x = "Birth year", y = y_label) +
      scale_x_reverse(breaks = seq(i - (i %% 5), 1860, by = -20),
                      sec.axis = sec_axis(transform = ~.,
                                          breaks = seq(i, i-108, by = -20),
                                          labels = c(0, 20, 40, 60, 80, 100),
                                          name = "Age")) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = sigma),
                         breaks = breaks) +
      theme_bw() +
      theme(text = element_text(size=16, family = "Helvetica"),
            plot.title = element_text(vjust = -16, size = 17),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.title=element_text(size=14),
            legend.key.size = unit(0.7, "cm"),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 18),
            #legend.box.spacing = unit(0, "cm"), 
            plot.margin = unit(c(-0.75,0.15,0.2,-0.6), "cm")) +
      guides(color = guide_legend(override.aes = list(size=2.5))) +
      coord_cartesian(ylim = ylim)
  })
  
  get_legend <- function(p) {
    g <- ggplotGrob(p)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  }
  
  legend <- get_legend(plot_list[[3]] + theme(legend.position = "right"))
  
  plot_list <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
  
  grid <- do.call(cowplot::plot_grid, c(plot_list, list(legend, nrow = 3, ncol = 3, align = "hv", axis = "tblr")))
  
  fig <- ggdraw() + draw_plot(grid, x = 0.03, y = 0, width = 0.97, height = 1)
  
  ggsave(filename = file, plot = fig, height = 10, width = 10, bg='#ffffff', dpi = 800)
}



fitted_plot <- function(data, year, x_breaks, sec_axis_name = "", xlab = "", ylab = "", 
                        label_text, label_x, margin = c(0, 0, -0.5, 0)) {
  data %>%
    filter(year_recode == year) %>%
    ggplot() +
    geom_ribbon(aes(x = birth_year, ymin = lower_flu_mx_p100k, ymax = upper_flu_mx_p100k), 
                fill = "gray75", alpha = 0.8) +
    geom_line(aes(x = birth_year, y = fitted_flu_mx_p100k), na.rm = TRUE) +
    geom_point(aes(x = birth_year, y = obs_flu_mx_p100k, color = imprinted_strain), 
               alpha = 0.6, size = 2.25, na.rm = TRUE) +
    scale_color_manual(name = "Imprinted strain", 
                       values = c(unknown = "gray35", pH1N1 = "#ff6698", H1N1_alpha = "#b70000", 
                                  H1N1_beta = "#ff8900", H1N1_gamma = "#ECC905", H2N2 = "#95C90F", 
                                  H3N2 = "#3683C3", "mixed/naive" = "gray60"),
      labels = c("unknown" = "Unknown", "mixed/naive" = "Mixed/Naïve", pH1N1 = bquote(pH1N1["\u03b1"]), 
                 H1N1_alpha = bquote(H1N1["\u03b1"]), H1N1_beta = bquote(H1N1["\u03b2"]), 
                 H1N1_gamma = bquote(H1N1["\u03b3"]))) +
    scale_x_reverse(breaks = x_breaks,
                    sec.axis = sec_axis(transform = ~.,
                                        breaks = seq(year, year - 108, by = -20),
                                        labels = c(0, 20, 40, 60, 80, 100),
                                        name = sec_axis_name)) +
    scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = 0.01),
                       breaks = c(0, 0.1, 1, 10, 100),
                       labels = c(0, 0.1, 1, 10, 100)) +
    coord_cartesian(ylim = c(0, 500)) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(text = element_text(size = 18, family = "Helvetica"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          plot.margin = unit(margin, "lines"), 
          axis.title = element_text(size = 14)) +
    annotate("label", x = label_x, y = 90, label = label_text, size = 4.5)
}
