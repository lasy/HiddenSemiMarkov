

compute_dtemp = function(E){
  colnames_to_keep = colnames(E)

  E = E %>%
    dplyr::group_by(user_id) %>%
    dplyr::mutate(mean_temp_7d = data.table::frollmean(temp, n = 7, align = "right", na.rm = TRUE, hasNA = TRUE),
                  mean_temp_7d_lagged = dplyr::lag(mean_temp_7d, 13),
                  dtemp_num = mean_temp_7d - mean_temp_7d_lagged,
                  dtemp_char = dplyr::case_when(
                    is.na(dtemp_num) ~ "missing",
                    dtemp_num >= 0.5 ~ "increase",
                    dtemp_num <= -0.5 ~ "decrease",
                    TRUE ~ "no_change"
                  ),
                  dtemp = dtemp_char %>% factor(.,levels = c("missing","increase","no_change","decrease"))
    ) %>%
    dplyr::select(tidyselect::all_of(colnames_to_keep), dtemp) %>%
    dplyr::ungroup()

  E
}



compute_bleeding_score = function(E){

  bleeding_dict = data.frame(code =  c("none","spotting" ,"light"  ,"medium"   ,"heavy"),
                             score = c(0,0.25,0.5,0.75,1))
  E = E  %>%
    dplyr::mutate(bleeding_score = bleeding_dict$score[match(bleeding, bleeding_dict$code)] %>%  tidyr::replace_na(0))
  E
}



compute_f.bleeding.last.week= function(E){

  colnames_to_keep = colnames(E)

  E = compute_bleeding_score(E)

  E = E  %>%
    dplyr::group_by(user_id) %>%
    dplyr::mutate(f.bleeding.last.week =
                    data.table::frollmean(bleeding_score>0, n = 7, align = "right", fill = 0, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(tidyselect::all_of(colnames_to_keep), f.bleeding.last.week)

  E
}




FAM_data_augmentation = function(X, features = c("dtemp","f.bleeding.last.week"), get_var_names_only = FALSE, get_var_types_only = FALSE){

  if(get_var_names_only) return(features)

  if(get_var_types_only) return(
    list(
      dtemp = list(
        type = "cat",
        values = c("missing","increase","no_change","decrease")
      ),
      f.bleeding.last.week = list(
        type = "continuous"
      )#,
      # bleeding_density_5d = list(
      #   type = "continuous"
      # )
    )
  )

  # the functions for each of these features are in "Scripts/00_functions_feature_engineering.R"
  ordered_features = features

  if("seq_id" %in% colnames(X)){X = X %>% dplyr::rename(user_id = seq_id, rel_date = t)}

  if(any(diff(order(X$user_id, X$rel_date)) != 1)) stop("X must be arranged by seq_id (user_id) and t (rel_date)")

  E = X
  for(f in ordered_features){
    fun = eval(parse(text = stringr::str_c("compute_",f)))
    E = fun(E = E)
  }

  if(any(E$user_id != X$user_id) | any(E$rel_date != X$rel_date)) stop("Order of E is not the same as order of X")

  E = E %>%  dplyr::select(all_of(features)) %>%  dplyr::ungroup()
  E
}


library(magrittr)
library(tictoc)

load(file = paste0("../../Papers/semiM/semiM-Data/Kindara/tmp_data/FAM_full_spec.Rdata"), verbose = TRUE)

# tic()
# FAM_full_init = initialize_hsmm(model = FAM_full_spec, verbose = TRUE)
# toc()

# save(FAM_full_init, file = "FAM_full_init.Rdata")
load(file = paste0("FAM_full_init.Rdata"), verbose = TRUE)

Xsim = purrr::map_dfr(
  .x = 1:6,
  .f = function(i) simulate_hsmm(model = FAM_full_spec, n_state_transitions = 100,
                                 seed = i, seq_id = stringr::str_c("sim_seq_",i))
)

tic()
FAM_full_fitted = fit_hsmm(model = FAM_full_init, X = Xsim, n_iter = 10, lock.transition = TRUE, lock.sojourn = TRUE, verbose = TRUE)
toc()





