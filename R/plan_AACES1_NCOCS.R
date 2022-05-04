plan <- drake_plan(

  clinical_data = data_import(fs::path("","Volumes","Peres_Research")),
  ancestry_data = ancestry_import(fs::path("","Volumes","Peres_Research")),
  tx_data = tx_import(fs::path("","Volumes","Peres_Research")),
  # tx_data_ncocs = tx_import2(fs::path("","Volumes","Peres_Research")),
  #-----------------------------------------------------------------------------------------------------------------
  ROI_tumor_ = roit_import(fs::path("","Volumes","Peres_Research")),
  
  ROI_stroma_ = rois_import(fs::path("","Volumes","Peres_Research")),
  
  ROI_total_ = roi_import(fs::path("","Volumes","Peres_Research")),
  
  # ROI_remove = roir_import(fs::path("","Volumes","Peres_Research")),
  #-----------------------------------------------------------------------------------------------------------------
  TMA1_tumor = tmat_import(fs::path("","Volumes","Peres_Research")),
  TMA1_stroma = tmas_import(fs::path("","Volumes","Peres_Research")),
  TMA1_total = tma17_import(fs::path("","Volumes","Peres_Research")),
  
  TMA2_tumor = tma2t_import(fs::path("","Volumes","Peres_Research")),
  TMA2_stroma = tma2s_import(fs::path("","Volumes","Peres_Research")),
  TMA2_total = tma18_import(fs::path("","Volumes","Peres_Research")),
  #-----------------------------------------------------------------------------------------------------------------
  TMAcases_remove = case_remove_import(fs::path("","Volumes","Peres_Research")),
  ROIcases_remove = caseROI_remove_import(fs::path("","Volumes","Peres_Research")),
  common_ROITMA_IDs = common_ROITMA_import(fs::path("","Volumes","Peres_Research")),
  #-----------------------------------------------------------------------------------------------------------------
  cases_match = match_cases_import(fs::path("","Volumes","Peres_Research")),
  
  # II  ### Data Cleaning
  
  # IIa ### TMA data
  
  ## 1-verify that core removed ARE removed from TMA data  #-------------- All good so not run anymore
  # uid <- paste(unique(TMAremove_tumor[1]), collapse = "|")
  # TMA_tumor <-
  #   TMA_tumor[(!grepl(uid, TMA_tumor$image_tag)), ]
  # TMA_stroma <-
  #   TMA_stroma[(!grepl(uid, TMA_stroma$image_tag)),]
  # 
  # uid <- paste(unique(TMA2remove_tumor[1]), collapse = "|")
  # TMA2_tumor <-
  #   TMA2_tumor[(!grepl(uid, TMA2_tumor$image_tag)),]
  # TMA2_stroma <-
  #   TMA2_stroma[(!grepl(uid, TMA2_stroma$image_tag)),]
  
  
  ## 2-bind TMA together
  TMA_tumor_ = binding(TMA1_tumor,TMA2_tumor),
  TMA_stroma_ = binding(TMA1_stroma,TMA2_stroma), 
  TMA_total_ = binding(TMA1_total,TMA2_total), 
  # Update variable names
  TMA_tumor = var_names(TMA_tumor_),
  TMA_stroma = var_names(TMA_stroma_),
  TMA_total = var_names(TMA_total_),
  ROI_tumor = var_names(ROI_tumor_),
  ROI_stroma = var_names(ROI_stroma_),
  ROI_total = var_names(ROI_total_)

)