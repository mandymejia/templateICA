GroupICATv4.0b Updates (May 01, 2017):

spm mexmaci64 files are updated to work on R2017a. 

GroupICATv4.0b Updates (April 28, 2017):

PCA paper reference S. Rachakonda, R. F. Silva, J. Liu and V. D. Calhoun, "Memory Efficient PCA Methods for Large Group ICA", Frontiers in Neuroscience, 2016 is mentioned in the 
icatb/icatb_analysis_functions/icatb_calculate_pca.m and icatb/icatb_parallel_files/icatb_parCalculatePCA.m files.

GroupICATv4.0b Updates (April 17, 2017):

icatb/icatb_display_functions/icatb_plot_connectogram.m file is updated to include:

	1. By default, only 12 network colors are provided instead of 14. Line number 244 is fixed.
	2. An option is provided to set connectogram figure size in icatb_defaults.m file. Please see CONNECTOGRAM_FIG_POS variable.
	3. Optional flag is provided not to reorder the connectivity matrix when all the components are selected.

GroupICATv4.0b Updates (April 13, 2017):

1. icatb/icatb_io_data_functions/icatb_read_gzip_nii.m is fixed to load data correctly when slices are entered as an input in the function.
2. Components (and its labels) were not updated in the connectogram after removing zeros in the FNC matrix (April 11th update). icatb/icatb_display_functions/icatb_plot_connectogram function 
is now fixed. 


GroupICATv4.0b Updates (April 12, 2017):

Network color is set to a default value if uisetcolor window is closed without selecting a value. Function icatb/icatb_display_functions/icatb_plot_connectogram.m is updated.

GroupICATv4.0b Updates (April 11, 2017):

1. GIFT now reads NIFTI Gzip (*.nii.gz) format. For very large files (multiband data) it is recommended to set NIFTI_GZ = 1 in icatb_defaults.m (files will be un-archived) for reading 
Gzip files faster.
2. Despike utility computational performance (speed) is improved.
3. We now provide various options to estimate clusters using elbow, BIC, AIC and Dunns index. Optimal clusters is also provided as a stand alone tool.
4. Connectogram is now added as a stand alone tool in GIFT -> Tools -> Display Tools -> Misc menu. Options to use rendered images, user defined RGB images and component labels are provided 
in the connectogram function. 
5. Generate mask utility is now added in the SBM toolbox. To access the tool, use SBM -> Utilities -> Generate Mask.
6. Deprecated legend parameter is removed in icatb_orthoViewer.m.

The following files are added or updated:	

	1. icatb/dfnc_toolbox.fig      
	2. icatb/dfnc_toolbox.m        
	3. icatb/gift.fig              
	4. icatb/icatb_defaults.m      
	5. icatb/icatb_displayGUI.m    
	6. icatb/icatb_setup_analysis.m
	7. icatb/sbm.fig   
	8. icatb/icatb_analysis_functions/icatb_parameterInitialization.m
	9. icatb/icatb_analysis_functions/icatb_runAnalysis.m   
	10. icatb/icatb_batch_files/icatb_read_batch_file.m
	11. icatb/icatb_display_functions/icatb_display_composite.m
	12. icatb/icatb_display_functions/icatb_orthoViewer.m
	13. icatb/icatb_display_functions/icatb_plot_connectogram.m
	14. icatb/icatb_helper_functions/icatb_createMask.m         
	15. icatb/icatb_helper_functions/icatb_despike.m            
	16. icatb/icatb_helper_functions/icatb_despike_tc.m         
	17. icatb/icatb_helper_functions/icatb_generateMask.m       
	18. icatb/icatb_helper_functions/icatb_get_countTimePoints.m
	19. icatb/icatb_helper_functions/icatb_loadAndInterpTC.m    
	20. icatb/icatb_helper_functions/icatb_optimal_clusters.m   
	21. icatb/icatb_helper_functions/icatb_post_process_dfnc.m  
	22. icatb/icatb_helper_functions/icatb_read_hdr.m           
	23. icatb/icatb_helper_functions/icatb_rename_4d_file.m     
	24. icatb/icatb_helper_functions/icatb_returnHInfo.m        
	25. icatb/icatb_helper_functions/icatb_update_mask.m        
	26. icatb/icatb_helper_functions/icatb_utilities.m  
	27. icatb/icatb_io_data_functions/icatb_dataSelection.m
	28. icatb/icatb_io_data_functions/icatb_gz.jar         
	29. icatb/icatb_io_data_functions/icatb_inputdlg2.m    
	30. icatb/icatb_io_data_functions/icatb_read_data.m    
	31. icatb/icatb_io_data_functions/icatb_read_gzip_nii.m
	32. icatb/icatb_mancovan_files/icatb_mysplinefun.m 
	33. icatb/icatb_mancovan_files/icatb_run_mancovan.m
	

GroupICATv4.0b Updates (March 01, 2017):

icatb/icatb_mancovan_files/icatb_run_mancovan.m file is fixed to handle NaNs in the nuisance covariates.