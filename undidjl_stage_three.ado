/*------------------------------------*/
/*undidjl_stage_three*/
/*written by Eric Jamieson */
/*version 0.6.4 2025-09-22 */
/*------------------------------------*/
version 14.1


cap program drop undidjl_stage_three
program define undidjl_stage_three, rclass

	syntax , folder(string) /// 
	[agg(string) covariates(string) /// 
	save_csv(string) interpolation(string) ///
	weights(string) seed(int 0) nperm(int 1001)]
	
	// Declare usage of Undid and start up Julia
	jl: using Undid
	jl: using DataFrames
	
	
	// Set seed for RI procedure
	if `seed' == 0 {
		qui jl: seed = rand(1:10000)
	}
	else if `seed' > 0 {
		qui jl: seed = `seed'
	}
	else {
		di as error "seed must be set to a value > 0."
		exit 222
	}
	
	// Allow variables to be passed to Julia
	global folder = subinstr("`folder'", "\", "/", .)

	// Parse agg 
	if "`agg'" == "" | "`agg'" == "g"{
		local agg "g"
		qui jl: agg = "g"
	}
	else if "`agg'" == "gt"{
		qui jl: agg = "gt"
	}
	else if "`agg'" == "silo"{
		qui jl: agg = "silo"
	}
	else if "`agg'" == "sgt" {
		qui jl: agg = "sgt"
	}
	else if "`agg'" == "none" {
		qui jl: agg = "none"
	}
	else {
		disp as error "Please set aggregration to silo, g, gt, sgt, or none."
	}
	
	// Parse covariates
	if "`covariates'" == "TRUE" | 	"`covariates'" == "true" | "`covariates'" == "T" | "`covariates'" == "True" {
		qui jl: covariates = true
	}
	else if "`covariates'" == "" | "`covariates'" == "FALSE" | "`covariates'" == "false" | "`covariates'" == "F" | "`covariates'" == "False" {
		qui jl: covariates = false
	}
	else { 
		display as error "Error: set covariates to true or false or omit the argument (defaults to false)"
	}
	
	// Parse save_all_csvs
	if "`save_csv'" == "TRUE" | "`save_csv'" == "true" | "`save_csv'" == "T" | "`save_csv'" == "True" {
		qui jl: save_all_csvs = true
		di as result "Saving combined_diff_data.csv to " "`c(pwd)'"
	}
	else if "`save_csv'" == "" | "`save_csv'" == "FALSE" | "`save_csv'" == "false" | "`save_csv'" == "F" | "`save_csv'" == "False" {
		qui jl: save_all_csvs = false
	}
	else { 
		display as error "Error: set save_csv to true or false or omit the argument (defaults to false)"
	}
	
	// Parse interpolation
	if "`interpolation'" == "" | "`interpolation'" == "FALSE" | "`interpolation'" == "false" | "`interpolation'" == "F" | "`interpolation'" == "False" {
		qui jl: interpolation = false
	}
	else if "`interpolation'" == "linear_function"{
		qui jl: interpolation = "linear_function"
	}
	else{
		display as error "Error: currently the only supported options for interpolation and extrapolation for missing diff_estimates are: linear_function"
	}
	
	// Parse weights
	if "`weights'" == "" {
		qui jl: weights = "both"
	}
	else {
		qui jl: weights = "`weights'"
	}
	
	qui jl: results = undid_stage_three("$folder", agg = agg, covariates = covariates, save_diff_data = save_all_csvs, interpolation = interpolation, weighting = weights, seed = seed, nperm = `nperm')

	    qui jl: if "att_g" in DataFrames.names(results) ///
                results.labels = string.(results.gvar); ///
			elseif "att_sgt" in DataFrames.names(results) ///
				results.labels = string.(results.sgt); ///
			elseif "att_s" in DataFrames.names(results) ///
				results.labels = string.(results.silos); ///
            elseif "att_gt" in DataFrames.names(results) ///
				results.labels = string.(results.gt) ///
            end		

				
	tempname result_frame
    qui cap frame drop `result_frame'
    qui frame create `result_frame'
    qui frame change `result_frame'
    qui jl use results
	
	qui cap confirm variable labels
	if !_rc {
		qui jl: st_local("rowlabels", join(string.(results.labels), " "))
		qui tostring labels, replace
		local counter = 1
		foreach rowlabel in `rowlabels' {
			qui replace labels = "`rowlabel'" in `counter'
			local counter = `counter' + 1
		}
	}
	

	local found 0
	foreach v in att_g att_gt att_s att_sgt {
		capture confirm variable `v'
		if !_rc {                
			local found 1        
			continue, break      
		}
	}
	if `found' == 0 {                
		local agg "none"
	}
	
	
	local N = _N
	if "`agg'" == "silo" {
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "                                       UnDiD.jl Results                    "
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "Silo                      | " as text "ATT             | SE     | p-val  | JKNIFE SE  | JKNIFE p-val | RI p-val|"
		di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
		
		// Initialize a temporary matrix to store the numeric results
        tempname table_matrix
        local num_rows = _N
        local num_cols = 7
        matrix `table_matrix' = J(`num_rows', `num_cols', .)
		local state_names ""
		
		forvalues i = 1/`N' {
			di as text %-25s "`=labels[`i']'" as text " |" as result %-16.7f att_s[`i'] as text " | " as result  %-7.3f att_s_se[`i'] as text "| " as result %-7.3f att_s_pval[`i'] as text "| " as result  %-11.3f att_s_se_jackknife[`i'] as text "| " as result %-13.3f att_s_jknife_pval[`i'] as text "|" as result %-9.3f ri_pval_att_s[`i'] as text "|"
    
			di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
			
		// Store the state name
        local state_name = labels[`i']
        local state_names `state_names' `state_name'
            
        // Fill the matrix with numeric values
            matrix `table_matrix'[`i', 1] = att_s[`i']
            matrix `table_matrix'[`i', 2] = att_s_se[`i']
            matrix `table_matrix'[`i', 3] = att_s_pval[`i']
            matrix `table_matrix'[`i', 4] = att_s_se_jackknife[`i']
            matrix `table_matrix'[`i', 5] = att_s_jknife_pval[`i']
            matrix `table_matrix'[`i', 6] = ri_pval_att_s[`i']
			matrix `table_matrix'[`i', 7] = weights[`i']
		}
		di as text "Aggregation: " as result "silo"
		
		// Set column names for the matrix
        matrix colnames `table_matrix' = ATT SE pval JKNIFE_SE JKNIFE_pval RI_pval W
        
        // Set row names for the matrix using the state names
        matrix rownames `table_matrix' = `state_names'
        
        // Store the matrix in r()
        return matrix restab = `table_matrix'
        

		
	} 
	else if "`agg'" == "gt" {
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "                                       UnDiD.jl Results                    "
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "gt                        | " as text "ATT             | SE     | p-val  | JKNIFE SE  | JKNIFE p-val | RI p-val|"
		di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
		
		// Initialize a temporary matrix to store the numeric results
        tempname table_matrix
        local num_rows = _N
        local num_cols = 7
        matrix `table_matrix' = J(`num_rows', `num_cols', .)
		
		forvalues i = 1/`N' {
			di as text %-25s "`=labels[`i']'" as text " |" as result %-16.7f att_gt[`i'] as text " | " as result  %-7.3f att_gt_se[`i'] as text "| " as result %-7.3f att_gt_pval[`i'] as text "| " as result  %-11.3f att_gt_se_jackknife[`i'] as text "| " as result %-13.3f att_gt_jknife_pval[`i'] as text "|" as result %-9.3f ri_pval_att_gt[`i'] as text "|"
    
			di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
			
			// Store the gt
            local gt_name = labels[`i']
            local gt_names `gt_names' `gt_name'
			
			// Fill the matrix with numeric values
            matrix `table_matrix'[`i', 1] = att_gt[`i']
            matrix `table_matrix'[`i', 2] = att_gt_se[`i']
            matrix `table_matrix'[`i', 3] = att_gt_pval[`i']
            matrix `table_matrix'[`i', 4] = att_gt_se_jackknife[`i']
            matrix `table_matrix'[`i', 5] = att_gt_jknife_pval[`i']
            matrix `table_matrix'[`i', 6] = ri_pval_att_gt[`i']
			matrix `table_matrix'[`i', 7] = weights[`i']
		}
		// Set column names for the matrix
        matrix colnames `table_matrix' = ATT SE pval JKNIFE_SE JKNIFE_pval RI_pval W
        
        // Set row names for the matrix using the gt names
        matrix rownames `table_matrix' = `gt_names'
        
        // Store the matrix in r()
        return matrix restab = `table_matrix'
		di as text "Aggregation: " as result "gt"

	} 
	else if "`agg'" == "g" {
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "                                       UnDiD.jl Results                    "
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "gvar                      | " as text "ATT             | SE     | p-val  | JKNIFE SE  | JKNIFE p-val | RI p-val|"
		di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
		
		// Initialize a temporary matrix to store the numeric results
        tempname table_matrix
        local num_rows = _N
        local num_cols = 7
        qui matrix `table_matrix' = J(`num_rows', `num_cols', .)
		
		forvalues i = 1/`N' {
			di as text %-25s "`=labels[`i']'" as text " |" as result %-16.7f att_g[`i'] as text " | " as result  %-7.3f att_g_se[`i'] as text "| " as result %-7.3f att_g_pval[`i'] as text "| " as result  %-11.3f att_g_se_jackknife[`i'] as text "| " as result %-13.3f att_g_jknife_pval[`i'] as text "|"  as result %-9.3f ri_pval_att_g[`i'] as text "|"
    
			di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
			
			// Store the gvar
            local g_name = labels[`i']
            local g_names `g_names' `g_name'
			
			// Fill the matrix with numeric values
            matrix `table_matrix'[`i', 1] = att_g[`i']
            matrix `table_matrix'[`i', 2] = att_g_se[`i']
            matrix `table_matrix'[`i', 3] = att_g_pval[`i']
            matrix `table_matrix'[`i', 4] = att_g_se_jackknife[`i']
            matrix `table_matrix'[`i', 5] = att_g_jknife_pval[`i']
            matrix `table_matrix'[`i', 6] = ri_pval_att_g[`i']
			matrix `table_matrix'[`i', 7] = weights[`i']
		}
		// Set column names for the matrix
        matrix colnames `table_matrix' = ATT SE pval JKNIFE_SE JKNIFE_pval RI_pval W
        
        // Set row names for the matrix using the gt names
        matrix rownames `table_matrix' = `g_names'
        
        // Store the matrix in r()
        return matrix restab = `table_matrix'
		di as text "Aggregation: " as result "g"
	}
	else if "`agg'" == "sgt" {
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "                                       UnDiD.jl Results                    "
		di as text "-----------------------------------------------------------------------------------------------------"
		di as text "sgt                       | " as text "ATT             | SE     | p-val  | JKNIFE SE  | JKNIFE p-val | RI p-val|"
		di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
		
		// Initialize a temporary matrix to store the numeric results
        tempname table_matrix
        local num_rows = _N
        local num_cols = 7
        matrix `table_matrix' = J(`num_rows', `num_cols', .)
		
		forvalues i = 1/`N' {
			di as text %-25s "`=labels[`i']'" as text " |" as result %-16.7f att_sgt[`i'] as text " | " as result  %-7.3f att_sgt_se[`i'] as text "| " as result %-7.3f att_sgt_pval[`i'] as text "| " as result  %-11.3f att_sgt_se_jackknife[`i'] as text "| " as result %-13.3f att_sgt_jknife_pval[`i'] as text "|" as result %-9.3f ri_pval_att_sgt[`i'] as text "|"
    
			di as text "--------------------------|-----------------|--------|--------|------------|--------------|---------|"
			
			// Store the gt
            local sgt_name = labels[`i']
            local sgt_names `sgt_names' `sgt_name'
			
			// Fill the matrix with numeric values
            matrix `table_matrix'[`i', 1] = att_sgt[`i']
            matrix `table_matrix'[`i', 2] = att_sgt_se[`i']
            matrix `table_matrix'[`i', 3] = att_sgt_pval[`i']
            matrix `table_matrix'[`i', 4] = att_sgt_se_jackknife[`i']
            matrix `table_matrix'[`i', 5] = att_sgt_jknife_pval[`i']
            matrix `table_matrix'[`i', 6] = ri_pval_att_sgt[`i']
			matrix `table_matrix'[`i', 7] = weights[`i']
		}
		// Set column names for the matrix
        matrix colnames `table_matrix' = ATT SE pval JKNIFE_SE JKNIFE_pval RI_pval W
        
        // Set row names for the matrix using the gt names
        matrix rownames `table_matrix' = `sgt_names'
        
        // Store the matrix in r()
        return matrix restab = `table_matrix'
		di as text "Aggregation: " as result "sgt"

		
	} 
	
		di as text "Aggregate ATT: " as result agg_att[1]
		di as text "Standard error: " as result agg_att_se[1]
		di as text "p-value: " as result agg_att_pval[1]
		di as text "Jackknife SE: " as result jackknife_se[1]
		di as text "Jackknife p-value: " as result jknife_pval[1]
		di as text "RI p-value: " as result ri_pval_agg_att[1]
		di as text "nperm: " as result nperm[1]
		local linesize = c(linesize)
		if `linesize' < 103 & "`agg'" != "none" {
			di as text "Results table may be squished, try expanding Stata results window."
		}
	
	// Store aggregate results in r()
    return scalar att = agg_att[1]
    return scalar se = agg_att_se[1]
    return scalar p = agg_att_pval[1]
    return scalar jkse = jackknife_se[1]
    return scalar jkp = jknife_pval[1]
    return scalar rip = ri_pval_agg_att[1]
	
	qui frame change default
    qui frame drop `result_frame'
	
end 

/*--------------------------------------*/
/* Change Log */
/*--------------------------------------*/
*0.1.1 - changed column types from Any to Float64 to ensure data is passed to Stata properly and changed argument from save_all_csvs to save_csv and now defaults to true
*0.1.2 - display filepaths for saved .csv files
*0.1.3 - backslashes to forwardslashes fixes Julia-Stata compatability issue
*0.1.4 - added p-values and output displays
*0.2.0 - added weights parameter
*0.3.0 - added computation of jackknife p -value for common treatment with silos >= 3
*0.4.0 - added new se's and p-vals
*0.4.1 - made correction to dof
*0.4.2 - enhanced robustness
*0.5.0 - added RI pvals at sub-aggregate level, changed program to rclass
*0.6.0 - added new weighting arg
*0.6.1 - removed deprecated code for calculating pvals (just doing everything on JL side now)
*0.6.2 - added agg options of sgt and none
*0.6.3 - overwrite blank agg option to agg = "g"
*0.6.4 - changed the way that the results row labels are passed to Stata from Julia to try and work around a Stata-Julia interface bug

