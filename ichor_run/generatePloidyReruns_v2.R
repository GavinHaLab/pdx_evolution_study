## script to check for ploidy mismatch between pairs of samples and relaunches 

#Rscript logr_ichor.R -i /fh/fast/ha_g/projects/PDX/matrices_full_set_1mb_severe/PDX_armMatrix_bcm_snp_breast.rds --ploidy "c(2,3)" --cores 1 -o results/ &> logr_ichor_logs/bcm_snp_breast.txt


library(data.table)
library(stringr)


args <- commandArgs(TRUE)

paramsFile <- "ploidy_tfx/pairings_tfx_ploidy_comparisons.txt"
matrixDir <- "/fh/fast/ha_g/projects/PDX/matrices_full_set_1mb_severe/"
outScript <- "./reRun_ichorCNA_cmds.sh"
outRerunDir <- "results/"
outSampleLists <- "./reRun_scripts/"
Rscript <- "logr_ichor.R"
logDir <- "logr_ichor_logs_rerun"
dir.create(outSampleLists)
dir.create(logDir)
ploidy_diff_threshold <- 0.2
numCores <- 8

params <- fread(paramsFile)
cohorts <- params[, unique(cohort)]
cohortsWithSpecialIDs <- c("jax_snp_bladder", "jax_snp_colorectal", "jax_snp_breast", "jax_snp_gbm", "jax_snp_luad", "jax_snp_lusc", "jax_snp_melanoma", "jax_snp_otherlung", "jax_snp_othertumor", "jax_snp_ovarian", "jax_snp_sarcoma")
specialPtIDs <- params[cohort %in% cohortsWithSpecialIDs, str_match(sample1, "([0-9A-Z]+)[0-9A-Z_-]+(P[0-9T]+)")[, 2:3]]
params[cohort %in% cohortsWithSpecialIDs, c("patient", "sample1_type") := as.list(data.frame(specialPtIDs))]
params[cohort %in% cohortsWithSpecialIDs, c("sample2_type") := str_match(sample2, "([0-9A-Z]+)[0-9A-Z_-]+(P[0-9T]+)")[, 3]]

otherPtIDs <- params[!cohort %in% cohortsWithSpecialIDs, str_match(sample1, "([0-9A-Z_]+)[0-9A-Z_-]+(P[0-9T]+)")[, 2:3]]
params[!cohort %in% cohortsWithSpecialIDs, c("patient", "sample1_type") := as.list(data.frame(otherPtIDs))]
params[!cohort %in% cohortsWithSpecialIDs, c("sample2_type") := str_match(sample2, "([0-9A-Z_]+)[0-9A-Z_-]+(P[0-9T]+)")[, 3]]

params[sample1_type != "PT", sample1_type := "PDX"]
params[sample2_type != "PT", sample2_type := "PDX"]
numCohorts <- length(cohorts)

batch_cmds <- NULL
samples_2_rerun_allCohorts <- NULL
# for each cohort
for (i in 1:numCohorts){
	params_cohort <- params[cohort == cohorts[i]]
	params_cohort_diff <- params_cohort[ploidy_diff > ploidy_diff_threshold] # & tfx_diff != 0]	
	numPairs <- nrow(params_cohort_diff)
	if (numPairs == 0){
		next
	}
	patients <- unique(params_cohort_diff$patient)
	
	samples_2_rerun <- data.table(sample_to_rerun = character(), ploidy_to_use = numeric(), patient = character())
	# for each pairing
	for (j in 1:length(patients)){

		# flatten by sample within patient 
		# look from original params_cohort, since we want all pairings for patient
		params_cohort_diff_patient <- params_cohort[patient == patients[j]]# | grepl(patients[j], sample1) | grepl(patients[j], sample2)]
		patient_params <- unique(rbind(params_cohort_diff_patient[, .(sample1, sample1_tfx, sample1_ploidy, sample1_type)], params_cohort_diff_patient[, .(sample2, sample2_tfx, sample2_ploidy, sample2_type)], use.names = FALSE))
		
		sample.rerun <- NULL
		ploidy.rerun <- NULL
		## only 1 pair (2 total samples for patient)
		if (nrow(params_cohort_diff_patient) <= 1){
			if (params_cohort_diff_patient$sample1_type == "PT"){ # sample 1 is pt
				sample.rerun <- params_cohort_diff_patient[, sample1]
				ploidy.rerun <- params_cohort_diff_patient[, sample2_ploidy]
			}else if (params_cohort_diff_patient$sample2_type == "PT"){ #sample 2 is pt
				sample.rerun <- params_cohort_diff_patient[, sample2]
				ploidy.rerun <- params_cohort_diff_patient[, sample1_ploidy]
			}else{ # none are pt, but are pdx samples
				## find the sample that has the higher TFx
				if (params_cohort_diff_patient$tfx_diff != 0){
					ind_max <- params_cohort_diff_patient[, apply(cbind(sample1_tfx, sample2_tfx), 1, which.max)]	
				}else{ # both PDX have same TFx estimate
					# find the sample that has the lower ploidy
					ind_max <- params_cohort_diff_patient[, apply(cbind(sample1_ploidy, sample2_ploidy), 1, which.min)]
				}
				ind_min <- abs(ind_max - 2) + 1	
				colName.sample <- paste0("sample", ind_min)
				colName.ploidy <- paste0("sample", ind_max, "_ploidy")
				sample.rerun <- params_cohort_diff_patient[, get(colName.sample)]
				ploidy.rerun <- params_cohort_diff_patient[, get(colName.ploidy)]
			}
		}else{ ## > 1 pair (> 2 total samples for patient)
			## find the PDX sample with highest TFx and there are no ties for max TFx across all samples
			# if (patient_params[sample1_type == "PDX", ][sample1_tfx == max(sample1_tfx), .N] == 1){
			# 	ploidy.rerun <- patient_params[sample1_type == "PDX"][which.max(sample1_tfx), sample1_ploidy]
			# 	sample.rerun <- patient_params[round(sample1_ploidy) != round(ploidy.rerun), sample1]
			# }else{ # all PDX have same TFx estimate
				## find the most common ploidy among ALL samples (PDX and PT)
				ploidy.rerun <- as.numeric(names(sort(table(patient_params[, round(sample1_ploidy)]), decreasing = TRUE))[1])
				sample.rerun <- patient_params[round(sample1_ploidy) != ploidy.rerun, sample1]
			# }	 
		}
		# store the samples and ploidy that need to be re-run
		if (length(sample.rerun) > 0){
			samples_2_rerun <- rbind(samples_2_rerun, data.table(sample.rerun, ploidy.rerun, patients[j]), use.names = FALSE)
		}
	}
	

	samples_2_rerun[, ploidy_to_use := round(ploidy_to_use)]
	samples_2_rerun <- unique(samples_2_rerun)
	#samples_2_rerun <- samples_2_rerun[, floor(median(ploidy_to_use)), by=sample_to_rerun]
	names(samples_2_rerun)[2] <- "ploidy_to_use"
	samples_2_rerun_allCohorts <- rbind(samples_2_rerun_allCohorts, samples_2_rerun)

	#message(cohorts[i], " has ", samples_2_rerun[duplicated(sample_to_rerun), .N], " duplicated sample to rerun.")
	message(cohorts[i], " has ", numPairs, " / ", nrow(params_cohort) ," pairs with discrepant ploidy. ", nrow(samples_2_rerun), " samples to rerun.")

	# output file of samples/ploidy to rerun for this cohort
	outCohortFile <- paste0(outSampleLists, "/", cohorts[i], "_sampleRerunList.txt")
	fwrite(samples_2_rerun, file = outCohortFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

	# add command to batch script
	cohort_matrix_file <- paste0(matrixDir, "PDX_armMatrix_", cohorts[i],".rds")
	log_file <- paste0(logDir, "/", cohorts[i], ".log")
	cmd_str <- paste0("Rscript ", Rscript, " -i ", cohort_matrix_file, " --sampleFile ", outCohortFile, " --ploidy \"c()\" --cores ", numCores, " -o ", outRerunDir, " &> ", log_file, " 2> ", log_file) 

	batch_cmds <- c(batch_cmds, cmd_str)

}


write.table(batch_cmds, file = outScript, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")





