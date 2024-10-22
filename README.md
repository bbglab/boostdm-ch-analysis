# BoostDM-CH manuscript analyses

## Content

This repo contains the source code to reproduce the analysis and figures of the paper:

> **Identification of Clonal Hematopoiesis Driver Mutations through In Silico Saturation Mutagenesis**
  Santiago Demajo, Joan Enric Ramis-Zaldivar, Ferran Muinos, Miguel L Grau, Maria Andrianova, Nuria Lopez-Bigas,
  Abel Gonzalez-Perez; DOI: https://doi.org/10.1158/2159-8290.CD-23-1416

## Repository organisation:  
```
- BoostDM_ch_analyses
	| 1-Observed_mutations 
	| 2-Blueprints  
	| 3-Discovery_index
	| 4-Benchmark
	|   | - Other_cohort
   	| 5-Experimental_assays
	|   | - DNMT3A
	|   | - TP53
    	| 6-BoostDM_cancer_comparison
	|   | - scripts
    	| 7-Expert_curated_rules_comparison
    	| 8-UKB_analyses
	|   | 0_Clinical_phenotypes
	|   | 1_Post_processing_calling
	|   |   | - polN
	|   | 2_Mutation_overview
	|   | 3_Clinical_associations
	|   | 4_BoostDM_scores_tiers
	|   | 5_Rules_comparison
	|   | 6-Fitness_analysis
	|   | 7-New_exclusive_boostDM
-UKB_variant_calling
	| - crams
	| - part1
	| - part2
	| - refs
	| - U2AF1
- Paper_data
	| - BoostDM-cancer
	| - BoostDM-CH
	|   | - discovery	
	|   | - evaluation/CH
    |   | - model_selection
	|   | - prediction
	|   | - splitcv_meta/CH
	|   | - training_meta/CH
	|   | - Results_crossvalidation_50_iterations
	| - Experimental_data
	|   | - DNMT3A
	|   | - TP53
	| - Expert_curated_rules
	| - Intogen-Cancer
	| - Intogen-CH
	| - pfam
```

## Folders descriptions and allocation of jupyter-notebooks and figures:  
- ```BoostDM_ch_analyses```, Contains all code to reproduce all analyses and figures 
    - ```1-Observed_mutations```, Observed mutations from the three training cohorts
    	- 	**Observed_mutations.ipynb**
     		- 	*Figure 1B*
      		- 	*Supplementary Figure S3*
      		- 	*Supplementary Figure S8B*
    - ```2-Blueprints```,
    	-	**Blueprint_models.ipynb**
     		-	*Figure 2A-B*
      		-	*Figure 2C*
      		-	*Figure 1C*
    - ```3-Discovery_index```,
    	-	**Discovery_index.ipynb**
     		-	*Supplementary Figure S1*
    - ```4-Benchmark```, Cross-validations...
    	-	**Cross-validation.ipynb**
     		-	*Figure 1C*
      		-	*Supplemental Figure S1G*
      	-	**Fscore50_plots.ipynb**
      		-	*Figure 1D*
      		-	*Figure 1B*
      		-	*Supplemental Figure S2A*
      	-	**Precision-recall.ipynb**
      		-	*Figure 1E*
      		-	*Supplemental Figure S2B*
        -	```Other_cohort```, Benchmarking with other cohorts
        	-	**GaoExclusive_benchmark.ipynb**
         		-	*Supplemental Figure S3C*
          	-	**JapaneseBiobank_benchmark.ipynb**
           		-	*Supplemental Figure S3C*
        -	```Results_crossvalidation_50_iterations```,
    - ```5-Experimental_assays```, Benchmarking with experimental assays
        -	```DNMT3A```,
        	-	**DNMT3A_training_cohort.ipynb**
         		-	*Supplemental Figure S3A (left)*
          		-	*Figure 1F*
           	-	**DNMT3A_JapaneseBiobank.ipynb**
            		-	*Supplemental Figure S3A (right)*
        -	```TP53```,
        	-	**TP53_Giacomelli_JapaneseBiobank.ipynb**
        		-	*Supplemental Figure S3B*
        	-	**TP53_Kato_Giacomelli_Training_cohort.ipynb**
        		-	*Supplemental Figure S3B*
        	-	**TP53_Kato_JapaneseBiobank.ipynb**
        		-	*Supplemental Figure S3B*
        	-	**TP53_Kotler_JapaneseBiobank.ipynb**
        		-	*Supplemental Figure S3B*
        	-	**TP53_Kotler_Training_cohort.ipynb**
        		-	*Supplemental Figure S3B*
        	-	**TP53_Ursu_JapaneseBiobank.ipynb**
        		-	*Supplemental Figure S3B*
        	-	**TP53_Ursu_Training_cohort.ipynb**
        		-	*Supplemental Figure S3B*
    - ```6-BoostDM_cancer_comparison```,
        -	**BoostDM-CH_cancer_comparison.ipynb**
      		-	*Figure 3*
        	-	*Supplemental Figure S7*
        - 	```scripts```, Contain dependencies and functions needed for the plots in the folder
    - ```7-Expert_curated_rules_comparison```,
    	-	**Expert_curated_rules.ipynb**
    		-	*Supplemental Figure S7*
      		-	*Supplemental Figure S20C*
    	-	**False_positive.ipynb**, False positive mutations from Bick et. al.
    - ```8-UKB_analyses```,
    	- ```0_Clinical_phenotypes```,
        	-	**cancerdata_manually.ipynb**
        	-	**cancerdata_variables.ipynb**
        	-	**covariates.ipynb**
        	-	**cvddata_variables.ipynb**
			-	**Infectiondata_variables.ipynb**
    	- ```1_Post_processing_calling```,
        	-	**Postprocessing.ipynb**
        	-	**Postprocessing_U2AF1.ipynb**
        	-	**Generate_final_files.ipynb**
        	-	**FINAL_merge.ipynb**
         	-	```polN```, Generate pool of normals and excluding list.
            		-	**poN.ipynb**
            		-	**poN_U2AF1.ipynb**
        - ```2_Mutation_overview```,
        	-	**Mut_overview.ipynb**
        		-	*Figure 4B*
        		-	*Supplemental Figure S8A*
        	-	**MutSig.ipynb**
        		-	*Supplemental Figure S8B*
        - ```3_Clinical_associations```,
        	-	**Mut_overview.ipynb**
        		-	*Figure 4C-D*
        		-	*Figure 5*
        		-	*Supplemental Figure S12*
        		-	*Supplemental Figure S13A-E*
        		-	*Supplemental Figure S14*
        		-	*Supplemental Figure S15*
        		-	*Supplemental Figure S16*
        	-	**Clinical_associations_non_observed.ipynb**
        		-	*Supplemental Figure S10A-B left*
        	-	**Clinical_associations_one_observed.ipynb**
        		-	*Supplemental Figure S10A-B right*
        	-	**Ethnicity_associations.ipynb**
        		-	*Supplemental Figure S10C*
        	-	**Age_association_SNV.ipynb**, Age associaton per SNV CH mutation
        - ```4_BoostDM_scores_tiers```,
        	-	**BoostDM_score_tiers.ipynb**
        		-	*Supplemental Figure S9*
        		-	*Supplemental Figure S13F*
        - ```5_Rules_comparison```,
        	-	**BoostDM_Rules_comparisons.ipynb**
        		-	*Supplemental Figure S18*
        	-	**UKB_mut_BoostDM_Rules.ipynb**
        		-	*Supplemental Figure S17*
        		-	*Supplemental Figure S21A*
        - ```6-Fitness_analysis```,
        	-	**1_variants_distribution.ipynb**, Get limit of detection of the cohort per gene
        	-	**2_ML_fitness_R882H.ipynb**, DNMT3A R882H hotspot fitness calculation
        	-	**3_All_mutations_fitness.ipynb**
        		-	*Supplemental Figure S11*
        - ```7-New_exclusive_boostDM```,
        	-	**BoostDM_Rules_comparisons_overview.ipynb**
        		-	*Supplemental Figure S19A-C*
        		-	*Supplemental Figure S20A*
        		-	*Figure 6A*
        	-	**Methods_analysis_exclusive.ipynb**
        		-	*Supplemental Figure S19D*
        		-	*Supplemental Figure S20E*
        		-	*Figure 6E*
        	-	**Methods_analysis_intersect_rules.ipynb**
        		-	*Supplemental Figure S20B, D, F*
        		-	*Figure 6B-D*
        	-	**BoostDM-CH_ALL_CHEK2_tests_VAF.ipynb**
        		-	*Supplemental Figure S21B*
        	-	**BoostDM-CH_ALL_CHEK2_tests.ipynb**
        		-	*Supplemental Figure S21c-D*



## Complementary content

You can access to boostDM-CH source code and documentation in the [boostDM-CH pipeline repository](https://github.com/bbglab/boostdm-pipeline/tree/1.0.1-ch).<br>
You can explore and download the main outputs of boostDM-CH in the [boostDM-CH website](https://www.intogen.org/ch/boostdm/search).<br>

