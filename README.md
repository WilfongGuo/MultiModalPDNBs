# Instruction of MMPDNB package
This package includes Matlab scripts and several datasets for demo of MMPDNB approach:
‘main_MMPDNB.m’ is a Matlab function for the routine of experimental analysis. MMPDNB aims to identify personalized biomarkers (multi-modal PDNBs or PDNBs) for detecting early warning signal of individual patients in cancer and find the drug target genes which can provide effective information for the early treatment of cancer by analyzing the biomarkers.

The input (case: BRCA) include:
(1)Cancer: BRCA.
(2)Path: The path of the user where the 'Main,m' is located.
(3)Popnum: the population size of multi-modal evolution algorithm.
(4)Max_CalNum: the maximum number of function evaluation.
(5)Experiment_num: the number of multi-modal evolution algorithm runs

The output results:
(1)PGIN_BRCA: BRCA patients' personalized gene interaction network construct by SSN method.
i.‘BRCA_i_PGIN.mat’ indicates that the PGIN of the i-th BRCA patient which contains the subnetwork adjacency matrix and the name of the gene in the subnetwork.

(2)BRCA_result: Non-dominated solutions of patient samples obtained by MMPDNB.
i.‘BRCA_sample_i_MMPDNB_boxchart.mat’ stores non-dominated solutions for the i-th patient by running the evolutionary algorithm 30 times each time.
ii.‘BRCA_sample_i_MMPDNB_PF.mat’ indicates that a group Pareto front solutions of i-th patient obtained by performing non-dominated sort on boxchart.
iii.‘BRCA_sample_1_MMPDNB_PS.mat’ indicates the DNB corresponding to Pareto front solutions.

(3)BRCA_DNB_name.txt: Genes’ name of multi-modal PDNB or PDNB with the biggest score  of BRCA patients.
i.For example：i-th patient: patient’s ID. PDNB or multi-modal PDNB: the name of genes in PDNB or multi-modal PDNB.

Suggestions
(1)Hardware suggestions for running this package: Window 10 or above; Matlab 2016 or above; RAM 32G or above.

(2)When users analyzed running this package, please note that:
i.set the path in the program.
ii.Parameter setting of Popnum, Max_CalNum, and Experiment_num will affect the running time. With default parameters, MMPDNB takes about 40 minutes to identify a PDNB in a BRCA patient. Users can decrease running time by modifying above parameter.

%    $Id: Main.m Created at 2021-11-15$ 
%   $Copyright (c) 2021 by School of Electrical Engineering, Zhengzhou University, Zhengzhou 450001, China$; 
%    $If any problem, please contact guowf@zzu.edu.cn for help. $
