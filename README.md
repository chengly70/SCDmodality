# SCDmodality
Code for results in "Coding odor modality in piriform cortex efficiently with low-dimensional subspaces: a Shared Covariance Decoding approach" 
By: Delaney Selb, Andrea Barreiro, Shree Hari Gautam, Woodrow Shew, Cheng Ly
DOI: xx

plot_aPC_psth.m -- script to plot pop- and trial-avg PSTH segmented by modality within a drug prep (in paper, Fig 1 Aii, bottom row), but also segmented by drug prep within a modality (not in paper)

plot_aPC_pcolor.m -- script to plot trial-avg spike count for a particular rat, showing all simul. recorded neurons. Nonlinear colorbar to see differences.

getShow_DecodeAcc_PC.m -- load single neuron decoding accuracy results saved in mat file IndPCdecode_varyWin_EB_0sSpon.mat (created by running indDecode_PC.m). Shows results in Fig 1C for aPC.

cf_OBaPC.m -- script to show p-values and effect sizes for single neuron decoding: that LDA < SCD < Optim by drug prep. Results in Appendix C.

plot_aPC_sparse.m -- script to show how sparse responses are in a trial (avg) as a percentage of neurons spiking to stim, plotted as a function of time by modality and drug prep -- shown in Appendix B.

For all 2x2 (pairs) decoding results, the main funciton is: getCalc_AllPairs_decode.m. This function calculates the LDA, SCD, and Optim linear decoding on all pairs and saves relevant statistics. 
This file takes A LONG time to run since it loops through all possible pairs, so edit the for-loops accordingly for shorter runs (see comment in file). Saves results in Reslt_Pairs_All.mat. 
This function loads 2 mat files: All_Recs.mat, DS_stats.mat, and relies on helper functions: fitcdiscr (MATLAB), kfoldPredict (MATLAB), canoncorr (MATLAB), opt1DDecode.m (finds optimal threshold after data is projected onto hyperplane), bruteDecode.m for optim linear by exhaustively surveying all hyperplans (from our PNAS Nexus paper).

For decoding results on unique pairs (LDA, optimal linear), the main function is getCalc_AllPairs_decodeAppend.m that must be run AFTER getCalc_AllPairs_decode.m (or have viable file Reslt_Pairs_All.mat), this appends the cell structures: Dop_u, Dlda_u to Reslt_Pairs_All.mat. 

plots_cfLDcca.m -- pair-by-pair comparison of LDA with SCD (Fig 2B), also shows stat. signif. SCD better than LDA (not reported in paper because so obvious)

plots_4decodes.m -- shows all decoding accuracy distributions for pairs (and 2x2) using LDA, SCD, and Optimal Linear.  Also shows statistical test results (p-values, effect sizes) highlighted in Results section. Figure 3A. 
End of script now outputs: p-values and effect sizes comparing differences for all 2x2 (between LDA, SCA, Optim), results in Appendix C, ALSO outputs p-values and effect sizes for comparing differences of whole decoding accuracy distributions within an algorithm in the command line.

For Theory plots (E & rho positive correlated, E & R1 negatively correlated), run the script: displayTheory.m . This will show all pairs 2x2 by default, BUT can change the variable flag_pairs=0 to show triplet results. 
Figures 3B, 4B. 

For the 3x3 (triplets) decoding, analogous script is: getCalc_Trips_decode.m. Just like for pairs above, this calculates LDA, SCD, and Optim linear decoding on a very large subset of all 3x3 triplets and saves relevant statistics. 
This files takes A LONG time to run since we consider a large number of 3x3 networks.  Saves results in Reslt_Trips.mat

plots_4decodes_Trips.m -- analogous to plots_4decodes.m BUT for all 3x3 networks. Figure 4A. End of script now outputs: p-values and effect sizes comparing differences for all 2x2 (between LDA, SCA, Optim), results in Appendix C, ALSO outputs p-values and effect sizes for comparing differences of whole decoding accuracy distributions within an algorithm in the command line.

Spike Stats Comparisons.
allFStats_cfDrugs_in_aPC.m -- script to show all box plots (Fig B2) comparing drug effects within a modality, testing all 4 stats (mean,var,cov,correl) WITHIN aPC, and corresponding p-values & effect sizes for all tests (shown in command window).
allFStats_OrtRt_in_aPC.m -- script to show all box plots (Fig B2) comparing modality within a drug prep, testing all 4 stats (mean,var,cov,correl) WITHIN aPC, and corresponding p-values & effect sizes for all tests (shown in command window).
allCross_cfDrugs.m -- script to show box plots (Fig B2) comparing drug effects within a modality, testing cross-cov and cross-correl (between OB & aPC). 
allCross_OrtRt.m -- script to show box plots (Fig B2) comparing modality within a drug prep, testing cross-cov and cross-correl (between OB & aPC).

M-files that create components used in main scripts/functions to calculate and show results in Figures:
create_AllRecs.m -- scripts that creates file All_Recs.mat, aggregates data across rat recordings (loading Rat[j]_ indPC_[NoDrug/Bic/Mus].mat files, see below). 
All_Recs.mat -- has cell variables of summed spike counts in 1s odor evoked period, named [Ob/Pc]_ [or/rt]_[ND/Bic/Mus] that indicate region (Ob or aPC), modality (ortho or retro), and drug preparation (ND,Bic,Mus). Each cell array has size (number recordings)x 1, and each element of the cell array is a matrix of size 10x(number of neurons), where 10 is the number of trials.

Processed mat files: Rat[j]_ indPC_[NoDrug/Bic/Mus].mat where j indicates rat index (arbitrary) and [NoDrug/Bic/Mus] indicates the drug preparation. 
Each mat file has 4 varibles: setParam (struct), setParam.Twin=0.1 (100ms time windows), setParam.SponSampLength=3 (3sec of spont. befored evoked), setParam.Toverlap=0 (no overlapping windows). 
TimeVars (struct) has fields numEvok (number of Twin windows in evoked state), FirstSpon (time bin index of first spon time), LastSpon (time bin index of last spont tim), FirstEvok (time bin index of first evoked time), LastEvok (time bin index of last evoked time). 
The time bin index was shifted by TimeVars.StimShift=15 bins, to make correspondence with original data, but this is not necessary here because [First/Last][Spon/Evok] indices correspond correctly to sOR and sRET variables. 
sOR (matrix of size 310 x 10 x num_aPC). 310=number of time bins of length 100ms, 10=number of ortho trials, num_aPC=number of aPC neurons (varies). 
sRET (matrix of size 310 x 10 x num_aPC). 310=number of time bins of length 100ms, 10=number of retro trials, num_aPC=number of aPC neurons (varies). 
Note: the variables setParam and TimeVars are the SAME for all 15 mat files (redundant). 

IndPCdecode_varyWin_EB_0sSpon.mat -- contains all single neuron decoding accuracy when using best possible threshold. Created by script: indDecode_PC.m

dSize_PC_perRecord.mat -- contains number of aPC cells in given recording (same data in Rat[j]_indPC .mat files)
dSizeCells_perRecord.mat -- contains number of OB cells in a given recording (same data as in Rat[j]_indCell .mat files)
