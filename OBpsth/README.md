Sub-directory to get OB PSTH shown in Figure 1 Aii (top row)

plot_OB_psth.m -- script to plot PSTH segmented by modality (in paper), also within a modality segmented by drug preparation (not in paper but extensively in Craft et al. '23 JNeurophys)

plot_OB_pcolor.m -- script to plot trial-avg spike count for a particular rat, showing all simul. recorded neurons. Nonlinear colorbar to see differences. These cells MATCH same recording as plot_aPC_pcolor.m in main directory.

getShow_DecodeAcc_OB.m -- load single neuron decoding accuracy results saved in mat file IndOBdecode_varyWin_EB_0sSpon.mat (created by running indDecode_OB.m). Shows results in Fig 1C for OB.

IndOBdecode_varyWin_EB_0sSpon.mat -- contains all single neuron decoding accuracy when using best possible threshold. Created by script: indDecode_OB.m

Processed mat files of spike data in OB is analogous to main directory, imported from github.com/michellecraft64/OB/Modality/data_analysis/

Rat[j]_indCell _[NoDrug/Bic/Mus].mat where j indicates rat index (arbitrary) and [NoDrug/Bic/Mus] indicates the drug preparation. Each mat file has 4 varibles: setParam (struct), setParam.Twin=0.1 (100ms time windows), setParam.SponSampLength=3 (3sec of spont. befored evoked), setParam.Toverlap=0 (no overlapping windows). TimeVars (struct) has fields numEvok (number of Twin windows in evoked state), FirstSpon (time bin index of first spon time), LastSpon (time bin index of last spont tim), FirstEvok (time bin index of first evoked time), LastEvok (time bin index of last evoked time). The time bin index was shifted by TimeVars.StimShift=15 bins, to make correspondence with original data, but this is not necessary here because [First/Last][Spon/Evok] indices correspond correctly to sOR and sRET variables. sOR (matrix of size 310 x 10 x num_aPC). 310=number of time bins of length 100ms, 10=number of ortho trials, num_aPC=number of aPC neurons (varies). sRET (matrix of size 310 x 10 x num_aPC). 310=number of time bins of length 100ms, 10=number of retro trials, num_aPC=number of aPC neurons (varies). Note: the variables setParam and TimeVars are the SAME for all 15 mat files (redundant).
