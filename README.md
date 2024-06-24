# SCDmodality
Code for results in "Coding odor modality in piriform cortex efficiently with low-dimensional subspaces: a Shared Covariance Decoding approach" 
By: Delaney Selb, Andrea Barreiro, Shree Hari Gautam, Woodrow Shew, Cheng Ly
DOI: xx



Processed mat files: Rat[j]_ indPC_[NoDrug/Bic/Mus].mat where j indicates rat index (arbitrary) and [NoDrug/Bic/Mus] indicates the drug preparation. 
Each mat file has 4 varibles: setParam (struct), setParam.Twin=0.1 (100ms time windows), setParam.SponSampLength=3 (3sec of spont. befored evoked), setParam.Toverlap=0 (no overlapping windows). 
TimeVars (struct) has fields numEvok (number of Twin windows in evoked state), FirstSpon (time bin index of first spon time), LastSpon (time bin index of last spont tim), FirstEvok (time bin index of first evoked time), LastEvok (time bin index of last evoked time). 
The time bin index was shifted by TimeVars.StimShift=15 bins, to make correspondence with original data, but this is not necessary here because [First/Last][Spon/Evok] indices correspond correctly to sOR and sRET variables. 
sOR (matrix of size 310 x 10 x num_aPC). 310=number of time bins of length 100ms, 10=number of ortho trials, num_aPC=number of aPC neurons (varies). 
sRET (matrix of size 310 x 10 x num_aPC). 310=number of time bins of length 100ms, 10=number of retro trials, num_aPC=number of aPC neurons (varies). 
Note: the variables setParam and TimeVars are the SAME for all 15 mat files (redundant). 
