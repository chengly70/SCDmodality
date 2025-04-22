This sub-directory contains functions to calculate the aPC spike count statistics (mean, var, cov, correl).

The function create_DSstats.m calculates a lot of important stats in aPC (x OB too), loads All_Recs.mat (see main directory and description for how this .mat file is created). 
Saves result in DS_stats.mat

NOTE: to create All_Recs[j].mat for j=1 (100ms), 2 (200ms), .., 9 (900m), use the script create_Allrecs.m in main directory BUT change the Tevok variable to j.
To create DS_stats[j].mat for j=1,..,9, use the script create_DSstats.m here BUT change file loaded (line 3) and file saved (line 51) names, append appropriate j.

There are numerous helper functions with descriptions in the various files.
