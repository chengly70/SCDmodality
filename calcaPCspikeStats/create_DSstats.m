% main script that calcs aPC spike stats and saves results in DS_stats.mat 

load("All_Recs.mat");
% cell array where bic is first cell, mus is second cell, nd is third cell
%% for each ob/pc, calculate n cells for each drug (nObBic, nObMus,nObND, nPcBic, nPcMus, nPcND)
nObBic = findNumCells(Ob_or_Bic);
nObMus = findNumCells(Ob_or_Mus);
nObND = findNumCells(Ob_or_ND);
nPcBic = findNumCells(Pc_or_Bic);
nPcMus = findNumCells(Pc_or_Mus);
nPcND = findNumCells(Pc_or_ND);
%% Build 'ob_or_all' cell array
ob_or_all = cell(3,1);
ob_or_all{1,:} = buildRatMatrix(Ob_or_Bic);
ob_or_all{2,:} = buildRatMatrix(Ob_or_Mus);
ob_or_all{3,:} = buildRatMatrix(Ob_or_ND);
ob_or_mean = findMean(Ob_or_Bic,Ob_or_Mus,Ob_or_ND);
ob_or_var = findVar(Ob_or_Bic,Ob_or_Mus,Ob_or_ND);
%% Build 'ob_rt_all' cell array, find mean,
ob_rt_all = cell(3,1);
ob_rt_all{1,:} = buildRatMatrix(Ob_rt_Bic);
ob_rt_all{2,:} = buildRatMatrix(Ob_rt_Mus);
ob_rt_all{3,:} = buildRatMatrix(Ob_rt_ND);
ob_rt_mean = findMean(Ob_rt_Bic,Ob_rt_Mus,Ob_rt_ND);
ob_rt_var = findVar(Ob_rt_Bic,Ob_rt_Mus,Ob_rt_ND);
%% Build 'ob_rt_all' cell array
pc_or_all = cell(3,1);
pc_or_all{1,:} = buildRatMatrix(Pc_or_Bic);
pc_or_all{2,:} = buildRatMatrix(Pc_or_Mus);
pc_or_all{3,:} = buildRatMatrix(Pc_or_ND);
pc_or_mean = findMean(Pc_or_Bic,Pc_or_Mus,Pc_or_ND);
pc_or_var = findVar(Pc_or_Bic,Pc_or_Mus,Pc_or_ND);
%% Build 'pc_rt_all' cell array
pc_rt_all = cell(3,1);
pc_rt_all{1,:} = buildRatMatrix(Pc_rt_Bic);
pc_rt_all{2,:} = buildRatMatrix(Pc_rt_Mus);
pc_rt_all{3,:} = buildRatMatrix(Pc_rt_ND);
pc_rt_mean = findMean(Pc_rt_Bic,Pc_rt_Mus,Pc_rt_ND);
pc_rt_var = findVar(Pc_rt_Bic,Pc_rt_Mus,Pc_rt_ND);
%% Build 'ob_pc_or_all' cell array
ob_pc_or_all = cell(3,1);
ob_pc_or_all{1,:} = buildRatMatrixByObXPc(Ob_or_Bic,Pc_or_Bic);
ob_pc_or_all{2,:} = buildRatMatrixByObXPc(Ob_or_Mus,Pc_or_Mus);
ob_pc_or_all{3,:} = buildRatMatrixByObXPc(Ob_or_ND,Pc_or_ND);
%% Build 'ob_pc_rt_all' cell array
ob_pc_rt_all = cell(3,1);
ob_pc_rt_all{1,:} = buildRatMatrixByObXPc(Ob_rt_Bic,Pc_rt_Bic);
ob_pc_rt_all{2,:} = buildRatMatrixByObXPc(Ob_rt_Mus,Pc_rt_Mus);
ob_pc_rt_all{3,:} = buildRatMatrixByObXPc(Ob_rt_ND,Pc_rt_ND);

save DS_stats