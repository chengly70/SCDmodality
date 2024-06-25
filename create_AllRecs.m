%script to get ALL recordings in struct format, save results in All_Recs.mat

clear
% Parameters
odorName='EB';
numTrials=20;
hlfTrials=10;

Twin=0.1; %100ms windows, MUST change if Rat%d.. changes
Tevok=1; %1sec; can change
nmTwins=Tevok/Twin-1;

tme=(0:Twin:Tevok)'; %useTime? must match length, which is =LastEvok-FirstSpon+1

load dSizeCells_perRecord.mat %for OB cells
nmCs_ND_OB=nmCells_ND;
nmCs_Bic_OB=nmCells_Bic;
nmCs_Mus_OB=nmCells_Mus;
clear nmCells_*

load dSize_PC_perRecord.mat %for PC cells
nmCs_ND_PC=nmCells_ND;
nmCs_Bic_PC=nmCells_Bic;
nmCs_Mus_PC=nmCells_Mus;
clear nmCells_*

numRecrd_ND=0; %same # recordings for OB and PC
for j=1:8 %8 total rats for ND
    numRecrd_ND=numRecrd_ND+length(nmCs_ND_OB{j,1});
end
numRecrd_Bic=0; %same # recordings for OB and PC
for j=1:4 %4 total for Bic
    numRecrd_Bic=numRecrd_Bic+length(nmCs_Bic_OB{j,1});
end
numRecrd_Mus=0; %same # recordings for OB and PC
for j=1:3 %3 total for Mus
    numRecrd_Mus=numRecrd_Mus+length(nmCs_Mus_OB{j,1});
end


for drug_to_keep=0:2

switch drug_to_keep
    case 0
        ind=[1 2 6 7 8 9 10 11]; %actual good rats
        drugName='NoDrug';
        nmOBir=nmCs_ND_OB;
        nmPCir=nmCs_ND_PC;
        tmpOBor=cell(numRecrd_ND,1);
        tmpOBrt=cell(numRecrd_ND,1);
        tmpPCor=cell(numRecrd_ND,1);
        tmpPCrt=cell(numRecrd_ND,1);
    case 1
        ind=[1 2 8 11]; %actual good rats
        drugName='Bic';
        nmOBir=nmCs_Bic_OB;
        nmPCir=nmCs_Bic_PC;
        tmpOBor=cell(numRecrd_Bic,1);
        tmpOBrt=cell(numRecrd_Bic,1);
        tmpPCor=cell(numRecrd_Bic,1);
        tmpPCrt=cell(numRecrd_Bic,1);
    case 2
        ind=[6 9 10]; %actual good rats
        drugName='Mus';
        nmOBir=nmCs_Mus_OB;
        nmPCir=nmCs_Mus_PC;
        tmpOBor=cell(numRecrd_Mus,1);
        tmpOBrt=cell(numRecrd_Mus,1);
        tmpPCor=cell(numRecrd_Mus,1);
        tmpPCrt=cell(numRecrd_Mus,1);
end


cnt=1;
cntInRat=1;
for i=ind
    fileNameOB=sprintf('Rat%d_IndCell_%s_%s.mat',i,odorName,drugName);
    %Load file, set parameters
    load(fileNameOB) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
    numEvok=TimeVars.numEvok;
    StimShift=TimeVars.StimShift;
    FirstEvok=TimeVars.FirstEvok;
    LastEvok=FirstEvok+nmTwins; %TimeVars.LastEvok;
    LastSpon=TimeVars.LastSpon;
    FirstSpon=LastSpon-0; %TimeVars.FirstSpon; %HARD CODED 0
    lenEvok=length(FirstEvok:1:LastEvok);
    sRetTmp=sRET(FirstSpon:LastEvok,:,:);
    sOrTmp=sOR(FirstSpon:LastEvok,:,:);
    sRetTmp=squeeze(sum(sRetTmp)); %sum spkCnts in [0,Tevok], 10 x Tot#cells
    sOrTmp=squeeze(sum(sOrTmp)); %sum spkCnts in [0,Tevok], 10 x Tot#cells
    fileNamePC=sprintf('Rat%d_indPC_%s_%s.mat',i,odorName,drugName); %repeat for PC
    load(fileNamePC) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
    sRetTmpPC=sRET(FirstSpon:LastEvok,:,:);
    sOrTmpPC=sOR(FirstSpon:LastEvok,:,:);
    sRetTmpPC=squeeze(sum(sRetTmpPC)); %sum spkCnts in [0,Tevok], 10 x Tot#cells
    sOrTmpPC=squeeze(sum(sOrTmpPC)); %sum spkCnts in [0,Tevok], 10 x Tot#cells

    tmpNob=cumsum(nmOBir{cntInRat,1});
    tmpNpc=cumsum(nmPCir{cntInRat,1});

    tmpOBor{cnt,1}=sOrTmp(:,1:tmpNob(1)); tmpOBrt{cnt,1}=sRetTmp(:,1:tmpNob(1));
    tmpPCor{cnt,1}=sOrTmpPC(:,1:tmpNpc(1)); tmpPCrt{cnt,1}=sRetTmpPC(:,1:tmpNpc(1));
    cnt=cnt+1;
    for j=2:length(tmpNob)
        tmpOBor{cnt,1}=sOrTmp(:,tmpNob(j-1)+1:tmpNob(j));
        tmpOBrt{cnt,1}=sRetTmp(:,tmpNob(j-1)+1:tmpNob(j));
        tmpPCor{cnt,1}=sOrTmpPC(:,tmpNpc(j-1)+1:tmpNpc(j));
        tmpPCrt{cnt,1}=sRetTmpPC(:,tmpNpc(j-1)+1:tmpNpc(j));
        cnt=cnt+1;
    end

    cntInRat=cntInRat+1;
end % loop through each rat

switch drug_to_keep
    case 0
        Ob_or_ND=tmpOBor;
        Ob_rt_ND=tmpOBrt;
        Pc_or_ND=tmpPCor;
        Pc_rt_ND=tmpPCrt;
    case 1
        Ob_or_Bic=tmpOBor;
        Ob_rt_Bic=tmpOBrt;
        Pc_or_Bic=tmpPCor;
        Pc_rt_Bic=tmpPCrt;
    case 2
        Ob_or_Mus=tmpOBor;
        Ob_rt_Mus=tmpOBrt;
        Pc_or_Mus=tmpPCor;
        Pc_rt_Mus=tmpPCrt;
end

end

save All_Recs Ob_or* Ob_rt* Pc_or* Pc_rt*