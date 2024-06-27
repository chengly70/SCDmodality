% run through subset of data to get Dop, Dcc1, D_lda; ONLY focusing on Triplets
% this takes A LONG time to run through all of the pairs, change lines 78 & 100 to calc a subset

numTrials=20; %total # of trials
hlfTrls=10;   %# of trials in each of the 2 stim types
lbls=[-ones(hlfTrls,1); ones(hlfTrls,1)]; % -1=Ortho, 1=Retro

load All_Recs.mat

load DS_stats.mat

% -- outputs --- 
Dop_m=cell(3,1); %collapse across recordings; 1stRow=Bic, 2nd=Mus, 3rd=ND, , BUT redundant b/c diff OB trips
DOBop_m=cell(3,1); %optim for OB trips (1st col), 2ndCol=OBcell#, 3rdCol=recording#, all to identify
Dcc1_m=cell(3,1);    %CC1 for PC trips
Dlda_m=cell(3,1);     %lda for PC trips, BUT redundant b/c diff OB trips
ccaR_m=cell(3,1);     %OB,PC corr after project to CCA (2nd col is stores recording# to track tripsID_)
Corr_OBCX=cell(3,1);
CrrOBC=cell(3,1);  %cross-correl in detail (NOT taking avg of ortho vs retro)

TotRecrds=[12; 12; 29];  %1=Bic, 2=Mus, 3=ND
desirdNmbTrips=4050000; %desired # of triplets, will have less

for DrgState=1:3 %1=Bic, 2=Mus, 3=ND

    tic

cnt=1; %start/reset count
for recNmb=1:TotRecrds(DrgState)

    switch DrgState
        case 1 %Bic
            OB_or=Ob_or_Bic{recNmb};
            OB_rt=Ob_rt_Bic{recNmb};
            PC_or=Pc_or_Bic{recNmb};
            PC_rt=Pc_rt_Bic{recNmb};
        case 2 %Mus
            OB_or=Ob_or_Mus{recNmb};
            OB_rt=Ob_rt_Mus{recNmb};
            PC_or=Pc_or_Mus{recNmb};
            PC_rt=Pc_rt_Mus{recNmb};
        case 3 %ND
            OB_or=Ob_or_ND{recNmb};
            OB_rt=Ob_rt_ND{recNmb};
            PC_or=Pc_or_ND{recNmb};
            PC_rt=Pc_rt_ND{recNmb};
    end

    IDcurr=(ob_pc_or_all{DrgState}(:,6)==recNmb); %get indices for current recNmb

    nmOB=size(OB_or,2); %could also use nOb[Bic/Mus/ND](recNmb) or max(ob_pc_[or/rt]_all{DrgState}(IDcurr,7))
    nmPC=size(PC_or,2); %could also use nPc[Bic/Mus/ND](recNmb) or  max(ob_pc_[or/rt]_all{DrgState}(IDcurr,8))

    % to easily ID cells triplets in OB & PC
    oid=nchoosek(1:nmOB,3); %matrix size nmOB*(nmOB-1)(nmOB-2)/6 X 3, cols are cellIDs
    pid=nchoosek(1:nmPC,3); %matrix size nmPC*(nmPC-1)(nmPC-2)/6 X 3, cols are cellIDs

    crsCorr_or=ob_pc_or_all{DrgState}(IDcurr,2);
    crsCorr_rt=ob_pc_rt_all{DrgState}(IDcurr,2);
    cXY= .5*(ob_pc_or_all{DrgState}(IDcurr,2) + ob_pc_rt_all{DrgState}(IDcurr,2)); %proper way to do noise correl
    obi_cxy=ob_pc_or_all{DrgState}(IDcurr,7); %cell ID in OB, correspond. to cXY
    pci_cxy=ob_pc_or_all{DrgState}(IDcurr,8); %cell ID in PC, correspond. to cXY

    %     % loop through randomly chosen tripss, or ALL tripss
    targNumb_PC=round(sqrt(desirdNmbTrips/(TotRecrds(DrgState)))); %target # of trips in each region to get total=^2
    targNum_OB=targNumb_PC;
    if(nmOB*(nmOB-1)*(nmOB-2)/6 > targNum_OB) %only choose rand subset if larger
        ob_prs=sort(randperm( nmOB*(nmOB-1)*(nmOB-2)/6, targNum_OB )); %# OB trips (random w/out repeating)
    else
        ob_prs=1:(nmOB*(nmOB-1)*(nmOB-2)/6);
    end
    if(nmPC*(nmPC-1)*(nmPC-2)/6 > targNumb_PC) %only choose rand subset if larger
        pc_prs=sort(randperm( nmPC*(nmPC-1)*(nmPC-2)/6, targNumb_PC )); %# PC trips (random w/out repeating)
    else
        pc_prs=1:(nmPC*(nmPC-1)*(nmPC-2)/6);
    end

   for tmpPC=pc_prs %forces different PCs
        pc_id=pid(tmpPC,:); %indices of the 3 PC cells
        cxdat=[PC_or(:,[pc_id(1) pc_id(2) pc_id(3)]); PC_rt(:,[pc_id(1) pc_id(2) pc_id(3)])];
            
        if( rank(cxdat)~=3 ) %nonresponsive or not full rank trips
            continue %skip this trips of PC
        end
        %whitening to deal with identical responses
        cxdat=cxdat+randn(numTrials,1)*(1e-6);
        
        % only calc/save optim OB decoding if (OB,PC) valid
        Dop_PC=bruteDecode(cxdat,lbls);

        % do linear fit, 10-fold cross-validation
        MdLinCv=fitcdiscr(cxdat,lbls,'DiscrimType','Linear','KFold',10);
        categ=kfoldPredict(MdLinCv);
        D_lda=sum(categ==lbls)/numTrials;
        % without cross-valid
        MdLin=fitcdiscr(cxdat,lbls);
        categ=predict(MdLin,cxdat);
        D_lda2=sum(categ==lbls)/numTrials;

        for tmpOB=ob_prs %forces different OBs
            ob_id=oid(tmpOB,:); %indices of the 3 OB cells
            obdat=[OB_or(:,[ob_id(1) ob_id(2) ob_id(3)]); OB_rt(:,[ob_id(1) ob_id(2) ob_id(3)])];

            id_crs=[find((obi_cxy==ob_id(1))&(pci_cxy==pc_id(1)));  find((obi_cxy==ob_id(1))&(pci_cxy==pc_id(2))); ...
                find((obi_cxy==ob_id(1))&(pci_cxy==pc_id(3)));  find((obi_cxy==ob_id(2))&(pci_cxy==pc_id(1))); ...
                find((obi_cxy==ob_id(2))&(pci_cxy==pc_id(2)));  find((obi_cxy==ob_id(2))&(pci_cxy==pc_id(3))); ...
                find((obi_cxy==ob_id(3))&(pci_cxy==pc_id(1)));  find((obi_cxy==ob_id(3))&(pci_cxy==pc_id(2))); ...
                find((obi_cxy==ob_id(3))&(pci_cxy==pc_id(3)))];

            crssCrl=cXY(id_crs); %9 values
            %tripss yield NaN for cross-correlation OR nonresponsive/rank-deficient
            if( sum(isnan(crssCrl)) > 0 || rank(obdat)~=3)    %withOUT pThres criteria
                continue %skip these 2x2 b/c no correl, or OB nonresponsive/rank-deficient
            end
            %whitening to deal with identical responses
            obdat=obdat+randn(numTrials,1)*(1e-6);

            %use W's code from PNAS Nexus
            Dop=bruteDecode(obdat,lbls);
            % our CCA method
            [~,~,R,~,V,~]=canoncorr(obdat,cxdat);
            Dcc1=opti1DDecode(V(:,1)',lbls);
            %save all Decodes and cross-correl
            Dop_m{DrgState,1}(cnt,1)=Dop_PC;  %redund(ineffic) BUT easier to track
            Dcc1_m{DrgState,1}(cnt,1)=Dcc1;
            Dlda_m{DrgState,1}(cnt,1:2)=[D_lda D_lda2]; %redund(ineffic) BUT easier to track, 1st col w/ 10-fold crossVal, 2nd col w/OUT
            ccaR_m{DrgState,1}(cnt,1:2)=[R(1) recNmb]; %save CCA max correl AND recording number (2nd col)
            Corr_OBCX{DrgState,1}(cnt,1:9)=[reshape(crssCrl,1,9)]; %save all 9 cXY 
            DOBop_m{DrgState,1}(cnt,1:7)=[Dop ob_id pc_id]; %optim OB, (2nd-7th col to find 3x3tripss)
            
            CrrOBC{DrgState,1}(cnt,1:18)=[reshape(crsCorr_or(id_crs),1,9) reshape(crsCorr_rt(id_crs),1,9)];

            cnt=cnt+1; %increment cnt; b/c if isnan(crssCrl) then dont save
        end
    end
    
    save([pwd,'/Reslt_Trips'],'Dop_m','Dcc1_m','Dlda_m','ccaR_m','Corr_OBCX','DOBop_m','CrrOBC')
end

toc
end
