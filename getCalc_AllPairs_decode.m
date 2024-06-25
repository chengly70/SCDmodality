% run through all data to get Dop, Dcc1, D_lda; ONLY focusing on pairs
% this takes A LONG time to run through all of the pairs, change lines 67 & 89 to calc a subset

numTrials=20; %total # of trials
hlfTrls=10;   %# of trials in each of the 2 stim types
lbls=[-ones(hlfTrls,1); ones(hlfTrls,1)]; % -1=Ortho, 1=Retro

load All_Recs.mat

load DS_stats.mat

% -- outputs --- 
Dop_m=cell(3,1); %collapse across recordings; 1stRow=Bic, 2nd=Mus, 3rd=ND, , BUT redundant b/c diff OB pair
DOBop_m=cell(3,1); %optim for OB pair (1st col), 2ndCol=OBcell#, 3rdCol=recording#, all to identify
Dcc1_m=cell(3,1);    %CC1 for PC pair
Dlda_m=cell(3,1);     %lda for PC pair, BUT redundant b/c diff OB pair
ccaR_m=cell(3,1);     %OB,PC corr after project to CCA (2nd col is stores recording# to track PairID_)
Corr_OBCX=cell(3,1); %Cross-correl taking the avg using homoscedasticity
CrrOBC=cell(3,1);  %cross-correl in detail (NOT taking avg of ortho vs retro)

TotRecrds=[12; 12; 29];  %1=Bic, 2=Mus, 3=ND

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

    % to easily ID cells/pairs in OB & PC
    nnM=reshape((1:nmOB*nmOB),nmOB,nmOB);
    all_in=nnM(triu(true([nmOB nmOB]),1)); %all indices for diff pairs
    [fid_ob,sid_ob]=ind2sub([nmOB nmOB],all_in);
    nnM=reshape((1:nmPC*nmPC),nmPC,nmPC);
    all_in=nnM(triu(true([nmPC nmPC]),1)); %all indices for diff pairs
    [fid_pc,sid_pc]=ind2sub([nmPC nmPC],all_in);

    crsCorr_or=ob_pc_or_all{DrgState}(IDcurr,2);
    crsCorr_rt=ob_pc_rt_all{DrgState}(IDcurr,2);
    cXY= .5*(crsCorr_or + crsCorr_rt); %proper way to do noise correl
    obi_cxy=ob_pc_or_all{DrgState}(IDcurr,7); %cell ID in OB, correspond. to cXY
    pci_cxy=ob_pc_or_all{DrgState}(IDcurr,8); %cell ID in PC, correspond. to cXY

   for tmpPC=1:(nmPC*(nmPC-1)/2) %forces different PCs
        pc_id=[fid_pc(tmpPC) ; sid_pc(tmpPC)]; %indices of the 2 PC cells
        cxdat=[PC_or(:,[pc_id(1) pc_id(2)]); PC_rt(:,[pc_id(1) pc_id(2)])];
            
        if( rank(cxdat)~=2 ) %nonresponsive or not full rank pair
            continue %skip this pair of PC
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

        for tmpOB=1:(nmOB*(nmOB-1)/2) %forces different OBs
            ob_id=[fid_ob(tmpOB) ; sid_ob(tmpOB)]; %indices of the 2 OB cells
            obdat=[OB_or(:,[ob_id(1) ob_id(2)]); OB_rt(:,[ob_id(1) ob_id(2)])];

            id_crs=[find((obi_cxy==ob_id(1))&(pci_cxy==pc_id(1)));  find((obi_cxy==ob_id(2))&(pci_cxy==pc_id(1))); ...
                find((obi_cxy==ob_id(1))&(pci_cxy==pc_id(2)));  find((obi_cxy==ob_id(2))&(pci_cxy==pc_id(2)))];

            crssCrl=cXY(id_crs);
            %pairs yield NaN for cross-correlation OR nonresponsive/rank-deficient
            if( sum(isnan(crssCrl)) > 0 || rank(obdat)~=2)    %withOUT pThres criteria
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
            Corr_OBCX{DrgState,1}(cnt,1:4)=[reshape(crssCrl,1,4)]; %save all 4 cXY 
            DOBop_m{DrgState,1}(cnt,1:5)=[Dop ob_id' pc_id']; %optim OB, (2nd-5th col to find 2x2pairs)

            CrrOBC{DrgState,1}(cnt,1:8)=[reshape(crsCorr_or(id_crs),1,4) reshape(crsCorr_rt(id_crs),1,4)];

            cnt=cnt+1; %increment cnt; b/c if isnan(crssCrl) then dont save
        end
    end
    
end
save([pwd,'/Reslt_Pairs_All'],'Dop_m','Dcc1_m','Dlda_m','ccaR_m','Corr_OBCX','DOBop_m','CrrOBC')

toc
end
