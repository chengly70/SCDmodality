% AFTER running getCalc_AllPairs_decode.m, running again to we can just
% append 2 new vars we want, the unique pairs for decoding (LDA, optim)

numTrials=20; %total # of trials
hlfTrls=10;   %# of trials in each of the 2 stim types
lbls=[-ones(hlfTrls,1); ones(hlfTrls,1)]; % -1=Ortho, 1=Retro

load All_Recs.mat

load DS_stats.mat

% -- outputs --- 
Dop_u=cell(3,1); %unique set of PC pairs used, optimal decoding (brute force)
Dlda_u=cell(3,1); %unique set of PC pairs used, LDA

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
        
        % save results
        Dop_u{DrgState,1}(cnt,1)=Dop_PC;
        Dlda_u{DrgState,1}(cnt,1:2)=[D_lda D_lda2];

        cnt=cnt+1; %increment cnt; b/c if isnan(crssCrl) then dont save
   end
   
    
end

save([pwd,'/Reslt_Pairs_All'],'Dop_u','Dlda_u','-append')

toc
end
