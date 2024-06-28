% all cross stats (cross-cov & cross-correl) between OB, aPC using box plots, COMPARING across drugs WITHIN a modality

load DS_stats.mat

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;

cov_or_ND=ob_pc_or_all{3}(:,1);
cov_rt_ND=ob_pc_rt_all{3}(:,1);
cov_or_Bic=ob_pc_or_all{1}(:,1);
cov_rt_Bic=ob_pc_rt_all{1}(:,1);
cov_or_Mus=ob_pc_or_all{2}(:,1);
cov_rt_Mus=ob_pc_rt_all{2}(:,1);

crr_or_ND=ob_pc_or_all{3}(:,2);
crr_or_ND=crr_or_ND(~isnan(crr_or_ND));
crr_rt_ND=ob_pc_rt_all{3}(:,2);
crr_rt_ND=crr_rt_ND(~isnan(crr_rt_ND));
crr_or_Bic=ob_pc_or_all{1}(:,2);
crr_or_Bic=crr_or_Bic(~isnan(crr_or_Bic));
crr_rt_Bic=ob_pc_rt_all{1}(:,2);
crr_rt_Bic=crr_rt_Bic(~isnan(crr_rt_Bic));
crr_or_Mus=ob_pc_or_all{2}(:,2);
crr_or_Mus=crr_or_Mus(~isnan(crr_or_Mus));
crr_rt_Mus=ob_pc_rt_all{2}(:,2);
crr_rt_Mus=crr_rt_Mus(~isnan(crr_rt_Mus));

% box plots with whiskers 
    figure
    hold on
    xgroupdat=[]; ydat=[]; mn_vl=zeros(6,1);
    for j=1:2
        switch j
            case 1 %ORTHO first
                tmp=cov_or_ND';
                tmp1=cov_or_Bic';
                tmp2=cov_or_Mus';
            case 2
                tmp=cov_rt_ND';
                tmp1=cov_rt_Bic';
                tmp2=cov_rt_Mus';
        end
        sz=length(tmp);
        sz1=length(tmp1);
        sz2=length(tmp2);
        xgroupdat=[xgroupdat; ones(sz,1)*(3*j-2) ; ones(sz1,1)*(3*j-1); ones(sz2,1)*3*j];
        ydat=[ydat tmp tmp1 tmp2];
        boxchart(ones(1,sz)*(3*j-2),tmp,'MarkerStyle','none','BoxEdgeColor',ccD(1,:),'BoxFaceColor',ccD(1,:))
        boxchart(ones(1,sz1)*(3*j-1),tmp1,'MarkerStyle','none','BoxEdgeColor',ccD(2,:),'BoxFaceColor',ccD(2,:))
        boxchart(ones(1,sz2)*3*j,tmp2,'MarkerStyle','none','BoxEdgeColor',ccD(3,:),'BoxFaceColor',ccD(3,:))
        tmp=tmp(~isnan(tmp)); tmp1=tmp1(~isnan(tmp1)); tmp2=tmp2(~isnan(tmp2));
        mn_vl(3*j-2:3*j)=[mean(tmp); mean(tmp1) ; mean(tmp2)];
        plot(3*j-2:3*j,mn_vl(3*j-2:3*j),'k.-','MarkerSize',22)
    end
    set(gca,'YLim',[-2 2]) %cov
    figure
    hold on
    xgroupdat=[]; ydat=[]; mn_vl=zeros(6,1);
    for j=1:2
        switch j
            case 1 %ORTHO first
                tmp=crr_or_ND';
                tmp1=crr_or_Bic';
                tmp2=crr_or_Mus';
            case 2
                tmp=crr_rt_ND';
                tmp1=crr_rt_Bic';
                tmp2=crr_rt_Mus';
        end
        sz=length(tmp);
        sz1=length(tmp1);
        sz2=length(tmp2);
        xgroupdat=[xgroupdat; ones(sz,1)*(3*j-2) ; ones(sz1,1)*(3*j-1); ones(sz2,1)*3*j];
        ydat=[ydat tmp tmp1 tmp2];
        boxchart(ones(1,sz)*(3*j-2),tmp,'MarkerStyle','none','BoxEdgeColor',ccD(1,:),'BoxFaceColor',ccD(1,:))
        boxchart(ones(1,sz1)*(3*j-1),tmp1,'MarkerStyle','none','BoxEdgeColor',ccD(2,:),'BoxFaceColor',ccD(2,:))
        boxchart(ones(1,sz2)*3*j,tmp2,'MarkerStyle','none','BoxEdgeColor',ccD(3,:),'BoxFaceColor',ccD(3,:))
        tmp=tmp(~isnan(tmp)); tmp1=tmp1(~isnan(tmp1)); tmp2=tmp2(~isnan(tmp2));
        mn_vl(3*j-2:3*j)=[mean(tmp); mean(tmp1) ; mean(tmp2)];
        plot(3*j-2:3*j,mn_vl(3*j-2:3*j),'k.-','MarkerSize',22)
    end
    set(gca,'YLim',[-1 1]) %corr

% --- for Ortho, Corr ---

% p-vals for significance for crr-retro ONLY
[~,pT_NDbic]=ttest2(crr_or_ND,crr_or_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(crr_or_ND,crr_or_Mus,'VarType','unequal');
[~,pT_BICmus]=ttest2(crr_or_Bic,crr_or_Mus,'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(crr_or_ND,crr_or_Bic);
[pW_NDmus,~,zW_NDmus]=ranksum(crr_or_ND,crr_or_Mus);
[pW_BICmus,~,zW_BICmus]=ranksum(crr_or_Bic,crr_or_Mus);
%One-way ANOVA
g_NDbic=[ones(length(crr_or_ND),1);2*ones(length(crr_or_Bic),1)];
g_NDmus=[ones(length(crr_or_ND),1);2*ones(length(crr_or_Mus),1)];
g_BICmus=[ones(length(crr_or_Bic),1);2*ones(length(crr_or_Mus),1)];
[pA_NDbic,tblA_NDbic] = anova1([crr_or_ND ; crr_or_Bic],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([crr_or_ND ; crr_or_Mus],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([crr_or_Bic ; crr_or_Mus],g_BICmus,'off');

% table to show ND diff
drgComp={'Bic~=ND';'ND~=Mus';'Bic~=Mus'}; 
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(crr_or_ND); numBic=length(crr_or_Bic); numMus=length(crr_or_Mus);

EffSiz_ttst=[ abs(mean(crr_or_ND)-mean(crr_or_Bic))/sqrt( ((numND-1)*var(crr_or_ND)+(numBic-1)*var(crr_or_Bic))/(numND+numBic-2) ) ; ...
    abs(mean(crr_or_ND)-mean(crr_or_Mus))/sqrt( ((numND-1)*var(crr_or_ND)+(numMus-1)*var(crr_or_Mus))/(numND+numMus-2) ); ...
    abs(mean(crr_or_Bic)-mean(crr_or_Mus))/sqrt( ((numBic-1)*var(crr_or_Bic)+(numMus-1)*var(crr_or_Mus))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

% p-vals for significance for cov-retro ONLY
[~,pT_NDbic]=ttest2(cov_or_ND,cov_or_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(cov_or_ND,cov_or_Mus,'VarType','unequal');
[~,pT_BICmus]=ttest2(cov_or_Bic,cov_or_Mus,'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(cov_or_ND,cov_or_Bic);
[pW_NDmus,~,zW_NDmus]=ranksum(cov_or_ND,cov_or_Mus);
[pW_BICmus,~,zW_BICmus]=ranksum(cov_or_Bic,cov_or_Mus);
%One-way ANOVA
g_NDbic=[ones(length(cov_or_ND),1);2*ones(length(cov_or_Bic),1)];
g_NDmus=[ones(length(cov_or_ND),1);2*ones(length(cov_or_Mus),1)];
g_BICmus=[ones(length(cov_or_Bic),1);2*ones(length(cov_or_Mus),1)];
[pA_NDbic,tblA_NDbic] = anova1([cov_or_ND ; cov_or_Bic],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([cov_or_ND;  cov_or_Mus],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([cov_or_Bic ; cov_or_Mus],g_BICmus,'off');

% table to show ND diff
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(cov_or_ND); numBic=length(cov_or_Bic); numMus=length(cov_or_Mus);

EffSiz_ttst=[ abs(mean(cov_or_ND)-mean(cov_or_Bic))/sqrt( ((numND-1)*var(cov_or_ND)+(numBic-1)*var(cov_or_Bic))/(numND+numBic-2) ) ; ...
    abs(mean(cov_or_ND)-mean(cov_or_Mus))/sqrt( ((numND-1)*var(cov_or_ND)+(numMus-1)*var(cov_or_Mus))/(numND+numMus-2) ); ...
     abs(mean(cov_or_Bic)-mean(cov_or_Mus))/sqrt( ((numBic-1)*var(cov_or_Bic)+(numMus-1)*var(cov_or_Mus))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)


%% --- REPEAT for Retro ---
% p-vals for significance for crr-retro ONLY
[~,pT_NDbic]=ttest2(crr_rt_ND,crr_rt_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(crr_rt_ND,crr_rt_Mus,'VarType','unequal');
[~,pT_BICmus]=ttest2(crr_rt_Bic,crr_rt_Mus,'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(crr_rt_ND,crr_rt_Bic);
[pW_NDmus,~,zW_NDmus]=ranksum(crr_rt_ND,crr_rt_Mus);
[pW_BICmus,~,zW_BICmus]=ranksum(crr_rt_Bic,crr_rt_Mus);
%One-way ANOVA
g_NDbic=[ones(length(crr_rt_ND),1);2*ones(length(crr_rt_Bic),1)];
g_NDmus=[ones(length(crr_rt_ND),1);2*ones(length(crr_rt_Mus),1)];
g_BICmus=[ones(length(crr_rt_Bic),1);2*ones(length(crr_rt_Mus),1)];
[pA_NDbic,tblA_NDbic] = anova1([crr_rt_ND ; crr_rt_Bic],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([crr_rt_ND ; crr_rt_Mus],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([crr_rt_Bic ; crr_rt_Mus],g_BICmus,'off');

% table to show ND diff
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(crr_rt_ND); numBic=length(crr_rt_Bic); numMus=length(crr_rt_Mus);

EffSiz_ttst=[ abs(mean(crr_rt_ND)-mean(crr_rt_Bic))/sqrt( ((numND-1)*var(crr_rt_ND)+(numBic-1)*var(crr_rt_Bic))/(numND+numBic-2) ) ; ...
    abs(mean(crr_rt_ND)-mean(crr_rt_Mus))/sqrt( ((numND-1)*var(crr_rt_ND)+(numMus-1)*var(crr_rt_Mus))/(numND+numMus-2) ); ...
    abs(mean(crr_rt_Bic)-mean(crr_rt_Mus))/sqrt( ((numBic-1)*var(crr_rt_Bic)+(numMus-1)*var(crr_rt_Mus))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

% p-vals for significance for cov-retro ONLY
[~,pT_NDbic]=ttest2(cov_rt_ND,cov_rt_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(cov_rt_ND,cov_rt_Mus,'VarType','unequal');
[~,pT_BICmus]=ttest2(cov_rt_Bic,cov_rt_Mus,'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(cov_rt_ND,cov_rt_Bic);
[pW_NDmus,~,zW_NDmus]=ranksum(cov_rt_ND,cov_rt_Mus);
[pW_BICmus,~,zW_BICmus]=ranksum(cov_rt_Bic,cov_rt_Mus);
%One-way ANOVA
g_NDbic=[ones(length(cov_rt_ND),1);2*ones(length(cov_rt_Bic),1)];
g_NDmus=[ones(length(cov_rt_ND),1);2*ones(length(cov_rt_Mus),1)];
g_BICmus=[ones(length(cov_rt_Bic),1);2*ones(length(cov_rt_Mus),1)];
[pA_NDbic,tblA_NDbic] = anova1([cov_rt_ND ; cov_rt_Bic],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([cov_rt_ND;  cov_rt_Mus],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([cov_rt_Bic ; cov_rt_Mus],g_BICmus,'off');

% table to show ND diff
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(cov_rt_ND); numBic=length(cov_rt_Bic); numMus=length(cov_rt_Mus);

EffSiz_ttst=[ abs(mean(cov_rt_ND)-mean(cov_rt_Bic))/sqrt( ((numND-1)*var(cov_rt_ND)+(numBic-1)*var(cov_rt_Bic))/(numND+numBic-2) ) ; ...
    abs(mean(cov_rt_ND)-mean(cov_rt_Mus))/sqrt( ((numND-1)*var(cov_rt_ND)+(numMus-1)*var(cov_rt_Mus))/(numND+numMus-2) ); ...
     abs(mean(cov_rt_Bic)-mean(cov_rt_Mus))/sqrt( ((numBic-1)*var(cov_rt_Bic)+(numMus-1)*var(cov_rt_Mus))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)


