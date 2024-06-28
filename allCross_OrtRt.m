% all cross stats (cross-cov & cross-correl) between OB, aPC using box plots, COMPARING across modality WITHIN a drugPrep

load DS_stats.mat

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

ccD1=[255 0 0; 0 0 255]./255;


%% box plots with whiskers

% or vs retr, within a drug prep
figure
hold on
xgroupdat=[]; ydat=[]; mn_vl=zeros(6,1);
for j=1:3
    switch j
        case 1
            tmp=cov_or_ND';
            tmp1=cov_rt_ND';
        case 2
            tmp=cov_or_Bic';
            tmp1=cov_rt_Bic';
        case 3
            tmp=cov_or_Mus';
            tmp1=cov_rt_Mus';
    end
    sz=length(tmp); sz1=length(tmp1); 
    xgroupdat=[xgroupdat; ones(sz,1)*(2*j-1) ; ones(sz1,1)*2*j];
    ydat=[ydat tmp tmp1];
    boxchart(ones(1,sz)*(2*j-1),tmp,'MarkerStyle','none','BoxEdgeColor','b','BoxFaceColor','b')
    boxchart(ones(1,sz1)*2*j,tmp1,'MarkerStyle','none','BoxEdgeColor','r','BoxFaceColor','r')
    tmp=tmp(~isnan(tmp)); tmp1=tmp1(~isnan(tmp1));
    mn_vl(2*j-1:2*j)=[mean(tmp); mean(tmp1)];
    plot([2*j-1 2*j],mn_vl(2*j-1:2*j),'k.-','MarkerSize',22)
end
set(gca,'YLim',[-2 2]) %cov

figure
hold on
xgroupdat=[]; ydat=[]; mn_vl=zeros(6,1);
for j=1:3
    switch j
        case 1
            tmp=crr_or_ND';
            tmp1=crr_rt_ND';
        case 2
            tmp=crr_or_Bic';
            tmp1=crr_rt_Bic';
        case 3
            tmp=crr_or_Mus';
            tmp1=crr_rt_Mus';
    end
    sz=length(tmp); sz1=length(tmp1); 
    xgroupdat=[xgroupdat; ones(sz,1)*(2*j-1) ; ones(sz1,1)*2*j];
    ydat=[ydat tmp tmp1];
    boxchart(ones(1,sz)*(2*j-1),tmp,'MarkerStyle','none','BoxEdgeColor','b','BoxFaceColor','b')
    boxchart(ones(1,sz1)*2*j,tmp1,'MarkerStyle','none','BoxEdgeColor','r','BoxFaceColor','r')
    tmp=tmp(~isnan(tmp)); tmp1=tmp1(~isnan(tmp1));
    mn_vl(2*j-1:2*j)=[mean(tmp); mean(tmp1)];
    plot([2*j-1 2*j],mn_vl(2*j-1:2*j),'k.-','MarkerSize',22)
end
set(gca,'YLim',[-1 1]) %corr


%% p-vals and ES (Effect Sizes)
[~,pT_ND]=ttest2(crr_or_ND,crr_rt_ND,'VarType','unequal');
[~,pT_Bic]=ttest2(crr_or_Bic,crr_rt_Bic,'VarType','unequal');
[~,pT_Mus]=ttest2(crr_or_Mus,crr_rt_Mus,'VarType','unequal');
%WCM rank sum
[pW_ND,~,zW_ND]=ranksum(crr_or_ND,crr_rt_ND);
[pW_Bic,~,zW_Bic]=ranksum(crr_or_Bic,crr_rt_Bic);
[pW_Mus,~,zW_Mus]=ranksum(crr_or_Mus,crr_rt_Mus);
%One-way ANOVA
g_ND=[ones(length(crr_or_ND),1);2*ones(length(crr_rt_ND),1)];
g_Bic=[ones(length(crr_or_Bic),1);2*ones(length(crr_rt_Bic),1)];
g_Mus=[ones(length(crr_or_Mus),1);2*ones(length(crr_rt_Mus),1)];
[pA_ND,tblA_ND] = anova1([crr_or_ND; crr_rt_ND],g_ND,'off');
[pA_Bic,tblA_Bic] = anova1([crr_or_Bic; crr_rt_Bic],g_Bic,'off');
[pA_Mus,tblA_Mus] = anova1([crr_or_Mus; crr_rt_Mus],g_Mus,'off');

% table to show Orth is stat signif diff than Retro
drgComp={'ND';'Bic';'Mus'}; 
Ttest=[pT_ND; pT_Bic; pT_Mus]; 
WRankSum=[pW_ND; pW_Bic; pW_Mus];
owanova=[pA_ND; pA_Bic; pA_Mus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(crr_or_ND); numBic=length(crr_or_Bic); numMus=length(crr_or_Mus);

EffSiz_ttst=[ abs(mean(crr_or_ND)-mean(crr_rt_ND))/sqrt( ((numND-1)*var(crr_or_ND)+(numND-1)*var(crr_rt_ND))/(numND+numND-2) ) ; ...
    abs(mean(crr_or_Bic)-mean(crr_rt_Bic))/sqrt( ((numBic-1)*var(crr_or_Bic)+(numBic-1)*var(crr_rt_Bic))/(numBic+numBic-2) ) ; ...
    abs(mean(crr_or_Mus)-mean(crr_rt_Mus))/sqrt( ((numMus-1)*var(crr_or_Mus)+(numMus-1)*var(crr_rt_Mus))/(numMus+numMus-2) )];
EffSiz_wrst=abs([zW_ND.zval/sqrt(numND+numND); zW_Bic.zval/sqrt(numBic+numBic) ; zW_Mus.zval/sqrt(numMus+numMus)]);
EffSiz_owanova=[tblA_ND{6}/tblA_ND{8}; tblA_Bic{6}/tblA_Bic{8} ; tblA_Mus{6}/tblA_Mus{8}];
T_EffSiz_ORob=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

%-- repeat for Cov
[~,pT_ND]=ttest2(cov_or_ND,cov_rt_ND,'VarType','unequal');
[~,pT_Bic]=ttest2(cov_or_Bic,cov_rt_Bic,'VarType','unequal');
[~,pT_Mus]=ttest2(cov_or_Mus,cov_rt_Mus,'VarType','unequal');
%WCM rank sum
[pW_ND,~,zW_ND]=ranksum(cov_or_ND,cov_rt_ND);
[pW_Bic,~,zW_Bic]=ranksum(cov_or_Bic,cov_rt_Bic);
[pW_Mus,~,zW_Mus]=ranksum(cov_or_Mus,cov_rt_Mus);
%One-way ANOVA
g_ND=[ones(length(cov_or_ND),1);2*ones(length(cov_rt_ND),1)];
g_Bic=[ones(length(cov_or_Bic),1);2*ones(length(cov_rt_Bic),1)];
g_Mus=[ones(length(cov_or_Mus),1);2*ones(length(cov_rt_Mus),1)];
[pA_ND,tblA_ND] = anova1([cov_or_ND; cov_rt_ND],g_ND,'off');
[pA_Bic,tblA_Bic] = anova1([cov_or_Bic; cov_rt_Bic],g_Bic,'off');
[pA_Mus,tblA_Mus] = anova1([cov_or_Mus; cov_rt_Mus],g_Mus,'off');

% table to show Orth is stat signif diff than Retro
drgComp={'ND';'Bic';'Mus'}; 
Ttest=[pT_ND; pT_Bic; pT_Mus]; 
WRankSum=[pW_ND; pW_Bic; pW_Mus];
owanova=[pA_ND; pA_Bic; pA_Mus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(cov_or_ND); numBic=length(cov_or_Bic); numMus=length(cov_or_Mus);

EffSiz_ttst=[ abs(mean(cov_or_ND)-mean(cov_rt_ND))/sqrt( ((numND-1)*var(cov_or_ND)+(numND-1)*var(cov_rt_ND))/(numND+numND-2) ) ; ...
    abs(mean(cov_or_Bic)-mean(cov_rt_Bic))/sqrt( ((numBic-1)*var(cov_or_Bic)+(numBic-1)*var(cov_rt_Bic))/(numBic+numBic-2) ) ; ...
    abs(mean(cov_or_Mus)-mean(cov_rt_Mus))/sqrt( ((numMus-1)*var(cov_or_Mus)+(numMus-1)*var(cov_rt_Mus))/(numMus+numMus-2) )];
EffSiz_wrst=abs([zW_ND.zval/sqrt(numND+numND); zW_Bic.zval/sqrt(numBic+numBic) ; zW_Mus.zval/sqrt(numMus+numMus)]);
EffSiz_owanova=[tblA_ND{6}/tblA_ND{8}; tblA_Bic{6}/tblA_Bic{8} ; tblA_Mus{6}/tblA_Mus{8}];
T_EffSiz_ORob=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)