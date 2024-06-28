% all first and second order stats using histograms, COMPARING across drugs WITHIN a modality
% WITHIN aPC only
% loops through all stats BUT using the same variable name.
% !! Pauses between each stat, so must hit any key to advance through loop

load DS_stats.mat

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;

for whichStat=1:4 %1=mean, 2=var, 3=cov, 4=correl
    %these names mean_pc_[or/rt]_[ND/Bic/Mus] get overwritten at each loop step
    mean_pc_or_ND=[];
    mean_pc_rt_ND=[];
    mean_pc_or_Bic=[];
    mean_pc_rt_Bic=[];
    mean_pc_or_Mus=[];
    mean_pc_rt_Mus=[];
    switch whichStat
        case 1 %mean spike counts
            for j=1:29
                mean_pc_or_ND=[mean_pc_or_ND pc_or_mean{3}{j}];
                mean_pc_rt_ND=[mean_pc_rt_ND pc_rt_mean{3}{j}];
                if(j<=12) %only 12 recordings for Bic & Mus
                    mean_pc_or_Bic=[mean_pc_or_Bic pc_or_mean{1}{j}];
                    mean_pc_rt_Bic=[mean_pc_rt_Bic pc_rt_mean{1}{j}];
                    mean_pc_or_Mus=[mean_pc_or_Mus pc_or_mean{2}{j}];
                    mean_pc_rt_Mus=[mean_pc_rt_Mus pc_rt_mean{2}{j}];
                end
            end
        case 2
            for j=1:29
                mean_pc_or_ND=[mean_pc_or_ND pc_or_var{3}{j}];
                mean_pc_rt_ND=[mean_pc_rt_ND pc_rt_var{3}{j}];
                if(j<=12) %only 12 recordings for Bic & Mus
                    mean_pc_or_Bic=[mean_pc_or_Bic pc_or_var{1}{j}];
                    mean_pc_rt_Bic=[mean_pc_rt_Bic pc_rt_var{1}{j}];
                    mean_pc_or_Mus=[mean_pc_or_Mus pc_or_var{2}{j}];
                    mean_pc_rt_Mus=[mean_pc_rt_Mus pc_rt_var{2}{j}];
                end
            end
        case 3
            mean_pc_or_ND=pc_or_all{3}(:,1)';
            mean_pc_rt_ND=pc_rt_all{3}(:,1)';
            mean_pc_or_Bic=pc_or_all{1}(:,1)';
            mean_pc_rt_Bic=pc_rt_all{1}(:,1)';
            mean_pc_or_Mus=pc_or_all{2}(:,1)';
            mean_pc_rt_Mus=pc_rt_all{2}(:,1)';
        case 4
            mean_pc_or_ND=pc_or_all{3}(:,2)';
            mean_pc_rt_ND=pc_rt_all{3}(:,2)';
            mean_pc_or_Bic=pc_or_all{1}(:,2)';
            mean_pc_rt_Bic=pc_rt_all{1}(:,2)';
            mean_pc_or_Mus=pc_or_all{2}(:,2)';
            mean_pc_rt_Mus=pc_rt_all{2}(:,2)';
            %get rid of NaN's dividing by 0
            mean_pc_or_ND=mean_pc_or_ND(~isnan(mean_pc_or_ND));
            mean_pc_rt_ND=mean_pc_rt_ND(~isnan(mean_pc_rt_ND));
            mean_pc_or_Bic=mean_pc_or_Bic(~isnan(mean_pc_or_Bic));
            mean_pc_rt_Bic=mean_pc_rt_Bic(~isnan(mean_pc_rt_Bic));
            mean_pc_or_Mus=mean_pc_or_Mus(~isnan(mean_pc_or_Mus));
            mean_pc_rt_Mus=mean_pc_rt_Mus(~isnan(mean_pc_rt_Mus));
    end

    % box plots with whiskers
    figure
    hold on
    xgroupdat=[]; ydat=[]; mn_vl=zeros(6,1);
    for j=1:2
        switch j
            case 1
                tmp=mean_pc_or_ND;
                tmp1=mean_pc_or_Bic;
                tmp2=mean_pc_or_Mus;
            case 2
                tmp=mean_pc_rt_ND;
                tmp1=mean_pc_rt_Bic;
                tmp2=mean_pc_rt_Mus;
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
    switch whichStat
        case 1
            set(gca,'YLim',[0 6]) %mean
        case 2
            set(gca,'YLim',[0 12]) %var
        case 3
            set(gca,'YLim',[-1 2]) %covar
        case 4
            set(gca,'YLim',[-1 1]) %corre
    end

    % --- for Ortho only to start ---
    % p-vals for significance for aPC
    [~,pT_NDbic]=ttest2(mean_pc_or_ND,mean_pc_or_Bic,'VarType','unequal');
    [~,pT_NDmus]=ttest2(mean_pc_or_ND,mean_pc_or_Mus,'VarType','unequal');
    [~,pT_BICmus]=ttest2(mean_pc_or_Bic,mean_pc_or_Mus,'VarType','unequal');
    %WCM rank sum
    [pW_NDbic,~,zW_NDbic]=ranksum(mean_pc_or_ND,mean_pc_or_Bic);
    [pW_NDmus,~,zW_NDmus]=ranksum(mean_pc_or_ND,mean_pc_or_Mus);
    [pW_BICmus,~,zW_BICmus]=ranksum(mean_pc_or_Bic,mean_pc_or_Mus);
    %One-way ANOVA
    g_NDbic=[ones(length(mean_pc_or_ND),1);2*ones(length(mean_pc_or_Bic),1)];
    g_NDmus=[ones(length(mean_pc_or_ND),1);2*ones(length(mean_pc_or_Mus),1)];
    g_BICmus=[ones(length(mean_pc_or_Bic),1);2*ones(length(mean_pc_or_Mus),1)];
    [pA_NDbic,tblA_NDbic] = anova1([mean_pc_or_ND mean_pc_or_Bic],g_NDbic,'off');
    [pA_NDmus,tblA_NDmus] = anova1([mean_pc_or_ND mean_pc_or_Mus],g_NDmus,'off');
    [pA_BICmus,tblA_BICmus] = anova1([mean_pc_or_Bic mean_pc_or_Mus],g_BICmus,'off');

    % table to show ND is stat signif better than Bic/Mus
    drgComp={'ND~=Bic';'ND~=Mus';'Bic~=Mus'};
    Ttest=[pT_NDbic; pT_NDmus; pT_BICmus];
    WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
    owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
    T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

    numND=length(mean_pc_or_ND); numBic=length(mean_pc_or_Bic); numMus=length(mean_pc_or_Mus);

    EffSiz_ttst=[ abs(mean(mean_pc_or_ND)-mean(mean_pc_or_Bic))/sqrt( ((numND-1)*var(mean_pc_or_ND)+(numBic-1)*var(mean_pc_or_Bic))/(numND+numBic-2) ) ; ...
        abs(mean(mean_pc_or_ND)-mean(mean_pc_or_Mus))/sqrt( ((numND-1)*var(mean_pc_or_ND)+(numMus-1)*var(mean_pc_or_Mus))/(numND+numMus-2) ) ; ...
        abs(mean(mean_pc_or_Bic)-mean(mean_pc_or_Mus))/sqrt( ((numBic-1)*var(mean_pc_or_Bic)+(numMus-1)*var(mean_pc_or_Mus))/(numBic+numMus-2) )];
    EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
    EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
    T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

    % --- REPEAT for Retro ---
    % p-vals for significance for aPC ONLY
    [~,pT_NDbic]=ttest2(mean_pc_rt_ND,mean_pc_rt_Bic,'VarType','unequal');
    [~,pT_NDmus]=ttest2(mean_pc_rt_ND,mean_pc_rt_Mus,'VarType','unequal');
    [~,pT_BICmus]=ttest2(mean_pc_rt_Bic,mean_pc_rt_Mus,'VarType','unequal');
    %WCM rank sum
    [pW_NDbic,~,zW_NDbic]=ranksum(mean_pc_rt_ND,mean_pc_rt_Bic);
    [pW_NDmus,~,zW_NDmus]=ranksum(mean_pc_rt_ND,mean_pc_rt_Mus);
    [pW_BICmus,~,zW_BICmus]=ranksum(mean_pc_rt_Bic,mean_pc_rt_Mus);
    %One-way ANOVA
    g_NDbic=[ones(length(mean_pc_rt_ND),1);2*ones(length(mean_pc_rt_Bic),1)];
    g_NDmus=[ones(length(mean_pc_rt_ND),1);2*ones(length(mean_pc_rt_Mus),1)];
    g_BICmus=[ones(length(mean_pc_rt_Bic),1);2*ones(length(mean_pc_rt_Mus),1)];
    [pA_NDbic,tblA_NDbic] = anova1([mean_pc_rt_ND mean_pc_rt_Bic],g_NDbic,'off');
    [pA_NDmus,tblA_NDmus] = anova1([mean_pc_rt_ND mean_pc_rt_Mus],g_NDmus,'off');
    [pA_BICmus,tblA_BICmus] = anova1([mean_pc_rt_Bic mean_pc_rt_Mus],g_BICmus,'off');

    % table to show ND is stat signif better than Bic/Mus
    Ttest=[pT_NDbic; pT_NDmus; pT_BICmus];
    WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
    owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
    T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

    numND=length(mean_pc_rt_ND); numBic=length(mean_pc_rt_Bic); numMus=length(mean_pc_rt_Mus);

    EffSiz_ttst=[ abs(mean(mean_pc_rt_ND)-mean(mean_pc_rt_Bic))/sqrt( ((numND-1)*var(mean_pc_rt_ND)+(numBic-1)*var(mean_pc_rt_Bic))/(numND+numBic-2) ) ; ...
        abs(mean(mean_pc_rt_ND)-mean(mean_pc_rt_Mus))/sqrt( ((numND-1)*var(mean_pc_rt_ND)+(numMus-1)*var(mean_pc_rt_Mus))/(numND+numMus-2) ) ; ...
        abs(mean(mean_pc_rt_Bic)-mean(mean_pc_rt_Mus))/sqrt( ((numBic-1)*var(mean_pc_rt_Bic)+(numMus-1)*var(mean_pc_rt_Mus))/(numBic+numMus-2) )];
    EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
    EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
    T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

    pause
end