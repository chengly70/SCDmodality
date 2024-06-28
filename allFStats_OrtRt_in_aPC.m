% all first and second order stats using histograms, COMPARING ortho/retro within a drug type
% WITHIN aPC only
% loops through all stats BUT using the same variable name.
% !! Pauses between each stat, so must hit any key to advance through loop

load DS_stats.mat

ccD1=[0 0 255; 255 0 0]./255;

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
    for j=1:3
        switch j
            case 1
                tmp=mean_pc_or_ND;
                tmp1=mean_pc_rt_ND;
            case 2
                tmp=mean_pc_or_Bic;
                tmp1=mean_pc_rt_Bic;
            case 3
                tmp=mean_pc_or_Mus;
                tmp1=mean_pc_rt_Mus;
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
   
    % p-vals for significance for aPC
    [~,pT_ND]=ttest2(mean_pc_or_ND,mean_pc_rt_ND,'VarType','unequal');
    [~,pT_Bic]=ttest2(mean_pc_or_Bic,mean_pc_rt_Bic,'VarType','unequal');
    [~,pT_Mus]=ttest2(mean_pc_or_Mus,mean_pc_rt_Mus,'VarType','unequal');
    %WCM rank sum
    [pW_ND,~,zW_ND]=ranksum(mean_pc_or_ND,mean_pc_rt_ND);
    [pW_Bic,~,zW_Bic]=ranksum(mean_pc_or_Bic,mean_pc_rt_Bic);
    [pW_Mus,~,zW_Mus]=ranksum(mean_pc_or_Mus,mean_pc_rt_Mus);
    %One-way ANOVA
    g_ND=[ones(length(mean_pc_or_ND),1);2*ones(length(mean_pc_rt_ND),1)];
    g_Bic=[ones(length(mean_pc_or_Bic),1);2*ones(length(mean_pc_rt_Bic),1)];
    g_Mus=[ones(length(mean_pc_or_Mus),1);2*ones(length(mean_pc_rt_Mus),1)];
    [pA_ND,tblA_ND] = anova1([mean_pc_or_ND mean_pc_rt_ND],g_ND,'off');
    [pA_Bic,tblA_Bic] = anova1([mean_pc_or_Bic mean_pc_rt_Bic],g_Bic,'off');
    [pA_Mus,tblA_Mus] = anova1([mean_pc_or_Mus mean_pc_rt_Mus],g_Mus,'off');

    % table to show ND is stat signif better than Bic/Mus
    drgComp={'ND';'Bic';'Mus'};
    Ttest=[pT_ND; pT_Bic; pT_Mus];
    WRankSum=[pW_ND; pW_Bic; pW_Mus];
    owanova=[pA_ND; pA_Bic; pA_Mus];
    T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

    numND=length(mean_pc_or_ND); numBic=length(mean_pc_or_Bic); numMus=length(mean_pc_or_Mus);

    EffSiz_ttst=[ abs(mean(mean_pc_or_ND)-mean(mean_pc_rt_ND))/sqrt( ((numND-1)*var(mean_pc_or_ND)+(numND-1)*var(mean_pc_rt_ND))/(numND+numND-2) ) ; ...
        abs(mean(mean_pc_or_Bic)-mean(mean_pc_rt_Bic))/sqrt( ((numBic-1)*var(mean_pc_or_Bic)+(numBic-1)*var(mean_pc_rt_Bic))/(numBic+numBic-2) ) ; ...
        abs(mean(mean_pc_or_Mus)-mean(mean_pc_rt_Mus))/sqrt( ((numMus-1)*var(mean_pc_or_Mus)+(numMus-1)*var(mean_pc_rt_Mus))/(numMus+numMus-2) )];
    EffSiz_wrst=abs([zW_ND.zval/sqrt(numND+numND); zW_Bic.zval/sqrt(numBic+numBic) ; zW_Mus.zval/sqrt(numMus+numMus)]);
    EffSiz_owanova=[tblA_ND{6}/tblA_ND{8}; tblA_Bic{6}/tblA_Bic{8} ; tblA_Mus{6}/tblA_Mus{8}];
    T_EffSiz_ORpc=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

    pause
end