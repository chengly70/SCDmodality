%plots to show LDA vs SCD with pair-by-pair comparison (gray lines)

load Reslt_Pairs_All.mat

%color convention
ccD=[0 200 0; 140 0 255; 128.01 128.01 128.01]./255;

% apply filter, don't use all pairs
DaLda=cell(3,1); DaCc1=cell(3,1);

for drgC=1:3

    DaCc1{drgC}=Dcc1_m{drgC}; %no filter
    DaLda{drgC}=Dlda_m{drgC}(:,1); %no filter

    idRev=DaLda{drgC}<0.5; %find ones that are mis-classified
    DaLda{drgC}(idRev)=1-DaLda{drgC}(idRev); %map to 1-LDA
end

ubNm=5000; %upperbound of numbers shown, naturally

for drgC=1:3
    nmPrs=length(DaCc1{drgC});
    rndInd=sort(randperm(nmPrs,ubNm));
    idUsd=[];
    figure('Renderer', 'painters')
    for j=rndInd
        semilogy([1 2],[DaLda{drgC}(j,1) DaCc1{drgC}(j)],'.-','color',.8*ones(1,3),'MarkerSize',12)
        idUsd=[idUsd; j];
        if(j==idUsd(1))
            hold on
        end

    end

    semilogy([1 2],[mean(DaLda{drgC}) mean(DaCc1{drgC})],'.-','Color',ccD(drgC,:),'MarkerSize',18)
    semilogy([1 2],[median(DaLda{drgC}) median(DaCc1{drgC})],'.--','Color',ccD(drgC,:),'MarkerSize',18)

    axis([.9 2.1 0.5-eps 1+eps])


end

%% do statistical tests to verify p-vals

p_ttst=zeros(3,1); p_wrs=zeros(3,1);  p_oan=zeros(3,1); z_wrs=cell(3,1); tblan=cell(3,1);
for drgC=1:3

    %Two-sample T-test assuming unequal var
    [~,ptmp]=ttest2(DaLda{drgC},DaCc1{drgC},'VarType','unequal');
    p_ttst(drgC,1)=ptmp;

    %Wilcoxon rank sum
    [ptmp,~,ztmp]=ranksum(DaLda{drgC},DaCc1{drgC}); %check out help ranksum for other options just FYI
    p_wrs(drgC,1)=ptmp;
    z_wrs{drgC,1}=ztmp;

    %One-way ANOVA
    lbls=[ones(length(DaLda{drgC}),1);2*ones(length(DaCc1{drgC}),1)]; %create labels for diff size x/yPop, Group elements
    [ptmp,tblmp] = anova1([DaLda{drgC}; DaCc1{drgC}],lbls,'off');
    p_oan(drgC,1)=ptmp;
    tblan{drgC,1}=tblmp;
end
% do effect sizes
% table to show if CCA is stat signif better than LDA
drgComp={'ND';'Bic';'Mus'};
Ttest=p_ttst([3 1 2]);
WRankSum=p_wrs([3 1 2]);
owanova=p_oan([3 1 2]);
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(DaLda{3}); numBic=length(DaLda{1}); numMus=length(DaLda{2});

EffSiz_ttst=[ abs(mean(DaLda{3})-mean(DaCc1{3}))/sqrt( ((numND-1)*var(DaLda{3})+(numND-1)*var(DaCc1{3}))/(numND+numND-2) ) ; ...
    abs(mean(DaLda{1})-mean(DaCc1{1}))/sqrt( ((numBic-1)*var(DaLda{1})+(numBic-1)*var(DaCc1{1}))/(numBic+numBic-2) ) ; ...
    abs(mean(DaLda{2})-mean(DaCc1{2}))/sqrt( ((numMus-1)*var(DaLda{2})+(numMus-1)*var(DaCc1{2}))/(numMus+numMus-2) )];
EffSiz_wrst=abs([z_wrs{3}.zval/sqrt(numND+numND); z_wrs{1}.zval/sqrt(numBic+numBic) ; z_wrs{2}.zval/sqrt(numMus+numMus)]);
EffSiz_owanova=[tblan{3}{6}/tblan{3}{8}; tblan{1}{6}/tblan{1}{8} ; tblan{2}{6}/tblan{2}{8}];
T_EffSiz_ORob=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)
%% histograms
bw=0.05;
for drgC=1:3
    figure
    hold on
    histogram(DaCc1{drgC},'Normalization','probability','BinWidth',bw,'FaceColor',ccD(drgC,:),'EdgeColor',ccD(drgC,:),'BinWidth',.05,'EdgeAlpha',1,'FaceAlpha',1)
    histogram(DaLda{drgC},'Normalization','probability','BinWidth',bw,'FaceColor',ccD(drgC,:),'EdgeColor','k','BinWidth',.05,'EdgeAlpha',0.5,'FaceAlpha',0.5)

end
