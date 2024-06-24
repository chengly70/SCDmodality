%show indiv decode in OB ; save_flag=1, =>create Decoding Accuracy mat file: DecodAcc_OB.mat

load IndOBdecode_varyWin_EB_0sSpon.mat

save_flag=0;

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

indWind=10; %index from 1,..,20; using 9=900ms from JNeurophys BUT many are just as good

DecA_ND=[]; DecA_Bic=[]; DecA_Mus=[];
for j=1:8 %total # rats, for ND
    DecA_ND=[DecA_ND; optThrsh_ND{j,4}(indWind,:).'];
    if j<=3 % bic & mus
        DecA_Bic=[DecA_Bic; optThrsh_Bic{j,4}(indWind,:).'];
        DecA_Mus=[DecA_Mus; optThrsh_Mus{j,4}(indWind,:).'];
    end
    if j==4 % bic
        DecA_Bic=[DecA_Bic; optThrsh_Bic{j,4}(indWind,:).'];
    end
end

if(save_flag)
    save DecodAcc_OB DecA_ND DecA_Bic DecA_Mus
end

figure('Renderer', 'Painters');
hold on
histogram(DecA_ND,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(DecA_Bic,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(DecA_Mus,'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
plot(mean(DecA_ND),.23,'.','color',ccD(1,:),'MarkerSize',22)
plot(mean(DecA_Bic),.22,'.','color',ccD(2,:),'MarkerSize',22)
plot(mean(DecA_Mus),.245,'.','color',ccD(3,:),'MarkerSize',22)
box off
set(gca,'FontSize',18)
set(gca,'XLim',[.5 1])

%% p-vals for significance
[~,pT_NDbic]=ttest2(DecA_ND,DecA_Bic,'VarType','unequal');
[~,pT_NDmus]=ttest2(DecA_ND,DecA_Mus,'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(DecA_ND,DecA_Bic);
[pW_NDmus,~,zW_NDmus]=ranksum(DecA_ND,DecA_Mus);
%One-way ANOVA
g_NDbic=[ones(length(DecA_ND),1);2*ones(length(DecA_Bic),1)];
g_NDmus=[ones(length(DecA_ND),1);2*ones(length(DecA_Mus),1)];
[pA_NDbic,tblA_NDbic] = anova1([DecA_ND; DecA_Bic],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([DecA_ND; DecA_Mus],g_NDmus,'off');

% table to show ND is stat signif better than Bic/Mus
drgComp={'Bic<ND';'Mus<ND'}; 
Ttest=[pT_NDbic; pT_NDmus]; 
WRankSum=[pW_NDbic; pW_NDmus];
owanova=[pA_NDbic; pA_NDmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(DecA_ND); numBic=length(DecA_Bic); numMus=length(DecA_Mus);

EffSiz_ttst=[ abs(mean(DecA_ND)-mean(DecA_Bic))/sqrt( ((numND-1)*var(DecA_ND)+(numBic-1)*var(DecA_Bic))/(numND+numBic-2) ) ; ...
    abs(mean(DecA_ND)-mean(DecA_Mus))/sqrt( ((numND-1)*var(DecA_ND)+(numMus-1)*var(DecA_Mus))/(numND+numMus-2) )];
EffSiz_wrst=[zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus)];
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

% pause
% close all
% end