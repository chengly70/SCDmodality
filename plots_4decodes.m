%script to plot various decoding accurac (uncoupled from the theory, just pure decoding over ALL pairs in aPC and/or OB
% by 'all' we mean responsive, that is not rank deficient, for pairs! for Pairs

%%
load('Reslt_Pairs_all.mat','Dlda_u','Dop_u','Dcc1_m')

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

idRev=Dlda_u{1}(:,1)<0.5; %find ones that are mis-classified
Dlda_u{1}(idRev)=1-Dlda_u{1}(idRev); %map to 1-LDA
idRev=Dlda_u{2}(:,1)<0.5; %find ones that are mis-classified
Dlda_u{2}(idRev)=1-Dlda_u{2}(idRev); %map to 1-LDA
idRev=Dlda_u{3}(:,1)<0.5; %find ones that are mis-classified
Dlda_u{3}(idRev)=1-Dlda_u{3}(idRev); %map to 1-LDA

% LDA applied to 'all' pairs of aPC
figure('Renderer', 'Painters');
hold on
histogram(Dlda_u{3}(:,1),'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(Dlda_u{1}(:,1),'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(Dlda_u{2}(:,1),'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
yHgh=0.2;
plot(mean(Dlda_u{3}(:,1)),yHgh,'.','color',ccD(1,:))
plot(mean(Dlda_u{1}(:,1)),yHgh,'.','color',ccD(2,:))
plot(mean(Dlda_u{2}(:,1)),yHgh,'.','color',ccD(3,:))
set(gca,'XLim',[0.5 1])

% 'all' SCD algorithm without any filtering 
figure('Renderer', 'Painters');
hold on
histogram(Dcc1_m{3},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(Dcc1_m{1},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(Dcc1_m{2},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
yHgh=0.3;
plot(mean(Dcc1_m{3}(:,1)),yHgh,'.','color',ccD(1,:))
plot(mean(Dcc1_m{1}(:,1)),yHgh,'.','color',ccD(2,:))
plot(mean(Dcc1_m{2}(:,1)),yHgh,'.','color',ccD(3,:))
set(gca,'XLim',[0.5 1])

% BruteForce optimal applied to 'all' pairs of aPC
figure('Renderer', 'Painters');
hold on
histogram(Dop_u{3},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
histogram(Dop_u{1},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
histogram(Dop_u{2},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor','none','EdgeColor',ccD(3,:),'LineWidth',2)
yHgh=0.3;
plot(mean(Dop_u{3}(:,1)),yHgh,'.','color',ccD(1,:))
plot(mean(Dop_u{1}(:,1)),yHgh,'.','color',ccD(2,:))
plot(mean(Dop_u{2}(:,1)),yHgh,'.','color',ccD(3,:))
set(gca,'XLim',[0.5 1])

% LDA applied to 'all' pairs of aPC
yHgh=0.2;
figure('Renderer', 'Painters');
hold on
histogram(Dlda_u{3}(:,1),'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
plot(mean(Dlda_u{3}(:,1)),yHgh,'.','color',ccD(1,:),'MarkerSize',12)
plot(median(Dlda_u{3}(:,1)),yHgh+.05,'o','color',ccD(1,:),'MarkerSize',8)
axis([0.5 1 0 0.3])
figure('Renderer', 'Painters');
hold on
histogram(Dlda_u{1}(:,1),'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
plot(mean(Dlda_u{1}(:,1)),yHgh,'.','color',ccD(2,:),'MarkerSize',12)
plot(median(Dlda_u{1}(:,1)),yHgh+.05,'o','color',ccD(2,:),'MarkerSize',8)
axis([0.5 1 0 0.3])
figure('Renderer', 'Painters');
hold on
histogram(Dlda_u{2}(:,1),'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(3,:),'LineStyle','none')
plot(mean(Dlda_u{2}(:,1)),yHgh,'.','color',ccD(3,:),'MarkerSize',12)
plot(median(Dlda_u{2}(:,1)),yHgh+.05,'o','color',ccD(3,:),'MarkerSize',8)
axis([0.5 1 0 0.3])

yHgh=0.28;
figure('Renderer', 'Painters');
hold on
histogram(Dcc1_m{3},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
plot(mean(Dcc1_m{3}),yHgh,'.','color',ccD(1,:),'MarkerSize',12)
plot(median(Dcc1_m{3}),yHgh+.02,'o','color',ccD(1,:),'MarkerSize',8)
axis([0.5 1 0 0.3])
figure('Renderer', 'Painters');
hold on
histogram(Dcc1_m{1},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
plot(mean(Dcc1_m{1}),yHgh,'.','color',ccD(2,:),'MarkerSize',12)
plot(median(Dcc1_m{1}),yHgh+.02,'o','color',ccD(2,:),'MarkerSize',8)
axis([0.5 1 0 0.3])
figure('Renderer', 'Painters');
hold on
histogram(Dcc1_m{2},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(3,:),'LineStyle','none')
plot(mean(Dcc1_m{2}),yHgh,'.','color',ccD(3,:),'MarkerSize',12)
plot(median(Dcc1_m{2}),yHgh+.02,'o','color',ccD(3,:),'MarkerSize',8)
axis([0.5 1 0 0.3])

% optim PC applied to 'all' pairs of aPC
yHgh=0.25;
figure('Renderer', 'Painters');
hold on
histogram(Dop_u{3},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(1,:),'LineStyle','none')
plot(mean(Dop_u{3}),yHgh,'.','color',ccD(1,:),'MarkerSize',12)
plot(median(Dop_u{3}),yHgh+.05,'o','color',ccD(1,:),'MarkerSize',8)
axis([0.5 1 0 0.3])
figure('Renderer', 'Painters');
hold on
histogram(Dop_u{1},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(2,:),'LineStyle','none')
plot(mean(Dop_u{1}),yHgh,'.','color',ccD(2,:),'MarkerSize',12)
plot(median(Dop_u{1}),yHgh+.05,'o','color',ccD(2,:),'MarkerSize',8)
axis([0.5 1 0 0.3])
figure('Renderer', 'Painters');
hold on
histogram(Dop_u{2},'Normalization','probability','BinWidth',bw,'BinEdges',edg,'FaceColor',ccD(3,:),'LineStyle','none')
plot(mean(Dop_u{2}),yHgh,'.','color',ccD(3,:),'MarkerSize',12)
plot(median(Dop_u{2}),yHgh+.05,'o','color',ccD(3,:),'MarkerSize',8)
axis([0.5 1 0 0.3])

%% p-vals, effect sizes, for LDA, SCD, THEN Optim

% --- LDA first
[~,pT_NDbic]=ttest2(Dlda_u{3}(:,1),Dlda_u{1}(:,1),'VarType','unequal');
[~,pT_NDmus]=ttest2(Dlda_u{3}(:,1),Dlda_u{2}(:,1),'VarType','unequal');
[~,pT_BICmus]=ttest2(Dlda_u{1}(:,1),Dlda_u{2}(:,1),'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(Dlda_u{3}(:,1),Dlda_u{1}(:,1));
[pW_NDmus,~,zW_NDmus]=ranksum(Dlda_u{3}(:,1),Dlda_u{2}(:,1));
[pW_BICmus,~,zW_BICmus]=ranksum(Dlda_u{1}(:,1),Dlda_u{2}(:,1));
%One-way ANOVA
g_NDbic=[ones(length(Dlda_u{3}(:,1)),1);2*ones(length(Dlda_u{1}(:,1)),1)];
g_NDmus=[ones(length(Dlda_u{3}(:,1)),1);2*ones(length(Dlda_u{2}(:,1)),1)];
g_BICmus=[ones(length(Dlda_u{1}(:,1)),1);2*ones(length(Dlda_u{2}(:,1)),1)];
[pA_NDbic,tblA_NDbic] = anova1([Dlda_u{3}(:,1); Dlda_u{1}(:,1)],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([Dlda_u{3}(:,1) ; Dlda_u{2}(:,1)],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([Dlda_u{1}(:,1) ; Dlda_u{2}(:,1)],g_BICmus,'off');

% table to show ND is stat signif better than Bic/Mus
drgComp={'Bic<ND';'ND<Mus';'Bic<Mus'}; 
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(Dlda_u{3}(:,1)); numBic=length(Dlda_u{1}(:,1)); numMus=length(Dlda_u{2}(:,1));

EffSiz_ttst=[ abs(mean(Dlda_u{3}(:,1))-mean(Dlda_u{1}(:,1)))/sqrt( ((numND-1)*var(Dlda_u{3}(:,1))+(numBic-1)*var(Dlda_u{1}(:,1)))/(numND+numBic-2) ) ; ...
    abs(mean(Dlda_u{3}(:,1))-mean(Dlda_u{2}(:,1)))/sqrt( ((numND-1)*var(Dlda_u{3}(:,1))+(numMus-1)*var(Dlda_u{2}(:,1)))/(numND+numMus-2) ); ...
    abs(mean(Dlda_u{1}(:,1))-mean(Dlda_u{2}(:,1)))/sqrt( ((numBic-1)*var(Dlda_u{2}(:,1))+(numMus-1)*var(Dlda_u{2}(:,1)))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

% --- SCD 
[~,pT_NDbic]=ttest2(Dcc1_m{3},Dcc1_m{1},'VarType','unequal');
[~,pT_NDmus]=ttest2(Dcc1_m{3},Dcc1_m{2},'VarType','unequal');
[~,pT_BICmus]=ttest2(Dcc1_m{1},Dcc1_m{2},'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(Dcc1_m{3},Dcc1_m{1});
[pW_NDmus,~,zW_NDmus]=ranksum(Dcc1_m{3},Dcc1_m{2});
[pW_BICmus,~,zW_BICmus]=ranksum(Dcc1_m{1},Dcc1_m{2});
%One-way ANOVA
g_NDbic=[ones(length(Dcc1_m{3}),1);2*ones(length(Dcc1_m{1}),1)];
g_NDmus=[ones(length(Dcc1_m{3}),1);2*ones(length(Dcc1_m{2}),1)];
g_BICmus=[ones(length(Dcc1_m{1}),1);2*ones(length(Dcc1_m{2}),1)];
[pA_NDbic,tblA_NDbic] = anova1([Dcc1_m{3}; Dcc1_m{1}],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([Dcc1_m{3} ; Dcc1_m{2}],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([Dcc1_m{1} ; Dcc1_m{2}],g_BICmus,'off');

% table to show ND is stat signif better than Bic/Mus
drgComp={'Bic<ND';'ND<Mus';'Bic<Mus'}; 
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(Dcc1_m{3}); numBic=length(Dcc1_m{1}); numMus=length(Dcc1_m{2});

EffSiz_ttst=[ abs(mean(Dcc1_m{3})-mean(Dcc1_m{1}))/sqrt( ((numND-1)*var(Dcc1_m{3})+(numBic-1)*var(Dcc1_m{1}))/(numND+numBic-2) ) ; ...
    abs(mean(Dcc1_m{3})-mean(Dcc1_m{2}))/sqrt( ((numND-1)*var(Dcc1_m{3})+(numMus-1)*var(Dcc1_m{2}))/(numND+numMus-2) ) ; ...
    abs(mean(Dcc1_m{1})-mean(Dcc1_m{2}))/sqrt( ((numBic-1)*var(Dcc1_m{2})+(numMus-1)*var(Dcc1_m{2}))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)

% --- OPtim
[~,pT_NDbic]=ttest2(Dop_u{3},Dop_u{1},'VarType','unequal');
[~,pT_NDmus]=ttest2(Dop_u{3},Dop_u{2},'VarType','unequal');
[~,pT_BICmus]=ttest2(Dop_u{1},Dop_u{2},'VarType','unequal');
%WCM rank sum
[pW_NDbic,~,zW_NDbic]=ranksum(Dop_u{3},Dop_u{1});
[pW_NDmus,~,zW_NDmus]=ranksum(Dop_u{3},Dop_u{2});
[pW_BICmus,~,zW_BICmus]=ranksum(Dop_u{1},Dop_u{2});
%One-way ANOVA
g_NDbic=[ones(length(Dop_u{3}),1);2*ones(length(Dop_u{1}),1)];
g_NDmus=[ones(length(Dop_u{3}),1);2*ones(length(Dop_u{2}),1)];
g_BICmus=[ones(length(Dop_u{1}),1);2*ones(length(Dop_u{2}),1)];
[pA_NDbic,tblA_NDbic] = anova1([Dop_u{3}; Dop_u{1}],g_NDbic,'off');
[pA_NDmus,tblA_NDmus] = anova1([Dop_u{3} ; Dop_u{2}],g_NDmus,'off');
[pA_BICmus,tblA_BICmus] = anova1([Dop_u{1} ; Dop_u{2}],g_BICmus,'off');

% table to show ND is stat signif better than Bic/Mus
drgComp={'Bic<ND';'ND<Mus';'Bic<Mus'}; 
Ttest=[pT_NDbic; pT_NDmus; pT_BICmus]; 
WRankSum=[pW_NDbic; pW_NDmus; pW_BICmus];
owanova=[pA_NDbic; pA_NDmus; pA_BICmus];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(Dop_u{3}); numBic=length(Dop_u{1}); numMus=length(Dop_u{2});

EffSiz_ttst=[ abs(mean(Dop_u{3})-mean(Dop_u{1}))/sqrt( ((numND-1)*var(Dop_u{3})+(numBic-1)*var(Dop_u{1}))/(numND+numBic-2) ) ; ...
    abs(mean(Dop_u{3})-mean(Dop_u{2}))/sqrt( ((numND-1)*var(Dop_u{3})+(numMus-1)*var(Dop_u{2}))/(numND+numMus-2) ); ...
    abs(mean(Dop_u{1})-mean(Dop_u{2}))/sqrt( ((numBic-1)*var(Dop_u{2})+(numMus-1)*var(Dop_u{2}))/(numBic+numMus-2) )];
EffSiz_wrst=abs([zW_NDbic.zval/sqrt(numND+numBic); zW_NDmus.zval/sqrt(numND+numMus); zW_BICmus.zval/sqrt(numBic+numMus)]);
EffSiz_owanova=[tblA_NDbic{6}/tblA_NDbic{8}; tblA_NDmus{6}/tblA_NDmus{8}; tblA_BICmus{6}/tblA_BICmus{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)
