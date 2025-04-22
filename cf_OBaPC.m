% compare single cell OB vs aPC decoding accuracy.

%% get single-neuron decod
load IndPCdecode_varyWin_EB_0sSpon.mat

indWind=10; %index from 1,..,20; ALL are bad for PC

ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
edg=(.5-.025:.05:1.025)';
bw=0.1;

DecA_NDpc=[]; DecA_Bicpc=[]; DecA_Muspc=[];
for j=1:8 %total # rats, for ND
    DecA_NDpc=[DecA_NDpc; optThrsh_ND{j,4}(indWind,:).'];
    if j<=3 % bic & mus
        DecA_Bicpc=[DecA_Bicpc; optThrsh_Bic{j,4}(indWind,:).'];
        DecA_Muspc=[DecA_Muspc; optThrsh_Mus{j,4}(indWind,:).'];
    end
    if j==4 % bic
        DecA_Bicpc=[DecA_Bicpc; optThrsh_Bic{j,4}(indWind,:).'];
    end
end

load IndOBdecode_varyWin_EB_0sSpon.mat
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


%% p-vals & effect sizes; showing OB > aPC

% --- LDA first
[~,pT_ND]=ttest2(DecA_ND,DecA_NDpc,'VarType','unequal');
[~,pT_Mus]=ttest2(DecA_Mus,DecA_Muspc,'VarType','unequal');
[~,pT_Bic]=ttest2(DecA_Bic,DecA_Bicpc,'VarType','unequal');
%WCM rank sum
[pW_ND,~,zW_ND]=ranksum(DecA_ND,DecA_NDpc);
[pW_Mus,~,zW_Mus]=ranksum(DecA_Mus,DecA_Muspc);
[pW_Bic,~,zW_Bic]=ranksum(DecA_Bic,DecA_Bicpc);
%One-way ANOVA
g_ND=[ones(length(DecA_ND),1);2*ones(length(DecA_NDpc),1)];
g_Mus=[ones(length(DecA_Mus),1);2*ones(length(DecA_Muspc),1)];
g_Bic=[ones(length(DecA_Bic),1);2*ones(length(DecA_Bicpc),1)];
[pA_ND,tblA_ND] = anova1([DecA_ND;DecA_NDpc],g_ND,'off');
[pA_Mus,tblA_Mus] = anova1([DecA_Mus;DecA_Muspc],g_Mus,'off');
[pA_Bic,tblA_Bic] = anova1([DecA_Bic;DecA_Bicpc],g_Bic,'off');

drgComp={'ND';'Mus';'Bic'}; 
Ttest=[pT_ND; pT_Mus; pT_Bic]; 
WRankSum=[pW_ND; pW_Mus; pW_Bic];
owanova=[pA_ND; pA_Mus; pA_Bic];
T_ndBetter=table(drgComp,Ttest,WRankSum,owanova)

numND=length(DecA_ND); numBic=length(DecA_Bic); numMus=length(DecA_Mus);
numNDpc=length(DecA_NDpc); numBicpc=length(DecA_Bicpc); numMuspc=length(DecA_Muspc);

EffSiz_ttst=[ abs(mean(DecA_ND)-mean(DecA_NDpc))/sqrt( ((numND-1)*var(DecA_ND)+(numNDpc-1)*var(DecA_NDpc))/(numND+numNDpc-2) ) ; ...
    abs(mean(DecA_Mus)-mean(DecA_Muspc))/sqrt( ((numMus-1)*var(DecA_Mus)+(numMuspc-1)*var(DecA_Muspc))/(numMus+numMuspc-2) ); ...
    abs(mean(DecA_Bic)-mean(DecA_Bicpc))/sqrt( ((numBic-1)*var(DecA_Bic)+(numBicpc-1)*var(DecA_Bicpc))/(numBic+numBicpc-2) )];
EffSiz_wrst=abs([zW_ND.zval/sqrt(numND+numNDpc); zW_Mus.zval/sqrt(numMus+numMuspc); zW_Bic.zval/sqrt(numBic+numBicpc)]);
EffSiz_owanova=[tblA_ND{6}/tblA_ND{8}; tblA_Mus{6}/tblA_Mus{8}; tblA_Bic{6}/tblA_Bic{8}];
T_EffSiz_DA=table(drgComp,EffSiz_ttst,EffSiz_wrst,EffSiz_owanova)