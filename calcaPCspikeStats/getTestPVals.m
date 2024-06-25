% v is a matrix of p-values organized as follows:
    %             c1: p_ttst   c2: p_wrst   c3: p_owa
    % r1: NDxBic
    % r2: NDxMus
    % r3: BicxMus 

function v = getTestPVals(A)
    % Two-sample T-test assuming unequal var (for NDxBic, NDxMus, BicxMus)
    [~,p_ttstNDvBic]=ttest2(A{3},A{1},'VarType','unequal');
    [~,p_ttstNDvMus]=ttest2(A{3},A{2},'VarType','unequal');
    [~,p_ttstBicvMus]=ttest2(A{1},A{2},'VarType','unequal');

    % Wilcoxon rank sum
    p_wrstNDvBic=ranksum(A{3},A{1}); 
    p_wrstNDvMus=ranksum(A{3},A{2});
    p_wrstBicvMus=ranksum(A{1},A{2});

    % One-way ANOVA
    lblsNDvBic=[ones(length(A{3}),1);2*ones(length(A{1}),1)]; 
    p_owaNDvBic = anova1([A{3}; A{1}],lblsNDvBic,'off');
    lblsNDvMus=[ones(length(A{3}),1);2*ones(length(A{2}),1)];
    p_owaNDvMus = anova1([A{3}; A{2}],lblsNDvMus,'off');
    lblsBicvMus=[ones(length(A{1}),1);2*ones(length(A{2}),1)];
    p_owaBicvMus = anova1([A{1}; A{2}],lblsBicvMus,'off');
    
    v = [p_ttstNDvBic p_wrstNDvBic p_owaNDvBic;
        p_ttstNDvMus p_wrstNDvMus p_owaNDvMus;
        p_ttstBicvMus p_wrstBicvMus p_owaBicvMus];
end
