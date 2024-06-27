%script to plot results of getCalc_AllPairs_decode.m, pairs of decoding
% or getCalc_Trips_decode.m triplets, DEPENDIND ON flag_pairs variable below (line 10)

clear

flag_sv=0; % if=1, will save new mat file (else NO save)

flag_pairs=1; %if=1, consider all 2x2 pairs, else will show 3x3 results

%color convention
ccD=[0 200 0; 140 0 255; 128.01 128.01 128.01]./255;

if(flag_pairs)
    stdThres=0.2;
    load Reslt_Pairs_All.mat
else
    stdThres=0.23;
    load Reslt_Trips.mat
end

%new variables after 2 filterS
d_nsCorr=cell(3,1); %1stCol=optim PC (use Dop_m), 2nd-3rdCol=OBpair, 4-5Col=PCpairID
nsCorr=cell(3,1);
recnmRat_nsCorr=cell(3,1); %use ccaR_m{j}(:,2), then map to Rat #
err_nsCorr=cell(3,1);
% ---- 
d_Rhoc=cell(3,1); %1stCol=optim PC (use Dop_m), 2ndCol=OBpairID (DOBop_m), 3rdCol=PCpairID (DOBop_m)
Rhoc=cell(3,1);
recnmRat_Rhoc=cell(3,1);
err_Rhoc=cell(3,1);

errDiff=cell(3,1);
cXY_mn=cell(3,1);
stdNs=cell(3,1);
for j=1:3
    if(flag_pairs)
        cXY_mn{j}=mean(abs(CrrOBC{j}(:,1:8)),2);
        stdNs{j}=std(CrrOBC{j}(:,1:8),[],2);
    else
        cXY_mn{j}=mean(abs(CrrOBC{j}(:,1:18)),2);
        stdNs{j}=std(CrrOBC{j}(:,1:18),[],2);
    end
    errDiff{j}=(Dop_m{j}-Dcc1_m{j})./(Dop_m{j}-0.5); %as percent, [0, 1]
end

%only get networks that are appropriate
viaB_p=cell(3,1);
viaB_crr=cell(3,1);
for j=1:3
        viaB_p{j} =(ccaR_m{j}(:,1)>=0.75) &  (DOBop_m{j}(:,1)>=0.75) & stdNs{j}<stdThres;
        viaB_crr{j}=(cXY_mn{j}<=.15) &  (DOBop_m{j}(:,1)>=0.75) & stdNs{j}<stdThres;
end

bw=.05;
for j=1:3
    xtmp=cXY_mn{j}(viaB_p{j});
    id=Dop_m{j}(viaB_p{j})>=Dcc1_m{j}(viaB_p{j}); %a VERY VERY tiny fraction where Dop is numerically better than SCD
    xtmp=abs(xtmp(id));
    ytmp=errDiff{j}(viaB_p{j});
    ytmp=ytmp(id);
    nsCorr{j,1}=xtmp;
    err_nsCorr{j,1}=ytmp;
    d_nsCorr{j,1}=[Dop_m{j}(viaB_p{j}) DOBop_m{j}(viaB_p{j},2:5)];
    d_nsCorr{j,1}=d_nsCorr{j,1}(id,:); 
    recnmRat_nsCorr{j,1}=ccaR_m{j,1}(viaB_p{j},2); %recording Number
    recnmRat_nsCorr{j,1}=recnmRat_nsCorr{j,1}(id);
    figure
    plot(xtmp,ytmp+(2e-3)*randn(size(xtmp)),'.','color',ccD(j,:),'MarkerSize',14)
    xlabel('Noise Correl')
    ylabel('Error')
    

    figure
    hold on
    [z,xedg,yedg]=histcounts2(xtmp,ytmp,'Normalization','probability','XBinEdges',(0-bw/2:bw:1+bw/2),'YBinEdges',(0-bw/2:bw:1+bw/2));
    pcolor(xedg(1:end-1)+bw/2,yedg(1:end-1)+bw/2,z')
    shading flat
    xlabel('Noise Correl')
    ylabel('Error')

    %table of stat signif results
    [tmp,ptmp,ltmp,utmp]=corrcoef([xtmp ytmp]);
    [tmp(1,2) ptmp(1,2) ltmp(1,2) utmp(1,2)]

end

for j=1:3
    xtmp=ccaR_m{j}((viaB_crr{j}),1);
    id=Dop_m{j}(viaB_crr{j})>=Dcc1_m{j}(viaB_crr{j});
    xtmp=xtmp(id);
    ytmp=errDiff{j}(viaB_crr{j});
    ytmp=ytmp(id);
    
    Rhoc{j,1}=xtmp;
    err_Rhoc{j,1}=ytmp;
    d_Rhoc{j,1}=[Dop_m{j}(viaB_crr{j}) DOBop_m{j}(viaB_crr{j},2:5)];
    d_Rhoc{j,1}=d_Rhoc{j,1}(id,:);
    recnmRat_Rhoc{j,1}=ccaR_m{j,1}(viaB_crr{j},2); %recording Number
    recnmRat_Rhoc{j,1}=recnmRat_Rhoc{j,1}(id);
    figure%('Renderer', 'Painters');
    hold on
    plot(xtmp,ytmp+(2e-3)*randn(size(xtmp)),'.','color',ccD(j,:),'MarkerSize',14)
    xlabel('Rcc1')
    ylabel('Error')

    figure 
    hold on
    [z,xedg,yedg]=histcounts2(xtmp,ytmp,'Normalization','probability','XBinEdges',(0-bw/2:bw:1+bw/2),'YBinEdges',(0-bw/2:bw:1+bw/2));
    pcolor(xedg(1:end-1)+bw/2,yedg(1:end-1)+bw/2,z')
    shading flat
    xlabel('Rcc1')
    ylabel('Error')

    %table of stat signif results
    [tmp,ptmp,ltmp,utmp]=corrcoef([xtmp ytmp]);
    [tmp(1,2) ptmp(1,2) ltmp(1,2) utmp(1,2)]

end

if(flag_sv)
    % ---- map recording number (1st col) to rat (2nd col) ----
    NDmap=[(1:29)' [ones(3,1);2*ones(3,1);3*ones(3,1);4*ones(5,1);5*ones(4,1);6*ones(4,1);7*ones(4,1);8*ones(3,1)]];
    BICmap=[(1:12)' [ones(3,1);2*ones(3,1);5*ones(4,1);8*ones(2,1)]];
    MUSmap=[(1:12)' [3*ones(3,1);6*ones(3,1);7*ones(6,1)]];
    recnmRat_nsCorr{1}(:,2)=zeros(size(recnmRat_nsCorr{1})); %initalize
    recnmRat_Rhoc{1}(:,2)=zeros(size(recnmRat_Rhoc{1})); %initalize
    for j=1:29
        id=(recnmRat_nsCorr{1}(:,1)==j);
        recnmRat_nsCorr{1}(id,2)=NDmap(j,2);
        id=(recnmRat_Rhoc{1}(:,1)==j);
        recnmRat_Rhoc{1}(id,2)=NDmap(j,2);
    end
    recnmRat_nsCorr{2}(:,2)=zeros(size(recnmRat_nsCorr{2})); %initalize
    recnmRat_Rhoc{2}(:,2)=zeros(size(recnmRat_Rhoc{2})); %initalize
    for j=1:12
        id=(recnmRat_nsCorr{2}(:,1)==j);
        recnmRat_nsCorr{2}(id,2)=BICmap(j,2);
        id=(recnmRat_Rhoc{2}(:,1)==j);
        recnmRat_Rhoc{2}(id,2)=BICmap(j,2);
    end
    recnmRat_nsCorr{3}(:,2)=zeros(size(recnmRat_nsCorr{3})); %initalize
    recnmRat_Rhoc{3}(:,2)=zeros(size(recnmRat_Rhoc{3})); %initalize
    for j=1:12
        id=(recnmRat_nsCorr{3}(:,1)==j);
        recnmRat_nsCorr{3}(id,2)=MUSmap(j,2);
        id=(recnmRat_Rhoc{3}(:,1)==j);
        recnmRat_Rhoc{3}(id,2)=MUSmap(j,2);
    end
    
    if(flag_pairs)
        save filtPairsDecd err_Rhoc err_nsCorr Rhoc d_Rhoc nsCorr d_nsCorr recnmRat_Rhoc recnmRat_nsCorr 
    else
        save filtTripsDecd err_Rhoc err_nsCorr Rhoc d_Rhoc nsCorr d_nsCorr recnmRat_Rhoc recnmRat_nsCorr 
    end
    
end