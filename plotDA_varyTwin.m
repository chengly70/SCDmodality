%script to plot decoding accurac for ALL pairs, but varying time window
%from 100ms to 1sec


%flName='Reslt_Pairs_All';
flName='Reslt_Trips';
ccD=[128.01 128.01 128.01; 0 200 0; 140 0 255]./255;
scl=0.25; %scales STD for visualization purposes

tmWn=(0.1:0.1:1)'; 
lT=length(tmWn);
% -- imprt outputs --
mnLDA=zeros(lT,3); %col DrugPrep: ND=1, Bic=2, Mus=3
mnCCA=zeros(lT,3);
mnOpt=zeros(lT,3);
stdLDA=zeros(lT,3);
stdCCA=zeros(lT,3);
stdOpt=zeros(lT,3);
pValLvC=zeros(lT,3); %p-values of diff LDA vs CCA, using WRST
pValCvO=zeros(lT,3); %p-values of diff CCA vs Optim, using WRST
effSLvC=zeros(lT,3); %EffectSize of diff LDA vs CCA, using WRST
effSCvO=zeros(lT,3); %EffectSize of diff CCA vs Optim, using WRST

for j=1:lT
    load([flName,num2str(j)],'Dlda_u','Dop_u','Dcc1_m')

    idRev=Dlda_u{1}(:,1)<0.5; %find ones that are mis-classified
    Dlda_u{1}(idRev)=1-Dlda_u{1}(idRev); %map to 1-LDA
    idRev=Dlda_u{2}(:,1)<0.5; %find ones that are mis-classified
    Dlda_u{2}(idRev)=1-Dlda_u{2}(idRev); %map to 1-LDA
    idRev=Dlda_u{3}(:,1)<0.5; %find ones that are mis-classified
    Dlda_u{3}(idRev)=1-Dlda_u{3}(idRev); %map to 1-LDA

    tmp=[mean(Dlda_u{3}(:,1)) mean(Dlda_u{1}(:,1)) mean(Dlda_u{2}(:,1))];
    mnLDA(j,:)=tmp; %save results
    tmp=[std(Dlda_u{3}(:,1)) std(Dlda_u{1}(:,1)) std(Dlda_u{2}(:,1))];
    stdLDA(j,:)=tmp; %save results

    tmp=[mean(Dcc1_m{3}) mean(Dcc1_m{1}) mean(Dcc1_m{2})];
    mnCCA(j,:)=tmp;
    tmp=[std(Dcc1_m{3}) std(Dcc1_m{1}) std(Dcc1_m{2})];
    stdCCA(j,:)=tmp;

    tmp=[mean(Dop_u{3}) mean(Dop_u{1}) mean(Dop_u{2})];
    mnOpt(j,:)=tmp;
    tmp=[std(Dop_u{3}) std(Dop_u{1}) std(Dop_u{2})];
    stdOpt(j,:)=tmp;

    % calc p-values & effect sizes
    numNDl=length(Dlda_u{3}(:,1)); numBicl=length(Dlda_u{1}(:,1)); numMusl=length(Dlda_u{2}(:,1));
    numNDc=length(Dcc1_m{3}); numBicc=length(Dcc1_m{1}); numMusc=length(Dcc1_m{2});
    numNDo=length(Dop_u{3}); numBico=length(Dop_u{1}); numMuso=length(Dop_u{2}); %calc optim here
    %WilCx rank sum
    [pW_ND1,~,zW_ND1]=ranksum(Dlda_u{3}(:,1),Dcc1_m{3});
    [pW_Bic1,~,zW_Bic1]=ranksum(Dlda_u{1}(:,1),Dcc1_m{1});
    [pW_Mus1,~,zW_Mus1]=ranksum(Dlda_u{2}(:,1),Dcc1_m{2});
    pValLvC(j,:)=[pW_ND1 pW_Bic1 pW_Mus1];
    tmp=abs([zW_ND1.zval/sqrt(numNDl+numNDc); zW_Bic1.zval/sqrt(numBicl+numBicc); zW_Mus1.zval/sqrt(numMusl+numMusc)]);
    effSLvC(j,:)=tmp';

    [pW_ND1,~,zW_ND1]=ranksum(Dop_u{3},Dcc1_m{3});
    [pW_Bic1,~,zW_Bic1]=ranksum(Dop_u{1},Dcc1_m{1});
    [pW_Mus1,~,zW_Mus1]=ranksum(Dop_u{2},Dcc1_m{2});
    pValCvO(j,:)=[pW_ND1 pW_Bic1 pW_Mus1];
    tmp=abs([zW_ND1.zval/sqrt(numNDo+numNDc); zW_Bic1.zval/sqrt(numBico+numBicc); zW_Mus1.zval/sqrt(numMuso+numMusc)]);
    effSCvO(j,:)=tmp';
end

for j=1:3
    figure
    hold on
    plot(tmWn,mnLDA(:,j),'-','color',ccD(j,:),'LineWidth',0.5)
    plot(tmWn,mnLDA(:,j)+scl*stdLDA(:,j),'--','color',ccD(j,:))
    plot(tmWn,mnLDA(:,j)-scl*stdLDA(:,j),'--','color',ccD(j,:))

    plot(tmWn,mnCCA(:,j),'-','color',ccD(j,:),'LineWidth',1.5)
    plot(tmWn,mnCCA(:,j)+scl*stdCCA(:,j),'--','color',ccD(j,:))
    plot(tmWn,mnCCA(:,j)-scl*stdCCA(:,j),'--','color',ccD(j,:))

    plot(tmWn,mnOpt(:,j),'-','color',ccD(j,:),'LineWidth',2)
    plot(tmWn,mnOpt(:,j)+scl*stdOpt(:,j),'--','color',ccD(j,:))
    plot(tmWn,mnOpt(:,j)-scl*stdOpt(:,j),'--','color',ccD(j,:))
    set(gca,'YLim',[.5 .85])
    set(gca,'FontSize',18)
    ylabel('Decoding Accuracy')
    xlabel('Time Window (s)')
end

for j=1:3
    figure
    subplot(2,1,1)
    hold on
    plot(tmWn,pValLvC(:,j),'-','color',ccD(j,:),'LineWidth',0.5)
    plot(tmWn,pValCvO(:,j),'-','color',ccD(j,:),'LineWidth',1.5)
    set(gca,'FontSize',18)
    ylabel('p-values')
    %set(gca,'YScale','log')
    set(gca,'YLim',[0 .01])

    subplot(2,1,2)
    hold on
    plot(tmWn,effSLvC(:,j),'-','color',ccD(j,:),'LineWidth',0.5)
    plot(tmWn,effSCvO(:,j),'-','color',ccD(j,:),'LineWidth',1.5)
    set(gca,'FontSize',18)
    ylabel('Effect Size')
    xlabel('Time Window (s)')
end
