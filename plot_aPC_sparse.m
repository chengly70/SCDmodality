% code to plot PSTH for aPC SPARSENESS, pooling over all recordings

clear
% Parameters
odorName='EB';

DrugStruct=cell(3,2); %Rows 1) Or 2) Ret / Columns 1) ND 2) Bic 3) Mus
% Loop over all drugs
for drug_to_keep = 0:2 % DRUG: no drug = 0; Bicu = 1 (less GABA_a inhib); Musc = 2 (more GABA_a inhib)
    %Load good rat data
    if drug_to_keep==0
        ind=[1 2 6 7 8 9 10 11]; %actual good rats
        drugName='NoDrug';
    elseif drug_to_keep==1
        ind=[1 2 8 11]; %actual good rats
        drugName='Bic';
    elseif drug_to_keep == 2
        ind=[6 9 10]; %actual good rats
        drugName='Mus';
    end
    numRats=length(ind);
    AllRat_sRET=[];
    AllRat_sOR=[];
    for i=ind
        fileName=sprintf('Rat%d_IndPC_%s_%s.mat',i,odorName,drugName);
        %Load file, set parameters
        load(fileName) %Ortho/Retro raw spike counts (lenTime, numTrials, nID)
        numTrials=size(sOR,2)+size(sRET,2); %per or/ret cell
        numEvok=TimeVars.numEvok;
        StimShift=TimeVars.StimShift;
        FirstEvok=TimeVars.FirstEvok;
        LastEvok=FirstEvok+49; %TimeVars.LastEvok; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        LastSpon=TimeVars.LastSpon;
        FirstSpon=LastSpon-20; %TimeVars.FirstSpon; %HARD CODED for opt twin-NEEDS CHANGING IF TWIN CHANGED
        lenEvok=length(FirstEvok:1:LastEvok);
        lenSpon=length(FirstSpon:1:LastSpon);
        nOB=size(sOR,3);
        sRetTmp=sRET(FirstSpon:LastEvok,:,:);
        sOrTmp=sOR(FirstSpon:LastEvok,:,:);
        AllRat_sRET=cat(3,AllRat_sRET,sRetTmp);
        AllRat_sOR=cat(3,AllRat_sOR,sOrTmp);
    end %all individual rats
    DrugStruct{drug_to_keep+1,1}=AllRat_sOR;
    DrugStruct{drug_to_keep+1,2}=AllRat_sRET;
end %all drugs
PSTH_Orth_ND=DrugStruct{1,1};
PSTH_Ret_ND=DrugStruct{1,2};
PSTH_Orth_Bic=DrugStruct{2,1};
PSTH_Ret_Bic=DrugStruct{2,2};
PSTH_Orth_Mus=DrugStruct{3,1};
PSTH_Ret_Mus=DrugStruct{3,2};
% calc percent of pop responding
nmND=size(PSTH_Orth_ND,3); nmBic=size(PSTH_Orth_Bic,3); nmMus=size(PSTH_Orth_Mus,3);
PSTH_Ret_ND(PSTH_Ret_ND>0)=1; PSTH_Ret_ND(PSTH_Orth_ND>0)=1; 
PSTH_Ret_ND(PSTH_Ret_Bic>0)=1; PSTH_Ret_ND(PSTH_Orth_Bic>0)=1; 
PSTH_Ret_ND(PSTH_Ret_Mus>0)=1; PSTH_Ret_ND(PSTH_Orth_Mus>0)=1; 
prSpk_Ret_ND=squeeze(sum(PSTH_Ret_ND,3))./nmND; prSpk_Or_ND=squeeze(sum(PSTH_Orth_ND,3))./nmND;
prSpk_Ret_Bic=squeeze(sum(PSTH_Ret_Bic,3))./nmBic; prSpk_Or_Bic=squeeze(sum(PSTH_Orth_Bic,3))./nmBic;
prSpk_Ret_Mus=squeeze(sum(PSTH_Ret_Mus,3))./nmMus; prSpk_Or_Mus=squeeze(sum(PSTH_Orth_Mus,3))./nmMus;

%Plot
tme=(-2:0.1:5)'; %must match with lines 32,34; length is =LastEvok-FirstSpon+1

ccG=.4*ones(1,3);

figure
hold on
plot(tme,mean(prSpk_Or_ND,2),'k','LineWidth',1)
plot(tme,mean(prSpk_Or_ND,2)+std(prSpk_Or_ND,0,2),'k','LineWidth',0.25)
plot(tme,mean(prSpk_Or_ND,2)-std(prSpk_Or_ND,0,2),'k','LineWidth',0.25)
ylabel('Percent of Population Spiking')
xlabel('Time (s)')
title('ND, Orth')
set(gca,'FontSize',18)
axis([-.25 5 0 .53])

figure
hold on
plot(tme,mean(prSpk_Ret_ND,2),'k','LineWidth',1)
plot(tme,mean(prSpk_Ret_ND,2)+std(prSpk_Ret_ND,0,2),'k','LineWidth',.25)
plot(tme,mean(prSpk_Ret_ND,2)-std(prSpk_Ret_ND,0,2),'k','LineWidth',.25)
ylabel('Percent of Population Spiking')
xlabel('Time (s)')
title('ND, Retr')
set(gca,'FontSize',18)
axis([-.25 5 0 .53])
hold off

figure
hold on
plot(tme,mean(prSpk_Or_Bic,2),'k','LineWidth',1)
plot(tme,mean(prSpk_Or_Bic,2)+std(prSpk_Or_Bic,0,2),'k','LineWidth',.25)
plot(tme,mean(prSpk_Or_Bic,2)-std(prSpk_Or_Bic,0,2),'k','LineWidth',.25)
ylabel('Percent of Population Spiking')
xlabel('Time (s)')
title('Bic, Orth')
set(gca,'FontSize',18)
axis([-.25 5 0 .53])
hold off
figure
hold on
plot(tme,mean(prSpk_Ret_Bic,2),'k','LineWidth',1)
plot(tme,mean(prSpk_Ret_Bic,2)+std(prSpk_Ret_Bic,0,2),'k','LineWidth',.25)
plot(tme,mean(prSpk_Ret_Bic,2)-std(prSpk_Ret_Bic,0,2),'k','LineWidth',.25)
ylabel('Percent of Population Spiking')
xlabel('Time (s)')
title('Bic, Retr')
set(gca,'FontSize',18)
axis([-.25 5 0 .53])
hold off

figure
hold on
plot(tme,mean(prSpk_Or_Mus,2),'k','LineWidth',1)
plot(tme,mean(prSpk_Or_Mus,2)+std(prSpk_Or_Mus,0,2),'k','LineWidth',.25)
plot(tme,mean(prSpk_Or_Mus,2)-std(prSpk_Or_Mus,0,2),'k','LineWidth',.25)
ylabel('Percent of Population Spiking')
xlabel('Time (s)')
title('Mus, Orth')
set(gca,'FontSize',18)
axis([-.25 5 0 .53])
hold off
figure
hold on
plot(tme,mean(prSpk_Ret_Mus,2),'k','LineWidth',1)
plot(tme,mean(prSpk_Ret_Mus,2)+std(prSpk_Ret_Mus,0,2),'k','LineWidth',.25)
plot(tme,mean(prSpk_Ret_Mus,2)-std(prSpk_Ret_Mus,0,2),'k','LineWidth',.25)
ylabel('Percent of Population Spiking')
xlabel('Time (s)')
title('Mus, Retr')
set(gca,'FontSize',18)
axis([-.25 5 0 .53])
hold off