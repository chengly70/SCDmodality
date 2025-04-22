% code to plot PSTH for OB trial avg frates in pcolor

clear
% Parameters
odorName='EB';

% Initialize structure
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
        fileName=sprintf('Rat%d_IndCell_%s_%s.mat',i,odorName,drugName);
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
% calc PSTH, avg over trials
mnPSTH_Ret_ND=squeeze(mean(PSTH_Ret_ND,2)); mnPSTH_Or_ND=squeeze(mean(PSTH_Orth_ND,2));
mnPSTH_Ret_Bic=squeeze(mean(PSTH_Ret_Bic,2)); mnPSTH_Or_Bic=squeeze(mean(PSTH_Orth_Bic,2));
mnPSTH_Ret_Mus=squeeze(mean(PSTH_Ret_Mus,2)); mnPSTH_Or_Mus=squeeze(mean(PSTH_Orth_Mus,2));

%Plot
tme=(-0.5:0.1:5)'; %must match with lines 32,34; length is =LastEvok-FirstSpon+1
mnPSTH_Ret_ND=mnPSTH_Ret_ND(16:end,:); mnPSTH_Or_ND=mnPSTH_Or_ND(16:end,:);
mnPSTH_Ret_Bic=mnPSTH_Ret_Bic(16:end,:); mnPSTH_Or_Bic=mnPSTH_Or_Bic(16:end,:);
mnPSTH_Ret_Mus=mnPSTH_Ret_Mus(16:end,:); mnPSTH_Or_Mus=mnPSTH_Or_Mus(16:end,:);
lT=length(tme);
%nmND=size(PSTH_Orth_ND,3); nmBic=size(PSTH_Orth_Bic,3);nmMus=size(PSTH_Orth_Mus,3); %all rats/recordings
nmND=94; nmBic=77; nmMus=140; %size(PSTH_Orth_Mus,3);

%% 
%mnFR=min([mnPSTH_Ret_ND(:);mnPSTH_Or_ND(:);mnPSTH_Ret_Bic(:);mnPSTH_Or_Bic(:);mnPSTH_Ret_Mus(:);mnPSTH_Or_Mus(:)]);
%mxFR=max([mnPSTH_Ret_ND(:);mnPSTH_Or_ND(:);mnPSTH_Ret_Bic(:);mnPSTH_Or_Bic(:);mnPSTH_Ret_Mus(:);mnPSTH_Or_Mus(:)]);
tmp=[reshape(mnPSTH_Ret_ND(:,1:nmND),lT*nmND,1); reshape(mnPSTH_Or_ND(:,1:nmND),lT*nmND,1); ...
    reshape(mnPSTH_Ret_Bic(:,1:nmBic),lT*nmBic,1); reshape(mnPSTH_Or_Bic(:,1:nmBic),lT*nmBic,1); ...
    reshape(mnPSTH_Ret_Mus(:,1:nmMus),lT*nmMus,1); reshape(mnPSTH_Or_Mus(:,1:nmMus),lT*nmMus,1)];
mnFR=min(tmp(:));
mxFR=max(tmp(:));
% get scaled colormap
x = 1:512;
x = x-(.9*mxFR-mnFR)*512/(mxFR-mnFR); %centerPoint=0.9*mxFR
x = .5*x/max(abs(x));  %scaleIntensity=5
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1;
newMap = interp1(x, jet(512), 1:512);
%newMap=turbo;
figure('Renderer', 'Painters');
% pcolor(tme,1:nmND,mnPSTH_Or_ND') %all recordings
subplot(3,2,1)
tmp=mnPSTH_Or_ND(:,1:nmND);
[~,idt]=sort(sum(tmp)); %max(tmp));%
pcolor(tme,1:nmND,tmp(:,idt)')
shading flat
colormap(newMap)
caxis([mnFR mxFR])
%colorbar
subplot(3,2,2)
tmp=mnPSTH_Ret_ND(:,1:nmND);
[~,idt]=sort(sum(tmp));
pcolor(tme,1:nmND,tmp(:,idt)')
shading flat
colormap(newMap)
caxis([mnFR mxFR])
subplot(3,2,3)
tmp=mnPSTH_Or_Bic(:,1:nmBic);
[~,idt]=sort(sum(tmp));
pcolor(tme,1:nmBic,tmp(:,idt)')
shading flat
colormap(newMap)
caxis([mnFR mxFR])
subplot(3,2,4)
tmp=mnPSTH_Ret_Bic(:,1:nmBic);
[~,idt]=sort(sum(tmp));
pcolor(tme,1:nmBic,tmp(:,idt)')
shading flat
colormap(newMap)
caxis([mnFR mxFR])
subplot(3,2,5)
tmp=mnPSTH_Or_Mus(:,1:nmMus);
[~,idt]=sort(sum(tmp));
pcolor(tme,1:nmMus,tmp(:,idt)')
shading flat
colormap(newMap)
caxis([mnFR mxFR])
subplot(3,2,6)
tmp=mnPSTH_Ret_Mus(:,1:nmMus);
[~,idt]=sort(sum(tmp));
pcolor(tme,1:nmMus,tmp(:,idt)')
shading flat
colormap(newMap)
caxis([mnFR mxFR])