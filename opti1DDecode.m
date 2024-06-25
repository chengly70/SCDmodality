function [D]=opti1DDecode(resp,stimtype)
%This function computes the decoding accuracy given some measured 1-dimensional responses
%to repeated trials of two types of stimulus 
% 'resp' is a vector of length 2n, where n is the number of trials for each
%       type of stimulus, each entry of 'resp' is a measured response 
% 'stimtype' is a vector of length 2n, each entry of 'stimtype' is either a
%       1 or a -1, these label the type of stimulus 
%'D' is the accuracy of the decoder.  It is the fraction of correctly 
%       predicted trials (averaged of all nshuf randomized train/test sets) 

ns=length(stimtype);

%compare order of distances to stimtype to get decoding accuracy (for any shift of the plane along its normal vector)
[~,sind]=sort(resp);
sortstim=stimtype(sind);
convtemp=abs(conv(sortstim,[ones(1,ns) -ones(1,ns)],'full'));
convtemp(:,[1:ns-1 2*ns:end])=[];

%calculate accuracy of decoding
D=max(convtemp)/2/ns+0.5;