function [bestD, randD, bestN, bestang]=bruteDecode(rp,stimtype)
%this function does brute force identification of optimal decoding using a
%plane to separate responses to multiple trials of two different stimuli

% input: 'rp' is a (num trials)x(num neurons) matrix.  Each matrix entry
% is one spike count response to one trial of stimulus.

% input: 'stimtype' is (num trials)x1 matrix.  Each matrix entry should be
% +1 or -1, which are labels for the type of stimulus presented on each
% trial.

% output: 'bestD' is the optimal decoding accuracy (fraction of trials
% correctly separated by a 

% output: 'randD' is the set of all decoding accuracies for all separation
% planes.  It can be helpful for comparing to 'bestD' as a control

% output: 'bestN' this is the normal vector for the best decoding plane

% output: 'bestang' is the angle defining the best decoding plane. (same
% information as the normal vector, just in a different form.)

ns=size(rp,1);

%center of mass (com) response
comr=mean(rp,1);

%responses relative to com
relr=rp-comr;

%whats the dimension of the data?
totdim=size(rp,2);

%dimension of the classifiction plane
pdim=totdim-1;

%do brute-force 2d optimal
if totdim==2
    da=pi/500; %da is the angular resolution for list of all angles to try
else
    da=pi/40;
end
anglist=[da:da:pi]';
numang=length(anglist);

%list of all possible angles for a classification plane
angcomb=anglist(fullfact(numang*ones(1,pdim)));

% make normal vectors for all possible classification planes
N=[cos(angcomb(:,1)) zeros(numang^pdim,totdim-2) prod(sin(angcomb(:,1:pdim)),2)];
for d=2:totdim-1
    N(:,d)=prod([sin(angcomb(:,1:d-1)) cos(angcomb(:,d))],2);
end

%response distances from each plane
rdis=N*relr';

%compare order of distances to stimtype to get decoding accuracy (for any shift of the plane along its normal vector)
[~,sind]=sort(rdis,2);
sortstim=stimtype(sind);
convtemp=abs(convn(sortstim,[ones(1,ns) -ones(1,ns)],'full'));
convtemp(:,[1:ns-1 2*ns:end])=[];

%calculate accuracy of decoding
temp=max(convtemp,[],2)/2/ns+0.5;
[bestD, angind]=max(temp);

%also collect Ds from other angles for control 
randD=temp;

bestang=angcomb(angind,:);
bestN=N(angind,:);
