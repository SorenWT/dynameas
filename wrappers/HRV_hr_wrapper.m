function hrv = HRV_hr_wrapper(hrdat,type)
% Calculates heart rate variability for heart-rate data
%
% Required inputs; 
%   hrdat: the output of a heart-rate transform
% 
% Recommended inputs:
%   type: the heart-rate variability measure used (default = SDRR, the SD
%       of the IBIs of all beats) 

if ~exist('type','var')
   type = 'SDRR'; 
end

evtypes = {hrdat.event.type};
evlats = [hrdat.event.latency];

for i = 1:(length(evtypes)-1)
    ibis(i) = evlats(i+1)-evlats(i);
    ibi_type{i} = [evtypes{i} '-' evtypes{i+1}];
end

% take only the IBIs that correspond to heartbeat-heartbeat differences
ibis = ibis(strcmpi(ibi_type,'heartbeat-heartbeat'));
ibis = ibis/hrdat.srate;

% remove unphysiological gaps
ibis(ibis>3) = [];

switch type
    case 'SDRR'
    hrv = std(ibis);
    case 'RMSSD'
    hrv = sqrt(mean(diff(ibis).^2));
end

