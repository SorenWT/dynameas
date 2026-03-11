function hr = Heartrate_hr_wrapper(hrdat)
% Computes the heart rate for heart-rate data
%
% Required inputs:
%   hrdat: the output of a heart-rate transform


evtypes = {hrdat.event.type};
evlats = [hrdat.event.latency];

% take only the IBIs that correspond to heartbeat-heartbeat differences -
% ignore if there's a boundary event in between so that we don't track IBIs
% across excluded segments
for i = 1:(length(evtypes)-1)
    ibis(i) = evlats(i+1)-evlats(i);
    ibi_type{i} = [evtypes{i} '-' evtypes{i+1}];
end
%beat_ind = evlats(strcmpi(evtypes,'heartbeat'));

hr = 60/(median(ibis(strcmpi(ibi_type,{'heartbeat-heartbeat'}))./hrdat.srate));