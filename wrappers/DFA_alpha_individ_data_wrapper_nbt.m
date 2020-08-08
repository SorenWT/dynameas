function [DFAout] = DFA_alpha_individ_EEG_wrapper_nbt(EEG)

[psum,~,f] = restingIAF(EEG.data,EEG.nbchan,3,[1 40],EEG.srate,[5 15],11,5);

if ~isnan(psum.paf)
    bp = [f(psum.iaw(1)) f(psum.iaw(2))];
else
    bp = [8 13];
end


Signal = EEG.data;

Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function

SignalInfo = nbt_Info; %this initializes an Info Object
SignalInfo.converted_sample_frequency = EEG.srate;

bpfreqs = bp;
AmplitudeEnvelope = nbt_GetAmplitudeEnvelope(Signal, SignalInfo, bpfreqs(1), bpfreqs(2), 2*(1/bpfreqs(1)));

EEG.data = AmplitudeEnvelope';

DFAout = DFA_EEG_wrapper_nbt(EEG);
