function [DFAOut] = DFA_EEG_wrapper(EEG)

DFAOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing DFA exponent...')

Signal = EEG.data;

Signal = transpose(Signal); %now channels are columns, time is rows. Needed for the DFA function

SignalInfo = nbt_Info; %this initializes an Info Object
SignalInfo.converted_sample_frequency = EEG.srate;

DFA_obj = nbt_doDFA(Signal,SignalInfo,[8 round(0.1*length(EEG.data)/EEG.srate)], [0.5 round(length(EEG.data)/EEG.srate)], 0.5, 0, 0, []);

DFAOut = DFA_obj.MarkerValues;
DFAOut = DFAOut';
