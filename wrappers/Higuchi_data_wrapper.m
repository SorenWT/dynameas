function [HiguchiOut] = Higuchi_EEG_handle(EEG)

HiguchiOut = zeros(1,EEG.nbchan);

disp(' ')
disp("Computing Fractal dimension with Higuchi's method...")

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    HiguchiOut(c) = Higuchi_FD(EEG.data(c,:),44);
end