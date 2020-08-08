function [LZCOut] = LZC_data_wrapper(dat,nsymbols,window,overlap)
% Calculates the Lempel-Ziv complexity in a sliding window
% 
% Required inputs:
%   dat: the original data in EEGLAB format
% 
% Optional inputs: 
%   nsymbols: the number of symbols to decompose the signal into (default =
%       2, using median split)
%   window: the length of the sliding window used in computation (default =
%       2500 data points)
%   overlap: the amount of overlap between successive windows (default =
%       0.9)
%
% Lempel, A., & Ziv, J. (1976). On the Complexity of Finite Sequences. 
% IEEE Transactions on Information Theory, 22(1), 75?81. 


LZCOut = zeros(1,dat.nbchan);

if ~exist('nsymbols','var')
   nsymbols = 2; 
end

if ~exist('window','var')
   window = 2500;
end

if ~exist('overlap','var')
   overlap = 0.9; 
end

disp(' ')
disp('Computing Lempel-Ziv Complexity...')

for c = 1:dat.nbchan
    fprintf([num2str(c) ' ']);
    [~,LZCOut(c)] = lzcomplexity_tramas(dat.data(c,:),'mediana',nsymbols,window,overlap);
end