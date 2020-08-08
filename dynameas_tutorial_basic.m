%% 
% Welcome to the tutorial for the Dynameas toolbox. This tutorial will cover 
% the basic usage of the toolbox, such as how to apply a set of measures to a 
% dataset using dm_applymeasure, how to compare two groups using dm_measurestatistics, 
% and how to visualize the results using dm_measurestatplot. 
% 
% % 
% *Setting up*
% 
% _EEGLAB and Fieldtrip_
% 
% Dynameas requires EEGLAB and Fieldtrip to run. You can download them here
% 
% Fieldtrip: <http://www.fieldtriptoolbox.org/download/ http://www.fieldtriptoolbox.org/download/>
% 
% EEGLAB: <https://sccn.ucsd.edu/eeglab/download.php https://sccn.ucsd.edu/eeglab/download.php>
% 
% Follow their respective setup instructions to make sure EEGLAB and Fieldtrip 
% are properly set up and added to the path. 
% 
% Next, download Dynameas from (Github link). 
% 
% Make sure Dynameas has been added to your path with all its subfolders. To 
% do this, enter the following commands at the prompt, and selecting the Dynameas 
% folder from the pop-up dialog.

dynameas_path = uigetdir;
addpath(genpath(dynameas_path))
savepath
%% 
% Finally, make sure you've downloaded the tutorial dataset. The tutorial dataset 
% consists of EEG data in eyes-closed resting state from a group of 14 healthy 
% individuals and 14 individuals with schizophrenia. The data comes from a study 
% by Olejarczyk and Jernajcyk (2017). The data has been preprocessed with the 
% HAPPE pipeline (Gabard-Durnam et al., 2018), and can be downloaded here (https://drive.google.com/open?id=1kByFVBXKVCYjUitnaaGbuoqiYTRZr7DJ). 
% The preprocessed data are in EEGLAB format.
% 
% For the tutorial below to work, select the tutorial dataset folder using the 
% dialog that pops up below:

tutorial_path = uigetdir;
%% 
% *Getting started*
% 
% Dynameas is a toolbox to apply measures of neural dynamics to whole datasets 
% at once. It also provides functionality to do statistics on these datasets, 
% and visualize the results in a number of ways. Dynameas works with _configuration 
% structures_, as in Fieldtrip. These are Matlab structure arrays (like dictionaries 
% in Python), which contain named fields with parameters that specify the analysis 
% you're going to perform. These structures are usually called _cfg_, but you 
% can name them anything that you want. 
% 
% We'll start by pointing the config structure to the data, and inputing the 
% data format. Let's start with the healthy controls. 
% 
% First, we'll make an empty config structure. 

cfg = struct;
%% 
% Next, we input the directory where the files of interest are found. 

cfg.dir = fullfile(tutorial_path,'healthy');
%% 
% tutorial_path is the path to the folder with the tutorial dataset, which you 
% selected above.
% 
% Next, we set the format. Currently, Dynameas can handle EEGLAB or Fieldtrip 
% formatted datasets. The tutorial dataset is in EEGLAB format, so

cfg.format = 'eeglab';
%% 
% Note that if your dataset is epoched, the epochs will be automatically concatenated 
% into a single, continuous recording. This is the case for the tutorial dataset. 
% 
% Next, we want to select which files we're interested in. We do this using 
% the wildcard character *

cfg.files = '*.set';
%% 
% This will select only the files which end in '.set' (the EEGLAB file extension). 
% 
% % 
% *Extension: *Use of cfg.files and cfg.dir is the means by which one can flexibly 
% select a variety of datasets, even if the relevant files are not all in the 
% same folder. For example, one could select both healthy controls and schizophrenia 
% participants at once by setting cfg.dir = tutorial_path and cfg.files = '*/*.set' 
% - this extra */ in cfg.files will tell Dynameas to look for all files ending 
% in .set _in all subdirectories of _my_path_to_tutorial_dataset. This is a good 
% way to use Dynameas with BIDS-formatted datasets. Alternatively, one could have 
% had healthy controls and schizophrenia participants in one folder to start with. 
% If this was the case, you could set cfg.dir = tutorial_path and cfg.files = 
% 'h*.set' to select only the healthy controls (as all their file names start 
% with h and end in .set).
% 
% % 
% Next, we'll name the output file where we want the results to be saved. 

cfg.outfile = fullfile(tutorial_path,'tutorial_measures_healthy.mat');
%% 
% % 
% Finally, let's select a few measures to apply. All the available measures 
% are found in the 'wrappers' folder, and you can look at each one's documentation 
% by using the help command. For this tutorial, we'll use the power-law exponent 
% (PLE) and autocorrelation window (ACW).
% 
% Measures are defined in Dynameas using _fuction handles_. Function handles 
% are pointers to functions in Matlab (you can read more about them here: <https://www.mathworks.com/help/matlab/matlab_prog/creating-a-function-handle.html 
% https://www.mathworks.com/help/matlab/matlab_prog/creating-a-function-handle.html> 
% ). Function handles all start with the @ symbol. For example, a function handle 
% for ACW would look like:

@ACW_data_wrapper
%% 
% Some measures, like PLE, have multiple inputs - for example, the raw data, 
% plus a frequency range. With function handles, you can create a handle which 
% accepts only one input (in our case, the actual EEG data) and pre-specifies 
% the others. This is what we'll do with the PLE. Since we didn't lowpass the 
% data during preprocessing, we'll take the frequency range as 1 to 100 Hz. 

@(EEG)PLE_data_wrapper(EEG,[1 100])
%% 
% The @(EEG) part specifies that 'EEG' is the only input which will actually 
% be given to the function. The (EEG,[1 100]) part specifies that the second argument 
% (frequency range) will always be [1 100], while the EEG argument is the one 
% we input.
% 
% So finally, to complete the config structure, we make a cell array of these 
% two function handles:

cfg.measure = {@ACW_data_wrapper @(EEG)PLE_data_wrapper(EEG,[1 100])};
%% 
% That's it! We're now ready to run Dynameas. We use dm_applymeasure to apply 
% this set of measures to the healthy controls

healthy_outputs = dm_applymeasure(cfg);
%% 
% Note that the outputs structure healthy_outputs is also saved to the file 
% we specified in cfg.outfile. 
%% 
% *The Outputs structure*
% 
% Dynameas outputs another structure, which has a number of different useful 
% fields. Let's take a quick look at each one.

healthy_outputs
%% 
% The fields are dimord, meas, startsub, chan, chanlocs, elec, data, sub, and 
% cfg. 

healthy_outputs.data
%% 
% The most important part of the outputs structure is the 'data' field. This 
% contains the actual values of the measures you calculated. This is usually a 
% 2D or 3D matrix, depending on whether you've applied multiple measures at once. 

healthy_outputs.dimord
%% 
% dimord specifies the order of the dimensions of outputs.data. 'sub_chan_meas' 
% means that the first dimension is subjects, the second dimension is channels, 
% and the third dimension is measures. The order of the measures is the same as 
% we input them in cfg.measure - so ACW first, then PLE. For example, if I wanted 
% the PLE value for subject 10 at channel 4, I would index outputs.data like so:

healthy_outputs.data(10,4,2)
%% 
% 
healthy_outputs.meas
%% 
% This is essentially a carbon copy of cfg.measure. It's included for use by 
% the statistics and plotting functions.

healthy_outputs.chan
%% 
% This contains the channel labels of the input data. 

healthy_outputs.chanlocs
%% 
% This is the EEGLAB chanlocs structure for the data. It's included for use 
% by later functions. The 'elec' or 'grad' field is the equivalent Fieldtrip channel 
% locations structure. 

healthy_outputs.startsub
%% 
% Usually, this is just 1, unless you've had to break up your processing of 
% the data for some reason. Ignore it for now.

healthy_outputs.sub
%% 
% This contains the file names for each file loaded by Dynameas, in the same 
% order as the 'sub' dimension of outputs.data. This makes it easier to correlate 
% Dynameas data with external variables. 

healthy_outputs.cfg
%% 
% Finally, the original config structure is included with outputs for posterity 
% and debugging. 
%% 
% *Comparing two groups*
% 
% We've now successfully calculated PLE and ACW for the healthy controls of 
% the tutorial dataset. What we're really after, though, is seeing if there are 
% differences in these measures between schizophrenia patients and healthy controls. 
% With that in mind, let's calculate the same measures on the schizophrenia patients 
% and compare. 
% 
% First, we'll calculate the same measures on the schizophrenia patients. We 
% can use the same config structure as for the healthy controls, we just have 
% to change the directory and the output file.

cfg.dir = fullfile(tutorial_path,'schiz');
cfg.outfile = fullfile(tutorial_path,'tutorial_measures_schiz.mat');
%% 
% Now we can apply the same measures to the schizophrenic participants.

schiz_outputs = dm_applymeasure(cfg);
%% 
% Now we have two outputs structures, one for the healthy participants, and 
% another for the schizophrenic participants. We can now use dm_measurestatistics 
% to compare these two statistically. For that, we'll need to set up another config 
% structure.

cfg = struct;
%% 
% The two main things to choose when doing statistics in Dynameas are the statistical 
% test to use, and the method of correction for multiple comparisons. By default, 
% Dynameas takes a mass-univariate approach to statistical testing, meaning separate 
% univariate tests are performed at each electrode/sensor, and the multiple comparisons 
% are corrected for using some method. 
% 
% The options for cfg.test include both parametric and non-parametric methods. 
% For the tutorial data, since there are only 14 subjects, we will opt for non-parametric 
% tests. We will use a non-parametric Wilcoxon rank-sum test to compare the two 
% groups (this is the non-parametric equivalent of a two-sample t test). The other 
% options for cfg.test are documented with the help function. 

cfg.test = 'ranksum';
%% 
% There are several options for correcting for multiple comparisons. Dynameas 
% recommends using the _cluster-based permutation test_ (Maris and Oostenveld, 
% 2007) by default. This procedure involves summing the test statistics of adjacent 
% significant electrodes, and comparing the sum to a permutation distribution. 
% This is a flexible way to establish whether a significant difference between 
% two groups is present, without requiring an a priori selection of electrodes. 
% However, it is important to note that the test does NOT give useful spatial 
% information - that is, one cannot infer statistically from a significant cluster 
% of frontal electrodes that frontal regions are most relevant. The reasons for 
% this are beyond the scope of this tutorial, but the reader is referred to Sassenhagen 
% and Draschkow (2019), or to the Fieldtrip documentation on this point (<http://www.fieldtriptoolbox.org/faq/how_not_to_interpret_results_from_a_cluster-based_permutation_test/ 
% http://www.fieldtriptoolbox.org/faq/how_not_to_interpret_results_from_a_cluster-based_permutation_test/>). 
% Other options include false-discovery rate correction (Benjamini and Yekutieli, 
% 2001), Bonferroni-Holm correction, or simply testing the mean value of each 
% measure across all channels (eliminating the multiple comparisons problem altogether). 
% We will use the cluster-based procedure here. 

cfg.multcompare = 'cluster';
%% 
% This is all we really need to input for dm_measurestatistics - the rest is 
% determined automatically based on these two choices. You can look at the documentation 
% of dm_measurestatistics with the help command to see additional options. For 
% now, let's just do the tests. dm_measurestatistics takes two arguments: the 
% cfg structure we just defined, plus a cell array of outputs structures from 
% dm_applymeasure:

stats = dm_measurestatistics(cfg,{healthy_outputs schiz_outputs});
%% 
% *The stats structure*
% 
% stats is a cell array of structures - each cell corresponds to one of the 
% measures we tested. Let's take a brief look through the stats structure, like 
% we did for the outputs structure. 

stats{2}
%% 
% The fields are p, effsize, cluster, and cfg. 

stats{2}.p
%% 
% These are the uncorrected, univariate p values from each electrode. 

stats{2}.effsize
%% 
% Measures of effect size are calculated automatically in Dynameas using the 
% Measures of Effect Size toolbox (Hentschke and Stüttgen, 2011). These are the 
% outputs from mes.m for each channel. 

stats{2}.cluster
%% 
% This is the output from ft_statistics_montecarlo. You can read more of the 
% details in the Fieldtrip documentation. The most important fields are outlined 
% below:
% 
% stats.cluster.prob: The corrected p-values for each channel. 
% 
% stats.cluster.mask: A vector of zeros and ones which indicates whether a given 
% channel is significant at alpha = 0.05. 
% 
% stats.cluster.posclusters or stats.cluster.negclusters: the p values and test 
% statistics for each "cluster" (a cluster is defined as a set of adjacent channels 
% which are significant at the sensor level). 
% 
% *Note: *When doing a two-tailed cluster-based permutation test, Fieldtrip 
% requires specifying an alpha level of 0.025, rather than 0.05, to reflect the 
% fact that two tests are being done. This is rather unintuitive, however, as 
% most two-tailed tests are still compared with an alpha of 0.05, with the multiplication 
% by two happening internally. Dynameas fixes this by default in its cluster statistics, 
% so the probabilities in the cluster field can be interpreted normally (i.e. 
% p < 0.05 is significant at alpha = 0.05). 
%% 
% *Plotting the results*
% 
% Finally, we can visualize the results using dm_measurestatplot. Once again, 
% we'll set up a new config structure. 

cfg = struct;
%% 
% We need to specify the type of data that we're trying to plot - this can be 
% EEG data, MEG data, or source-level data. 

cfg.datatype = 'eeg';
%% 
% Next, we'll specify names for the conditions to be plotted. We're contrasting 
% healthy participants and schizophrenic participants, so

cfg.cond = {'Healthy','Schizophrenia'};
%% 
% Let's specify the name of the measures to plot. We chose the ACW and PLE, 
% so 

cfg.measname = {'ACW','PLE'};
%% 
% The most important thing to specify is the type of plot desired. This can 
% be 'topo', 'violin', or 'combined'. 'Topo' will plot a topographical representation 
% of the measure in each condition, plus a topographical plot of the difference, 
% with significant channels marked. 'Violin' will plot a violin plot where each 
% point is a subject's mean value over all channels. We can also plot both of 
% these together using the option 'combined' (which is specified by default). 
% We'll use the 'combined' option here. 

cfg.plotmode = 'combined';
%% 
% To plot the results, we input both the original data structures and the new 
% stats structure we calculated. dm_measurestatplot will output two plots, one 
% for each of the two measures we tested. 

dm_measurestatplot(cfg,{healthy_outputs schiz_outputs},stats)
%% 
% And that's it! Unfortunately there were no significant results with these 
% measures and this small dataset, but you can play around with the others and 
% see if anything works.
% 
%