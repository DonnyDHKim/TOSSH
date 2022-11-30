%% Add directories to path
% We navigate to the TOSSH directory and add it to the Matlab path. 
% If we already are in this directory, we can use the pwd command:
mydir = pwd;
% Alternatively, we can specify my_dir manually:
% mydir = 'D:/Sebastian/Documents/MATLAB/TOSSH';
cd(mydir)
addpath(genpath(mydir));

%% Load data
path = './data/'; % specify path

name_list = {
    'berkeley_west-virginia_5894384_01616500',
    'caldwell_texas_1631587_14190500',
    'douglas_colorado_191739_6709000',
    'polk_oregon_23762661_14190500',
    'stark_ohio_19389766_03118500',
    'travis_texas_5781369_08159000'
    };

%% Import data
src_data = readtable(string(append(path, name_list(2),'/USGS_6hrly.csv'))); Q = src_data{:,3}; T = src_data{:,2};
src_data2 = readtable(string(append(path, name_list(2),'/RawFlow_', 'AP.csv'))); Q2 = src_data2{:,7}; T2 = src_data2{:,9};


[MRC_num_segments_USGS , Segment_slopes_USGS] = sig_MRC_SlopeChanges(Q, T, 'plot_results', true, 'recession_length', 5, 'eps', median(Q)*0.5);
[MRC_num_segments_NWM , Segment_slopes_NWM] = sig_MRC_SlopeChanges(Q2, T2, 'plot_results', true, 'recession_length', 5, 'eps', min(Q2)*0.1);

%%
recession_length_list = {5, 5, 5, 5, 5, 5};
eps_list = {min(Q)*0.025, min(Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05}