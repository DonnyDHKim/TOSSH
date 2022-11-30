%% Add directories to path
% We navigate to the TOSSH directory and add it to the Matlab path. 
% If we already are in this directory, we can use the pwd command:
mydir = pwd;
% Alternatively, we can specify my_dir manually:
% mydir = 'D:/Sebastian/Documents/MATLAB/TOSSH';
cd(mydir);
addpath(genpath(mydir));

%% Load data
path = './data/'; % specify path

name_list = {
    'douglas_colorado_191739_6709000',
    'caldwell_texas_1631587_08173000',
    'berkeley_west-virginia_5894384_01616500',    
    'travis_texas_5781369_08159000',
    'stark_ohio_19389766_03118500',    
    'polk_oregon_23762661_14190500'
    };

recession_length_list = {5, 5, 5, 5, 5, 5};
%eps_list = {min(Q)*0.025, (Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05}
%eps_list = {min(Q)*0.025, (Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05}

function peak_df = peak_calc(Q, T)
    if 


%%
% AP_NWM comparison outputs.
% How much NWM is different from AP.
%out_combined = {};
AP_NWM_out_diff_tbl = [];
AP_NWM_out_tbl = [];
AP_NWM_rowname_list = [];

scheme_list = {
    'NWM',
    'AP'};

for i = 1:6
    name = string(name_list(i));
    
    Q_cell = {};
    T_cell = {};
    FDC_highslp_tbl = [];
    MRC_midslp_tbl = [];
    RLD_tbl = [];
    
    src_data = readtable(string(append(path, name,'/RawFlow_AP.csv')));
    %src_data = src_data(1:14611,:);
    src_data = src_data(14612:height(src_data),:);
    Q = src_data{:,7};
    eps_list = {min(Q)*0.025, min(Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05};
    
    for j = 1:length(scheme_list)
        %string(append(path, name,'/RawFlow_', scheme_list(j),'.csv'))
        src_data = readtable(string(append(path, name,'/RawFlow_', scheme_list(j),'.csv')));
        src_data = src_data(14612:height(src_data),:);

        Q = src_data{:,7};
        T = src_data{:,9};
        
        FDC_highslp =  sig_FDC_slope(Q,T,'slope_range',[0 0.33],'fitLogSpace',true);
        [~, Segment_slopes] = sig_MRC_SlopeChanges(Q, T, 'recession_length', 5, 'eps', eps_list{i}); %douglas colorado = problem
        MRC_midslp = Segment_slopes(2);
        RLD = sig_RisingLimbDensity(Q,T,'rising_limb_length',2);

        Q_cell{j,1} = Q;
        T_cell{j,1} = T;
        FDC_highslp_tbl = [FDC_highslp_tbl; FDC_highslp];
        MRC_midslp_tbl = [MRC_midslp_tbl; MRC_midslp];
        RLD_tbl = [RLD_tbl; RLD];
        name_short = split(name, "_");
        rowname = append(name_short(1), "_", scheme_list(j))
        AP_NWM_rowname_list = [AP_NWM_rowname_list; rowname];
    end

    basic_out = calc_BasicSet(Q_cell,T_cell); % in struct format
    basic_out2 = struct2table(basic_out); % in table format, but messy
    basic_out3 = basic_out2(:,[1, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25 ,27]); % table format, tidied up.
    FDC_highslp_tbl = array2table(FDC_highslp_tbl);
    MRC_midslp_tbl = array2table(MRC_midslp_tbl);
    RLD_tbl = array2table(RLD_tbl);
    basic_out3 = [basic_out3,FDC_highslp_tbl,MRC_midslp_tbl,RLD_tbl];
    x = array2table(100*2*(basic_out3{1,:} - basic_out3{2,:})./(basic_out3{1,:}+basic_out3{2,:}), 'variablenames', basic_out3.Properties.VariableNames);
    y = array2table(100*2*(basic_out3{2,:} - basic_out3{2,:})./(basic_out3{2,:}+basic_out3{2,:}), 'variablenames', basic_out3.Properties.VariableNames);

    AP_NWM_out_tbl = [AP_NWM_out_tbl; basic_out3];
    AP_NWM_out_diff_tbl = [AP_NWM_out_diff_tbl; x; y];
    
end
AP_NWM_out_tbl.Properties.RowNames = AP_NWM_rowname_list;
AP_NWM_out_diff_tbl.Properties.RowNames = AP_NWM_rowname_list;




%%
% AP_HUC12L comparison outputs.
% How much NWM is different from AP.
%out_combined = {};
AP_HUC12L_out_diff_tbl = [];
AP_HUC12L_out_tbl = [];
AP_HUC12L_rowname_list = [];

scheme_list = {
    'HUC12L',
    'AP'};

for i = 1:6
    name = string(name_list(i));
    
    Q_cell = {};
    T_cell = {};
    FDC_highslp_tbl = [];
    MRC_midslp_tbl = [];
    RLD_tbl = [];
    
    src_data = readtable(string(append(path, name,'/RawFlow_AP.csv')));
    src_data = src_data(14612:height(src_data),:);
    Q = src_data{:,7};
    eps_list = {min(Q)*0.025, min(Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05};
    
    for j = 1:length(scheme_list)
        %string(append(path, name,'/RawFlow_', scheme_list(j),'.csv'))
        src_data = readtable(string(append(path, name,'/RawFlow_', scheme_list(j),'.csv')));
        %src_data = src_data(1:14611,:);
        src_data = src_data(14612:height(src_data),:);
        
        Q = src_data{:,7};
        T = src_data{:,9};
                
        FDC_highslp =  sig_FDC_slope(Q,T,'slope_range',[0 0.33],'fitLogSpace',true);
        [~, Segment_slopes] = sig_MRC_SlopeChanges(Q, T, 'recession_length', 5, 'eps', eps_list{i}); %douglas colorado = problem
        MRC_midslp = Segment_slopes(2);
        RLD = sig_RisingLimbDensity(Q,T,'rising_limb_length',2);

        Q_cell{j,1} = Q;
        T_cell{j,1} = T;
        FDC_highslp_tbl = [FDC_highslp_tbl; FDC_highslp];
        MRC_midslp_tbl = [MRC_midslp_tbl; MRC_midslp];
        RLD_tbl = [RLD_tbl; RLD];
        name_short = split(name, "_");
        rowname = append(name_short(1), "_", scheme_list(j))
        AP_HUC12L_rowname_list = [AP_HUC12L_rowname_list; rowname];
    end

    basic_out = calc_BasicSet(Q_cell,T_cell); % in struct format
    basic_out2 = struct2table(basic_out); % in table format, but messy
    basic_out3 = basic_out2(:,[1, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25 ,27]); % table format, tidied up.
    FDC_highslp_tbl = array2table(FDC_highslp_tbl);
    MRC_midslp_tbl = array2table(MRC_midslp_tbl);
    RLD_tbl = array2table(RLD_tbl);
    basic_out3 = [basic_out3,FDC_highslp_tbl,MRC_midslp_tbl,RLD_tbl];
    x = array2table(100*2*(basic_out3{1,:} - basic_out3{2,:})./(basic_out3{1,:}+basic_out3{2,:}), 'variablenames', basic_out3.Properties.VariableNames);
    y = array2table(100*2*(basic_out3{2,:} - basic_out3{2,:})./(basic_out3{2,:}+basic_out3{2,:}), 'variablenames', basic_out3.Properties.VariableNames);

    AP_HUC12L_out_tbl = [AP_HUC12L_out_tbl; basic_out3];
    AP_HUC12L_out_diff_tbl = [AP_HUC12L_out_diff_tbl; x; y];
    
end
AP_HUC12L_out_tbl.Properties.RowNames = AP_HUC12L_rowname_list;
AP_HUC12L_out_diff_tbl.Properties.RowNames = AP_HUC12L_rowname_list;





%%
% Allocation comparison outputs.
%out_combined = {};
alloc_out_diff_tbl = [];
alloc_out_tbl = [];
alloc_rowname_list = [];
alloc_diff_rowname_list = [];

scheme_list = {
    'NWM',
    'PERT',
    'AP',
    'SHUF'
    };

for i = 1:6
    name = string(name_list(i));
    
    Q_cell = {};
    T_cell = {};
    FDC_highslp_tbl = [];
    MRC_midslp_tbl = [];
    RLD_tbl = [];
    
    src_data = readtable(string(append(path, name,'/RawFlow_AP.csv')));
    src_data = src_data(14612:height(src_data),:);
    Q = src_data{:,7};
    eps_list = {min(Q)*0.025, min(Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05};
    
    for j = 1:length(scheme_list)
        %string(append(path, name,'/RawFlow_', scheme_list(j),'.csv'))
        src_data = readtable(string(append(path, name,'/RawFlow_', scheme_list(j),'.csv')));
        %src_data = src_data(1:14611,:);
        src_data = src_data(14612:height(src_data),:);
        
        Q = src_data{:,7};
        T = src_data{:,9};
                
        FDC_highslp =  sig_FDC_slope(Q,T,'slope_range',[0 0.33],'fitLogSpace',true);
        [~, Segment_slopes] = sig_MRC_SlopeChanges(Q, T, 'recession_length', 5, 'eps', eps_list{i}); %douglas colorado = problem
        MRC_midslp = Segment_slopes(2);
        RLD = sig_RisingLimbDensity(Q,T,'rising_limb_length',2);

        
        Q_cell{j,1} = Q;
        T_cell{j,1} = T;
        FDC_highslp_tbl = [FDC_highslp_tbl; FDC_highslp];
        MRC_midslp_tbl = [MRC_midslp_tbl; MRC_midslp];
        RLD_tbl = [RLD_tbl; RLD];
        
        name_short = split(name, "_");
        rowname = append(name_short(1), "_", scheme_list(j))
        alloc_rowname_list = [alloc_rowname_list; rowname];

    end

    basic_out = calc_BasicSet(Q_cell,T_cell); % in struct format
    basic_out2 = struct2table(basic_out); % in table format, but messy
    basic_out3 = basic_out2(:,[1, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25 ,27]); % table format, tidied up.
    FDC_highslp_tbl = array2table(FDC_highslp_tbl);
    MRC_midslp_tbl = array2table(MRC_midslp_tbl);
    RLD_tbl = array2table(RLD_tbl);
    basic_out3 = [basic_out3,FDC_highslp_tbl,MRC_midslp_tbl,RLD_tbl];
    
    x = array2table(100*2*(basic_out3{2,:} - basic_out3{1,:})./(basic_out3{2,:}+basic_out3{1,:}), 'variablenames', basic_out3.Properties.VariableNames);
    y = array2table(100*2*(basic_out3{4,:} - basic_out3{3,:})./(basic_out3{4,:}+basic_out3{3,:}), 'variablenames', basic_out3.Properties.VariableNames);

    diff_rowname = [append(name_short(1), "_diff_NWM_PERT"); append(name_short(1), "_diff_AP_SHUF")];

    alloc_diff_rowname_list = [alloc_diff_rowname_list; diff_rowname];
        
    %out_combined{i,1} = basic_out;
    alloc_out_diff_tbl = [alloc_out_diff_tbl; x; y];
    alloc_out_tbl = [alloc_out_tbl; basic_out3];     
    
end
alloc_out_tbl.Properties.RowNames = alloc_rowname_list;
alloc_out_diff_tbl.Properties.RowNames = alloc_diff_rowname_list;
%%
%a=[1 11; 2 22; 3 33; 4 44; 5 55]
%b=mean(a)

% save table
writetable(alloc_out_tbl, 'alloc_pt2.csv', 'WriteRowNames',true)
writetable(alloc_out_diff_tbl, 'alloc_diff_pt2.csv', 'WriteRowNames',true)

writetable(AP_NWM_out_tbl, 'AP_NWM_pt2.csv', 'WriteRowNames',true)
writetable(AP_NWM_out_diff_tbl, 'AP_NWM_diff_pt2.csv', 'WriteRowNames',true)

writetable(AP_HUC12L_out_tbl, 'AP_HUC12L_pt2.csv', 'WriteRowNames',true)
writetable(AP_HUC12L_out_diff_tbl, 'AP_HUC12L_diff_pt2.csv', 'WriteRowNames',true)





%%
% AP_HUC12L comparison outputs.
% How much NWM is different from AP.
%out_combined = {};
AP_SHUF_out_diff_tbl = [];
AP_SHUF_out_tbl = [];
AP_SHUF_rowname_list = [];

scheme_list = {
    'SHUF',
    'AP'};

for i = 1:6
    name = string(name_list(i));
    
    Q_cell = {};
    T_cell = {};
    FDC_highslp_tbl = [];
    MRC_midslp_tbl = [];
    RLD_tbl = [];
    
    src_data = readtable(string(append(path, name,'/RawFlow_AP.csv')));
    Q = src_data{:,7};
    eps_list = {min(Q)*0.025, min(Q)*0.1, min(Q)*0.025, min(Q)*0.25, min(Q)*0.001, min(Q)*0.05};
    
    for j = 1:length(scheme_list)
        %string(append(path, name,'/RawFlow_', scheme_list(j),'.csv'))
        src_data = readtable(string(append(path, name,'/RawFlow_', scheme_list(j),'.csv')));
        
        Q = src_data{:,7};
        T = src_data{:,9};
                
        FDC_highslp =  sig_FDC_slope(Q,T,'slope_range',[0 0.33],'fitLogSpace',true);
        [~, Segment_slopes] = sig_MRC_SlopeChanges(Q, T, 'recession_length', 5, 'eps', eps_list{i}); %douglas colorado = problem
        MRC_midslp = Segment_slopes(2);
        RLD = sig_RisingLimbDensity(Q,T,'rising_limb_length',2);

        Q_cell{j,1} = Q;
        T_cell{j,1} = T;
        FDC_highslp_tbl = [FDC_highslp_tbl; FDC_highslp];
        MRC_midslp_tbl = [MRC_midslp_tbl; MRC_midslp];
        RLD_tbl = [RLD_tbl; RLD];
        name_short = split(name, "_");
        rowname = append(name_short(1), "_", scheme_list(j))
        AP_SHUF_rowname_list = [AP_SHUF_rowname_list; rowname];
    end

    basic_out = calc_BasicSet(Q_cell,T_cell); % in struct format
    basic_out2 = struct2table(basic_out); % in table format, but messy
    basic_out3 = basic_out2(:,[1, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25 ,27]); % table format, tidied up.
    FDC_highslp_tbl = array2table(FDC_highslp_tbl);
    MRC_midslp_tbl = array2table(MRC_midslp_tbl);
    RLD_tbl = array2table(RLD_tbl);
    basic_out3 = [basic_out3,FDC_highslp_tbl,MRC_midslp_tbl,RLD_tbl];
    x = array2table(100*2*(basic_out3{2,:} - basic_out3{1,:})./(basic_out3{2,:}+basic_out3{1,:}), 'variablenames', basic_out3.Properties.VariableNames);
    y = array2table(100*2*(basic_out3{2,:} - basic_out3{2,:})./(basic_out3{2,:}+basic_out3{2,:}), 'variablenames', basic_out3.Properties.VariableNames);

    AP_SHUF_out_tbl = [AP_SHUF_out_tbl; basic_out3];
    AP_SHUF_out_diff_tbl = [AP_SHUF_out_diff_tbl; x; y];
    
end
AP_SHUF_out_tbl.Properties.RowNames = AP_SHUF_rowname_list;
AP_SHUF_out_diff_tbl.Properties.RowNames = AP_SHUF_rowname_list;




%%
% testing caldwell MRC
src_data = readtable(string(append(path, name_list(2),'/RawFlow_', 'NWM.csv')));
src_data = readtable(string(append(path, name_list(2),'/RawFlow_', 'AP.csv')));
src_data = readtable(string(append(path, name_list(2),'/RawFlow_', 'HUC12L.csv')));


Q = src_data{:,7}; T = src_data{:,9};

[MRC_num_segments , Segment_slopes] = sig_MRC_SlopeChanges(Q, T, 'plot_results', true, 'recession_length', 5);
%%

[MRC_num_segments , Segment_slopes] = sig_MRC_SlopeChanges(Q, T, 'recession_length', 11);
RLD = sig_RisingLimbDensity(Q,T,'plot_results',true,'rising_limb_length',2);

x = [split(name, "_")]
x = x(1)
append(x,"_",scheme_list(j))