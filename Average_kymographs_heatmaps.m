%%%=== Average_kymographs_heatmaps ===%%%

% This script loads several R_{tau} heatmaps and stores them in a cell
% array. The user then selects the centre of each kymograph. These are then cropped,
% and the aligned heatmaps are then averaged.

%% Load data structures 

clear variables
close all
clc


% Enter folder and data structures 
LoadFolder = 'Z:\Users\George\Documents\PhD\Data\Fast_Scanning_Outputs\NuPOD_48Nsp1\Not_Normalised';
Output_Folder = 'Z:\Users\George\Documents\PhD\Data\Fast_Scanning_Outputs\AveragedOutputs\Nsp1\NoNorm\Version4';

% Enter data structures to be averaged
FileNames = {'AutoCorrelationR_Output_DataStructure_Pore2_Nsp1_x48_cleaved_5Hz_retrace_not_norm_files_54to66_5Hz_channel_1',...
    'AutoCorrelationR_Output_DataStructure_Pore3_Nsp1_x48_cleaved_5Hz_retrace_not_norm_files_26to34_5Hz_channel_1',...
    'AutoCorrelationR_Output_DataStructure_Pore4_Nsp1_x48_cleaved_5Hz_retrace_not_norm_files_176to182_5Hz_channel_1'};

Unique_save_name = 'Nsp1_5Hz_NotNorm';

save_outputs = 1; %<-- !!
clims_R      = [0 8]; %<-- auto-correlation heatmap colour range
window_nm    = 30; %<-- set size of cropped and averaged window
colourmap    = linspecer(124); %<-- make nice colourmap

%% Load first data structure and save variables

FileName_1     = FileNames{1};
FullFileName_1 = fullfile(LoadFolder, strcat(FileName_1, '.mat'));

load(FullFileName_1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_axis_nm                         = AutoCorrelationR_Output_DataStructure.x_axis_nm;                       
time                              = AutoCorrelationR_Output_DataStructure.time;                               
Height_data_analyse               = AutoCorrelationR_Output_DataStructure.Height_data_analyse;                
z_shifted_std                     = AutoCorrelationR_Output_DataStructure.z_shifted_std;   
R_matrix                          = AutoCorrelationR_Output_DataStructure.R_normalised_matrix;                
l_range_t                         = AutoCorrelationR_Output_DataStructure.l_range_t;                                                            
clims_z                           = AutoCorrelationR_Output_DataStructure.clims_z; 

% ensure all time-lag axes begin with 1/f, and not 0.
dt_plot = time(end)/length(time);
l_range_n = length(l_range_t);
l_range_t = linspace(dt_plot, l_range_t(end), l_range_n);

dt_plot_array    = zeros(size(FileNames));
dt_plot_array(1) = dt_plot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear AutoCorrelationR_Output_DataStructure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load each data structure and save information into cell arrays

% pre-allocate
x_axis_nm_cell           = cell(1, length(FileNames));
time_cell                = cell(1, length(FileNames));
Height_data_analyse_cell = cell(1, length(FileNames));
z_shifted_std_cell       = cell(1, length(FileNames));
R_matrix_cell            = cell(1, length(FileNames));
l_range_t_cell           = cell(1, length(FileNames));

% store information from first data structure
x_axis_nm_cell{1}           = x_axis_nm;
time_cell{1}                = time;
Height_data_analyse_cell{1} = Height_data_analyse;
z_shifted_std_cell{1}       = z_shifted_std;
R_matrix_cell{1}            = R_matrix;
l_range_t_cell{1}           = l_range_t;

for i = 1:length(FileNames)-1 %<-- iterate through remaining data structures
    
    currentfilename = FileNames{i+1};
    fullfilename    = fullfile(LoadFolder, strcat(currentfilename, '.mat'));
    
    load(fullfilename)
    
    % store into cell arrays
    x_axis_nm_cell{i+1}           = AutoCorrelationR_Output_DataStructure.x_axis_nm; 
    Height_data_analyse_cell{i+1} = AutoCorrelationR_Output_DataStructure.Height_data_analyse;  
    z_shifted_std_cell{i+1}       = AutoCorrelationR_Output_DataStructure.z_shifted_std; 
    R_matrix_cell{i+1}            = AutoCorrelationR_Output_DataStructure.R_normalised_matrix; 
        
    time      = AutoCorrelationR_Output_DataStructure.time;
    l_range_t = AutoCorrelationR_Output_DataStructure.l_range_t;
    
    % ensure all time-lag axes begin with 1/f, and not 0.
    dt_plot = time(end)/length(time);
    l_range_n = length(l_range_t);
    l_range_t = linspace(dt_plot, l_range_t(end), l_range_n);
    
    time_cell{i+1}                = time;
    l_range_t_cell{i+1}           = l_range_t;
    dt_plot_array(i+1)            = dt_plot;
    
end

%% Pick the centre of a NuPOD and crop the R heatmap and kymograph

halfwindow_nm = window_nm/2; %<-- usually 15 nm

% pre-allocate cell array of cropped R heatmaps etc.
R_mat_cropped_cell       = cell(1, length(FileNames));
Height_data_cropped_cell = cell(1, length(FileNames));
x_axis_nm_cropped_cell   = cell(1, length(FileNames));

centre_nm_array                   = zeros(1, length(FileNames));
centre_minus_halfwindow_nm_array  = zeros(1, length(FileNames));
centres_plus_halfwindow_nm_array  = zeros(1, length(FileNames));
centre_idx_array                  = zeros(1, length(FileNames));
centre_minus_halfwindow_idx_array = zeros(1, length(FileNames));
centres_plus_halfwindow_idx_array = zeros(1, length(FileNames));

for n = 1:length(FileNames) %<-- for each data structure

    % pull-out kymograph, R heatmap etc
    Height_data = Height_data_analyse_cell{n};
    x_axis_nm   = x_axis_nm_cell{n};
    time        = time_cell{n};
    l_range_t   = l_range_t_cell{n};    
    R_heatmap   = R_matrix_cell{n}; %<-- not yet scaled by time

    % pick two positions (either inside edge of the DNA scaffold) and take the
    % middle position
    centres_idx     = ImageFuncs.InputCentres_Copper(Height_data);
    centres_col_idx = centres_idx(:,2);
    centre_col_idx  = round(sum(centres_col_idx)/2);
    
    close(gcf) %<-- close this figure

    % crop a window around this position
    centre_nm                  = x_axis_nm(centre_col_idx);
    centre_plus_halfwindow_nm  = centre_nm + halfwindow_nm;
    centre_minus_halfwindow_nm = centre_nm - halfwindow_nm;

    centre_plus_halfwindow_abs_nm  = abs(x_axis_nm - centre_plus_halfwindow_nm);
    centre_plus_halfwindow_idx     = find(centre_plus_halfwindow_abs_nm == min(centre_plus_halfwindow_abs_nm), 1);
    centre_plus_halfwindow_nm      = x_axis_nm(centre_plus_halfwindow_idx);
    
    centre_minus_halfwindow_abs_nm  = abs(x_axis_nm - centre_minus_halfwindow_nm);
    centre_minus_halfwindow_idx     = find(centre_minus_halfwindow_abs_nm == min(centre_minus_halfwindow_abs_nm), 1);
    centre_minus_halfwindow_nm      = x_axis_nm(centre_minus_halfwindow_idx);
    
    % save centre and window edge positions (idx and nm)
    centre_nm_array(n)                   = centre_nm;
    centre_minus_halfwindow_nm_array(n)  = centre_minus_halfwindow_nm;
    centres_plus_halfwindow_nm_array(n)  = centre_plus_halfwindow_nm;
    centre_idx_array(n)                  = centre_col_idx;
    centre_minus_halfwindow_idx_array(n) = centre_minus_halfwindow_idx;
    centres_plus_halfwindow_idx_array(n) = centre_plus_halfwindow_idx;

    % crop height data, R heatmaps, and x_axis_nm and save in cell arrays
    Height_data_cropped     = Height_data(:, centre_minus_halfwindow_idx:centre_plus_halfwindow_idx);
    R_mat_cropped           = R_heatmap(:, centre_minus_halfwindow_idx:centre_plus_halfwindow_idx);
    x_axis_nm_pos_cropped   = x_axis_nm(1:(centre_plus_halfwindow_idx - centre_minus_halfwindow_idx)+1);
            
    x_axis_nm_middle  = max(x_axis_nm_pos_cropped)/2;
    x_axis_nm_cropped = x_axis_nm_pos_cropped - x_axis_nm_middle;
        
    Height_data_cropped_cell{n} = Height_data_cropped;
    R_mat_cropped_cell{n}       = R_mat_cropped;
    x_axis_nm_cropped_cell{n}   = x_axis_nm_cropped;
    
    % linear plot
    fig_name = strcat(Unique_save_name, '_subplot_cropped_individual_Height_R_file_', num2str(n));
    figure_save_name = fullfile(Output_Folder, fig_name);
    R_tau_max = 80;
    MyPlots.SubplotKymographAutoCorrelationHeatmap_simplesave(clims_z, clims_R, halfwindow_nm, dt_plot, R_tau_max, colourmap,...
        x_axis_nm_cropped, time, Height_data_cropped, l_range_t, R_mat_cropped, figure_save_name, save_outputs)
        
    % logarithmic plot
    fig_name = strcat(Unique_save_name, '_subplot_cropped_individual_Height_log10_R_file_', num2str(n));
    figure_save_name = fullfile(Output_Folder, fig_name);
    MyPlots.SubPlotKymographAutoCorrelationHeatmapLogarithmic_simplesave(clims_z, clims_R, colourmap, halfwindow_nm, dt_plot, R_tau_max, x_axis_nm_cropped,...
        time, Height_data_cropped, l_range_t, R_mat_cropped, figure_save_name, save_outputs)
    
end

%% Normalise by time contribution

% 1st, find longest time contribution
time_max_array = zeros(size(FileNames));

for i = 1:length(FileNames)
    
    time_array        = time_cell{i};
    time_max_i        = max(time_array);
    time_max_array(i) = time_max_i;
    
end

time_max_array_sum       = sum(time_max_array');
time_normalisation_array = time_max_array ./ time_max_array_sum; %<-- normalisation scalar for each heatmap

%% now go through the R_mats and scale by their time_normalisation value

R_mat_crop_norm_by_t_cell_t = cell(size(R_mat_cropped_cell));
R_mat_crop_nan_counter      = cell(size(R_mat_cropped_cell));
R_norm_mat_running_ave      = zeros(size(R_mat_cropped_cell{1}));

for i = 1:length(FileNames) %<-- for each heatmap
    
    time_norm_val       = time_normalisation_array(i); %<-- scalar
    R_mat               = R_mat_cropped_cell{i};   
    R_mat_crop_t_norm   = R_mat * time_norm_val; %<-- scaled heatmap
    
    R_mat_crop_norm_by_t_cell_t{i} = R_mat_crop_t_norm; %<-- save individually scaled heatmaps
    R_norm_mat_running_ave         = R_norm_mat_running_ave + R_mat_crop_t_norm; %<-- make the running average
    
end

R_norm_mat_averaged = R_norm_mat_running_ave;

%% Using the centred, normalised (by time contribution) R heatmaps, radially bin data

centre_idx               = ceil(length(x_axis_nm_cropped)/2); %<-- get central index value
x_axis_nm_cropped_radial = x_axis_nm_cropped(centre_idx:end); %<-- make new x_axis_nm from centre (0 nm) outwards

% pre-allocate
R_mat_crop_norm_by_t_radial = zeros(length(l_range_t), length(x_axis_nm_cropped_radial));
R_mat_crop_norm_by_t_radial_averaged = zeros(length(l_range_t), length(x_axis_nm_cropped_radial));
R_mat_crop_norm_by_t_radial_cell = cell(size(R_mat_crop_norm_by_t_cell_t));

% make new R_tau matrices from 0nm to 15 nm by binning
for i = 1:length(FileNames) %<-- for all R_tau heatmaps
    
    Rmat_normalised = R_mat_crop_norm_by_t_cell_t{i}; %<-- pull-out heatmap
    
    Centre_nm_R_mat_column_vec    = Rmat_normalised(:, centre_idx); %<-- 0nm data
    Negative_nm_R_mat_column_vecs = fliplr(Rmat_normalised(:, 1:(centre_idx-1))); %<-- -15nm up to 0nm (flipped)
    Positive_nm_R_mat_column_vecs = Rmat_normalised(:, (centre_idx+1):end); %<-- after 0 nm to 15 nm
    
    % average the negative and positive arrays
    Average_nm_R_mat_column_vecs = (Negative_nm_R_mat_column_vecs + Positive_nm_R_mat_column_vecs) ./ 2;
    
    % put radially binned and averaged data into new matrix
    R_mat_crop_norm_by_t_radial(:, 1) = Centre_nm_R_mat_column_vec;
    R_mat_crop_norm_by_t_radial(:, 2:end) = Average_nm_R_mat_column_vecs;
    
    % save into cell array and keep the running average
    R_mat_crop_norm_by_t_radial_cell{i} = R_mat_crop_norm_by_t_radial;
    R_mat_crop_norm_by_t_radial_averaged = R_mat_crop_norm_by_t_radial_averaged + R_mat_crop_norm_by_t_radial;
    
end
    
    
%% Plot averaged results

fig_name = strcat(Unique_save_name, '_R_radiallybinned_8nm2', num2str(n));
figure_save_name = fullfile(Output_Folder, fig_name);
MyPlots.PlotAveragedAutoCorrelationHeatmap_simplesave(clims_R, halfwindow_nm, dt_plot, R_tau_max, colourmap, x_axis_nm_cropped_radial, time, l_range_t, R_mat_crop_norm_by_t_radial_averaged, figure_save_name, save_outputs)
        
fig_name = strcat(Unique_save_name, '_R_log_radiallybinned_8nm2', num2str(n));
figure_save_name = fullfile(Output_Folder, fig_name);
MyPlots.PlotAveragedAutoCorrelationHeatmapLogarithmic_simplesave(clims_R, colourmap, halfwindow_nm, dt_plot, R_tau_max, time, x_axis_nm_cropped_radial, l_range_t, R_mat_crop_norm_by_t_radial_averaged, figure_save_name, save_outputs)

total_minutes_scanning = time_max_array_sum/60;
%% Save output as data structure

FullFileOutput = fullfile(Output_Folder, strcat('AutoCorrelationR_Output_DataStructure_Averaged_', Unique_save_name));

AutoCorrelationR_Output_DataStructure_Averaged.AveragedData.x_axis_nm_cropped      = x_axis_nm_cropped;
AutoCorrelationR_Output_DataStructure_Averaged.AveragedData.l_range_t              = l_range_t;
AutoCorrelationR_Output_DataStructure_Averaged.AveragedData.R_norm_mat_running_ave = R_norm_mat_running_ave;
AutoCorrelationR_Output_DataStructure_Averaged.AveragedData.time_max_array_sum     = time_max_array_sum;
AutoCorrelationR_Output_DataStructure_Averaged.AveragedData.window_nm              = window_nm;
AutoCorrelationR_Output_DataStructure_Averaged.AveragedData.Unique_save_name       = Unique_save_name;

AutoCorrelationR_Output_DataStructure_Averaged.IndividualData.time_cell                = time_cell;
AutoCorrelationR_Output_DataStructure_Averaged.IndividualData.Height_data_cropped_cell = Height_data_cropped_cell;
AutoCorrelationR_Output_DataStructure_Averaged.IndividualData.R_norm_mat_cropped_cell  = R_mat_cropped_cell;
AutoCorrelationR_Output_DataStructure_Averaged.IndividualData.x_axis_nm_cropped_cell   = x_axis_nm_cropped_cell;
AutoCorrelationR_Output_DataStructure_Averaged.IndividualData.time_max_array           = time_max_array;

save(strcat(FullFileOutput, '.mat'), 'AutoCorrelationR_Output_DataStructure_Averaged');
