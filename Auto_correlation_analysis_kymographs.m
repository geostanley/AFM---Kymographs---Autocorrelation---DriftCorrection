%%%=== Auto_correlation_analysis_kymographs ===%%%

close all
clear variables
clc

%% Enter data variables

Data_directory    = 'Z:\Users\George\Documents\PhD\Data\AFM\Fast Scan\2017_07_25 - Fast Scan - FSD - NIM - Tapping\NPC2\';
File_name_generic = 'NPC2';

File_Nos                = 6:7;
Capture_Direction_First = 'down';
Down_Scan_Only          = 'yes';
Line_Rate_Hz            = 5.23;  
Line_Rate_Save          = round(Line_Rate_Hz); %<-- for string name, not calculation
NS_nine_one             = 'no';
Nupods                  = 'no'; %<-- is this a NuPOD or NPC? 1 for NuPOD
channel                 = 4; %<-- data channel as saved from Nanoscope
Samples_per_line        = 304;
image_width_nm          = 300; 

tau_max_secs                      = 100; %<-- maximum time-lag for calculating auto-coorelation factor, R (usually set in range 80-100)
autocorrelation_over_defined_secs = 1; %<-- for defining max tau, set to 1

% save figures and data to...
Unique_save_name      = 'NPC2_5Hz_retrace';
Save_figure_directory = 'Z:\Users\George\Documents\PhD\Data\Fast_Scanning_Outputs\2017_07_25';

save_outputs = 0; %<-- !!!

% set window for drift correction
drift_shift_window_nm = 10; %<-- refers to total window width

%% Script begins

% make some nice colours for plotting
colourmap = linspecer(124);
C         = linspecer(4);
c1        = C(1,:);

% colour limts for R in nm^2
clims_R = [0 6];
% colour limits for height scale in nm
if length(Nupods) == 3
    clims_z = [-5 15];
else
    clims_z = [-5 60];
end

% Create strings of the files to be loaded
[FileNames, ~] = ImageFuncs.Create_FileNames_Cell_3_spm_pfc(Data_directory, File_name_generic, File_Nos, File_Nos, NS_nine_one); %<-- file numbers twice for spm and pfc files


%%%=== will be getting rid of...
Apply_Gauss_Filt   = 0;
Apply_Fourier_Filt = 0;
ID_nm            = 5;
spatial_cut_off  = 1/ID_nm;
%%%===


%% Load the images and flatten

height_data = ImageFuncs.Load_NS_Height_Data_Kymographs(FileNames, channel, Down_Scan_Only, Capture_Direction_First);

% remove the slope from each image (twice)

BinWidth_nm      = 0.1; % the smaller the value the more accurate the XYZ indexing. 0.1nm should be fine.
Plane_fit_mask_1 = 0.55;
greater_than_1   = 0;

height_data_flat_initial = ImageFuncs.Plane_Subtraction_Images_in_Cell(height_data, BinWidth_nm, Plane_fit_mask_1, greater_than_1);

BinWidth_nm_2    = 0.1; % the smaller the value the more accurate the XYZ indexing. 0.1nm should be fine.
Plane_fit_mask_2 = 0.3;
greater_than_2   = 0;

height_data_flat         = ImageFuncs.Plane_Subtraction_Images_in_Cell(height_data, BinWidth_nm_2, Plane_fit_mask_2, greater_than_2);      

%% Vertically concatonate the images

z_data = height_data_flat{1};

for i = 1:length(height_data_flat)-1
    z_data_i = height_data_flat{i+1};
    z_data   = vertcat(z_data, z_data_i);
end
    
Height_data_before_final_planefit = z_data;

%% Apply plane flattening to contatonated data

% first, through the whole thing...

% transform data into XYZ array
[XYZ_array_entire] = ImageFuncs.Matrix_to_Nx3array(Height_data_before_final_planefit);

% find bottom 30% of height data for plane fitting
Plane_fit_mask = 0;
greater_than   = 1;
BinWidth_nm    = 0.1;
[XYZ_array_entire_for_plane_fit, ~] = ImageFuncs.XYZarray_indexed_by_percentage_height(XYZ_array_entire, BinWidth_nm, Plane_fit_mask, greater_than);

% find the plane of the bottom 30% of height data
[plane_entire_concatonated] = ImageFuncs.PlaneFit_XYZarray(Height_data_before_final_planefit, XYZ_array_entire_for_plane_fit);

% subtract plane from data
Height_data_entire_plane_fit = Height_data_before_final_planefit - plane_entire_concatonated;

% then to the bottom 30%...

% transform data into XYZ array
[XYZ_array] = ImageFuncs.Matrix_to_Nx3array(Height_data_entire_plane_fit);

% find bottom 30% of height data for plane fitting
Plane_fit_mask = 0.3;
greater_than   = 0;
BinWidth_nm    = 0.1;
[XYZ_array_for_plane_fit, XYZ_remaining_data] = ImageFuncs.XYZarray_indexed_by_percentage_height(XYZ_array, BinWidth_nm, Plane_fit_mask, greater_than);

% find the plane of the bottom 30% of height data
[plane_concatonated] = ImageFuncs.PlaneFit_XYZarray(Height_data_before_final_planefit, XYZ_array_for_plane_fit);


% subtract plane from data
Height_data = Height_data_entire_plane_fit - plane_concatonated;

%% create fast-scan axis in nm and slow-scan axis in s

[row_t, col_z] = size(Height_data);

% create array for image width
x_axis_nm = linspace(0, image_width_nm, col_z);

% create array for time axis
dt_perLine = 1/Line_Rate_Hz;
time       = linspace(dt_perLine, dt_perLine*row_t, row_t);

%% Plot plane-corrected kymographs

MyPlots.Height_Plane_Scatter(XYZ_remaining_data, XYZ_array_for_plane_fit, plane_concatonated)

MyPlots.SurfaceKymograph(clims_z, x_axis_nm, time, Height_data_before_final_planefit)
title('Origianl height data', 'FontSize', 13)

MyPlots.SurfaceKymograph(clims_z, x_axis_nm, time, Height_data)
title('Height data after final plane subtraction', 'FontSize', 13)

%% Correct for drift in kymograph

[Height_data_shifted, ~, ~] = ImageFuncs.Drift_Correction_Height_Data(Height_data, x_axis_nm, image_width_nm, drift_shift_window_nm);

% save original unshifted height data for plotting later, and change
% Height_data to the drift corrected version
Height_data_before_shift = Height_data;
Height_data_analyse      = Height_data_shifted;

%% Standard deviation of pixels with time

z_shifted_std   = zeros(1, col_z);
z_unshifted_std = zeros(1, col_z);


for i = 1:col_z
    
    pixel_t_z          = Height_data_analyse(:,i); %<-- drift-corrected
    z_shifted_std(1,i) = std(pixel_t_z);    
    
    pixel_t_z            = Height_data_before_shift(:,i); %<-- not drift-corrected
    z_unshifted_std(1,i) = std(pixel_t_z); 
    
end

% subplot kymographs before and after drift-correction with std plots
MyPlots.KymographBeforeAfterWithStdSubplots(clims_z, x_axis_nm, time, Height_data_before_shift, z_unshifted_std, Height_data_analyse,...
    z_shifted_std, Unique_save_name, Save_figure_directory, Line_Rate_Save, channel, save_outputs)

% plot drift-corrected kymograph only
MyPlots.KymographPlot(clims_z, x_axis_nm, time, Height_data_analyse)
title('Drift-corrected kymograph', 'FontSize', 13)

% Plot SD of drift-corrected kymograph only
ylims = [0 6]; %<-- plot limits for std
MyPlots.SimpleLinePlot_Ratio121(ylims, x_axis_nm, z_shifted_std, Unique_save_name, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
title('SD of drift-corrected kymograph', 'FontSize', 13)

% Subplot kymographs: original and drift-corrected
MyPlots.KymographSubplots(clims_z, x_axis_nm, time, Height_data_before_shift, Height_data_analyse, Unique_save_name, Save_figure_directory, Line_Rate_Save, channel, save_outputs)

%% Calculate the auto-correlation factor, R, as a function of time-lag, tau (in secs)

display('Calculating the auto-correlation factors...')
[R_not_normalised_matrix, l_range_t] = Auto_correlation_function_not_normalised(Height_data_analyse, Samples_per_line, Line_Rate_Hz, tau_max_secs, autocorrelation_over_defined_secs, dt_perLine);

% R_normalised_matrix_scaled = imresize(R_not_normalised_matrix, size(Height_data_analyse)); %<-- keep for plotting in another script

%% Height plot and R_{tau} plotted, separately

FigureSaveName = 'kymograph';
MyPlots.KymographPlotColourbar(clims_z, x_axis_nm, time, Height_data_analyse, Unique_save_name, FigureSaveName,...
    Save_figure_directory, Line_Rate_Save, channel, save_outputs)

FigureSaveName = 'linspecer_R';
MyPlots.AutoCorrelationHeatmapColourbar(clims_R, colourmap, x_axis_nm, l_range_t, R_not_normalised_matrix,...
    Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)

% Overlay height plot with R_{tau} intensity plot
FigureSaveName = 'linspecer_Overlay';
[l_range_t_plot_axis_scaled, R_normalised_matrix_scaled] = MyPlots.KymographAutoCorrelationHeatmapOverlay(clims_z, clims_R, colourmap, x_axis_nm, time, Height_data_analyse, l_range_t, R_not_normalised_matrix,...
    Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs);

% Subplot of drift-corrected height with R_{tau}
FigureSaveName = 'linspecer_subplot_height_R';
MyPlots.SubplotKymographAutoCorrelationHeatmap(clims_z, clims_R, colourmap, x_axis_nm, time, Height_data_analyse, l_range_t, R_not_normalised_matrix,...
    Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)


%% Plot R heatmaps on logarithmic scale

FigureSaveName = 'linspecer - log10_R';
log10_l_range_t = MyPlots.AutoCorrelationHeatmapLogarithmicPlot(clims_R, colourmap, dt_perLine, x_axis_nm, l_range_t, R_not_normalised_matrix,...
    Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs);

% sub-plot with kymograph
FigureSaveName = 'linspecer - SubPlot_Height_log10_R';
MyPlots.SubPlotKymographAutoCorrelationHeatmapLogarithmic(clims_z, clims_R, colourmap, dt_perLine, x_axis_nm, time, Height_data_analyse, l_range_t,...
    R_not_normalised_matrix, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)

%% Save out information as a data structure

FullFileOutput = fullfile(Save_figure_directory, strcat('AutoCorrelationR_Output_DataStructure_', Unique_save_name, '_', num2str(Line_Rate_Save), 'Hz_channel_', num2str(channel)));

AutoCorrelationR_Output_DataStructure.x_axis_nm                         = x_axis_nm;
AutoCorrelationR_Output_DataStructure.time                              = time; 
AutoCorrelationR_Output_DataStructure.Height_data_before_final_planefit = Height_data_before_final_planefit;
AutoCorrelationR_Output_DataStructure.Height_data_before_shift          = Height_data_before_shift;
AutoCorrelationR_Output_DataStructure.Height_data_analyse               = Height_data_analyse;
AutoCorrelationR_Output_DataStructure.z_shifted_std                     = z_shifted_std; % st of shifted data
AutoCorrelationR_Output_DataStructure.z_unshifted_std                   = z_unshifted_std; % st of shifted data
AutoCorrelationR_Output_DataStructure.R_normalised_matrix               = R_not_normalised_matrix;
AutoCorrelationR_Output_DataStructure.l_range_t                         = l_range_t;
AutoCorrelationR_Output_DataStructure.R_normalised_matrix_scaled        = R_normalised_matrix_scaled;
AutoCorrelationR_Output_DataStructure.l_range_t_plot_axis_scaled        = l_range_t_plot_axis_scaled;
AutoCorrelationR_Output_DataStructure.log10_l_range_t                   = log10_l_range_t;
AutoCorrelationR_Output_DataStructure.clims_R                           = clims_R;
AutoCorrelationR_Output_DataStructure.clims_z                           = clims_z;

save(strcat(FullFileOutput, '.mat'), 'AutoCorrelationR_Output_DataStructure');

