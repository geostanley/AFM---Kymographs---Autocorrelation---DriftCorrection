classdef ImageFuncs
% Contains all functions related to loading, flattening and manipulating
% images/kymographs

    methods(Static)
        
        
        %%%=== Create_FileNames_Cell_3_spm_pfc ===%%%

        % This function creates a cell containing strings of filenames. Input the
        % directory contatining the desired data files, the file suffix common to
        % all data files, and an array containing the numbers of all file ends.
        % It creates both spm and pfc file names.

        % Also tell the function which version of Nanoscope the spm files are from.
        % If 9.1, input the string 'yes'. If 9.2 or higher, input 'no'.

        % A cell array with all the file names will be returned.

        function [FileNames_spm, FileNames_pfc] = Create_FileNames_Cell_3_spm_pfc(Data_Directory, File_Suffix, File_Nos_spm, File_Nos_pfc, NS_nine_one)

        File_Numbers_spm = cell(1,length(File_Nos_spm));
        for i=1:length(File_Nos_spm);
            if File_Nos_spm(i) < 10
                File_Numbers_spm{i} = strcat('00', num2str(File_Nos_spm(i)));
            elseif File_Nos_spm(i) >=10 && File_Nos_spm(i) < 100
                File_Numbers_spm{i} = strcat('0', num2str(File_Nos_spm(i)));
            else
                File_Numbers_spm{i} = num2str(File_Nos_spm(i));
            end
        end

        File_Numbers_pfc = cell(1,length(File_Nos_pfc));
        for i=1:length(File_Nos_pfc);
            if File_Nos_pfc(i) < 10
                File_Numbers_pfc{i} = strcat('00', num2str(File_Nos_pfc(i)));
            elseif File_Nos_pfc(i) >=10 && File_Nos_pfc(i) < 100
                File_Numbers_pfc{i} = strcat('0', num2str(File_Nos_pfc(i)));
            else
                File_Numbers_pfc{i} = num2str(File_Nos_pfc(i));
            end
        end


        if length(NS_nine_one) == 3

            FileNames_spm = cell(1,length(File_Nos_spm));
            FileNames_pfc = cell(1,length(File_Nos_pfc));

            for i = 1:length(File_Nos_spm)

                FileNames_spm{i} = horzcat(Data_Directory, File_Suffix,'.', File_Numbers_spm{i});
                FileNames_pfc{i} = horzcat(Data_Directory, File_Suffix,'.', File_Numbers_pfc{i}, '.pfc');

            end

        else

            FileNames_spm = cell(1,length(File_Nos_spm));
            FileNames_pfc = cell(1,length(File_Nos_pfc));

            for i = 1:length(File_Nos_spm)

                FileNames_spm{i} = horzcat(Data_Directory, File_Suffix,'.0_00', File_Numbers_spm{i}, '.spm');
                FileNames_pfc{i} = horzcat(Data_Directory, File_Suffix,'.0_00', File_Numbers_pfc{i}, '.pfc');

            end

        end

        end
        
        
        %%%=== Load_NS_Height_Data_Kymographs ===%%%

        % Load NanoScope AFM height images, and order chronologically in time. If
        % Up and Down scanning, all up scans are flipped such that each lien scan
        % is chronological. This is for making kymographs.


        function [height_data] = Load_NS_Height_Data_Kymographs(FileNames, channel, Down_Scan_Only, Capture_Direction_First)

        % Load the images into MATLAB, and ensure they run from row 1 to row end, chronologically in time
        NSMU = NSMatlabUtilities();

        height_data = cell(1,length(FileNames));
        raw_data    = cell(size(FileNames));


        for i = 1:length(FileNames)

            PFQNM_FileName = FileNames{i};   
            NSMU.Open(PFQNM_FileName);

            % Channel 1 is Trace Height Sensor. Extract this data.
            [data, ~, ~] = NSMU.GetImageData(channel, NSMU.METRIC);
            raw_data{i} = data;

        end

        if length(Down_Scan_Only) == 2

            if length(Capture_Direction_First) == 2 % if first image was taken in Up direction

                % All images taken in Up direction must be flipped. Therefore, 
                % as first image taken in Up direction, all odds to be flipped, 
                % evens left as normal.

                for i = 1:length(raw_data)
                    Z = raw_data{i};
                    if mod(i, 2) == 1 % if odd number
                        Z = flipud(Z);
                    end
                    height_data{i} = Z;
                end
            else
                for i = 1:length(raw_data)
                    Z = raw_data{i};
                    if mod(i, 2) == 0 % if even number
                        Z = flipud(Z);
                    end
                    height_data{i} = Z;
                end
            end 
        else % if Down Scan Only = yes, do not flip any of the data  
            for i = 1:length(raw_data)
                    Z = raw_data{i};       
                    height_data{i} = Z;          
            end    
        end
        end
        
        
        
        
        %%%=== Matrix_to_Nx3array ===%%%
        % This function takes in a matrix (square or rectangular), and converts 
        % it into an Nx3 array. I.e., XYZ coordinates.

        function [XYZ_array] = Matrix_to_Nx3array(matrix)

            % get the dimensions of the cropped image. As square, r=c.
            [r,c] = size(matrix);

            % transform the square matrix into an Nx3 array
            XYZ_array = zeros(length(matrix(:)), 3);

            for i=1:c
                XYZ_array((((i-1)*r+1):i*r), 1) = i;
            end

            for j = 1:c
                XYZ_array((((j-1)*r+1):j*r), 2) = 1:r;
            end

            for i=1:length(XYZ_array)
                XYZ_array(i,3) = matrix(XYZ_array(i,2),XYZ_array(i,1));
            end

        end
               

    
        %%%=== XYZarray_indexed_by_percentage_height ===%%%

        % This function takes in an Nx3 array (i.e., XYZ coordinates), takes a
        % histogram of the data, and only keeps the data which is either above
        % (greater_than == 1), or below (greater_than == 0) the Planefit_mask,
        % which is a number between 0 and 1, and is a percentage. 

        % User must also input the BinWidth in nm. The smaller the better. Value
        % should probably be <0.5 nm.

        % I.e., greater_than = 1 & Plane_fit_mask = 0.4, means keep the top 60% of
        % the height data.

        function [XYZ_array_for_plane_fit, XYZ_array_not_for_plane_fit] = XYZarray_indexed_by_percentage_height(XYZ_array, BinWidth_nm, Plane_fit_mask, greater_than)

            X_array = XYZ_array(:,1);
            Y_array = XYZ_array(:,2);
            Z_array = XYZ_array(:,3);

            % Find the top x% of height data for plane fitting.
            % Take a histogram of data, find the desired cut-off point, and only take
            % that data forward for finding the 1st order plane fit - just like
            % Nanoscope Analysis.

            % take histogram of height data, binning into 1nm bins
            Z_histedges = [floor(min(Z_array)):BinWidth_nm:ceil(max(Z_array))];
            Z_hist_counts = histcounts(Z_array, Z_histedges);
            Z_hist_counts_sum = sum(Z_hist_counts(:));

            % create an array of bin centres in nm
            Bin_centres = zeros(1, length(Z_histedges)-1);
            for i = 1:length(Bin_centres)
                Bin_centres(i) = Z_histedges(i) + (BinWidth_nm/2);
            end

            % create an array with aggregate counts based on bins
            Z_aggregate_hist_counts = zeros(size(Z_hist_counts));
            Z_aggregate_hist_counts(1) = Z_hist_counts(1);
            for i=1:length(Z_hist_counts)-1
                Z_aggregate_hist_counts(i+1) = Z_aggregate_hist_counts(i) + Z_hist_counts(i+1);
            end

            % represent aggregate histcounts as running percent, then find the 
            % cut-off point as an index and in nm
            Z_aggregate_hist_counts_percentage = Z_aggregate_hist_counts ./ Z_hist_counts_sum;
            Z_aggregate_hist_counts_abs = abs(Z_aggregate_hist_counts_percentage - Plane_fit_mask);
            Z_cut_off_idx = find(Z_aggregate_hist_counts_abs == min(Z_aggregate_hist_counts_abs),1);  
            Z_cut_off_nm  = Bin_centres(Z_cut_off_idx);

            % Pull out the data points at, or above the cut off point (nm)
            if greater_than == 1
                Z_array_for_plane_fit_idx = find(Z_array >= Z_cut_off_nm);
                Z_array_not_for_plane_fit_idx = find(Z_array < Z_cut_off_nm);
            else
                Z_array_for_plane_fit_idx = find(Z_array < Z_cut_off_nm);
                Z_array_not_for_plane_fit_idx = find(Z_array >= Z_cut_off_nm);
            end

            % save desired data
            X_array_for_plane_fit = X_array(Z_array_for_plane_fit_idx);
            Y_array_for_plane_fit = Y_array(Z_array_for_plane_fit_idx);
            Z_array_for_plane_fit = Z_array(Z_array_for_plane_fit_idx);
            XYZ_array_for_plane_fit = [X_array_for_plane_fit Y_array_for_plane_fit Z_array_for_plane_fit];

            % save the rest of the data for plotting if wanted
            X_array_not_for_plane_fit = X_array(Z_array_not_for_plane_fit_idx);
            Y_array_not_for_plane_fit = Y_array(Z_array_not_for_plane_fit_idx);
            Z_array_not_for_plane_fit = Z_array(Z_array_not_for_plane_fit_idx);
            XYZ_array_not_for_plane_fit = [X_array_not_for_plane_fit Y_array_not_for_plane_fit Z_array_not_for_plane_fit];
        end
        

        
        function C = planefit(x,y,z)
        % A function to fit x,y,z, data to a plane in 3D space.
        %
        % z = x * C(1) + y*C(2) + C(3);
        %
        % Example:
        % --------
        % x = -10:10;
        % y = -10:10;
        % [xx yy] = meshgrid(x,y);
        % zz = C(1)*xx+C(2)*yy + C(3) + 2*randn(size(xx));
        % plot3(xx(:),yy(:),zz(:),'.')
        % C = planefit(xx(:),yy(:),zz(:));
        % zzft = C(1)*xx+C(2)*yy + C(3);
        % hold on;
        % surf(xx,yy,zzft,'edgecolor','none')
        % grid on
        %
        % Val Schmidt
        % Center for Coastal and Ocean Mapping
        % University of New Hampshire
        % copywrite 2012

        xx = x(:);
        yy = y(:);
        zz = z(:);
        N = length(xx);
        O = ones(N,1);

        C = [xx yy O]\zz;
        
        end
        
        %%%=== PlaneFit_XYZarray ===%%%

        % This function takes in an Nx3 array (i.e., XYZ coordinates), and fits a
        % plane to the data (output as a matrix the same size as 'matrix', from
        % which the XYZarray would have been created).


        function [plane] = PlaneFit_XYZarray(matrix, XYZ_array)

            X_array = XYZ_array(:,1);
            Y_array = XYZ_array(:,2);
            Z_array = XYZ_array(:,3);

            C = ImageFuncs.planefit(X_array, Y_array, Z_array);

            a = C(1);
            b = C(2);
            c = C(3);

            plane = zeros(size(matrix));

            [row,col]=size(matrix);

            for i=1:col
                for j=1:row
                    plane(j,i) = (i*a) + (j*b) + c;
                end
            end

        end
        
        
        %%%=== Plane_Subtraction_Images_in_Cell ===%%%

        % Perform a 1st order plane background subtraction from each image in a
        % cell array.

        function [height_data_flat] = Plane_Subtraction_Images_in_Cell(height_data, BinWidth_nm, Plane_fit_mask_1, greater_than_1)

        height_data_flat = cell(size(height_data));

        for i = 1:length(height_data)

            height_mat = height_data{i};

            % transform to XYZ array for plane fitting
            [XYZ_array_1] = ImageFuncs.Matrix_to_Nx3array(height_mat);

            % create new XYZ array of only bottom 50% of data
            [XYZ_array_for_plane_fit_1, ~] = ImageFuncs.XYZarray_indexed_by_percentage_height(XYZ_array_1, BinWidth_nm, Plane_fit_mask_1, greater_than_1);

            % find the plane of the indexed bottom 50% of height data
            [plane_1] = ImageFuncs.PlaneFit_XYZarray(height_mat, XYZ_array_for_plane_fit_1);
            % subtract plane from data. This operation also brings the
            % background = 0 nm.
            height_mat_1_flat = height_mat - plane_1;

            % save background flattened data into new cell array: height_data_flat
            height_data_flat{i} = height_mat_1_flat;



        end

        end
        
        
        %%%=== InputCentres_Copper ===%%%
        
        % This function is used to assign the coordinates for the central axes of
        % rotation for the NPCs. The input argument is the AFM image in greyscale.
        % A figure is then produced showing the image, and the user must manually
        % click on the picture to assign the central axes of rotation for each NPC.
        % The output is an Nx2 array of the coordinates, to the nearest pixel, of
        % the selected centres.

        function Centres = InputCentres_Copper(Image)

        % for saturating colour scale, first take histogram of height data, binning
        % into 1nm bins

        heightdata=Image;
        heightdata_array = heightdata(:);
        histedges = [floor(min(heightdata_array)):ceil(max(heightdata_array))];
        hist_counts = histcounts(heightdata_array, histedges);
        hist_counts_sum = sum(hist_counts(:));

        % create an array with aggregate counts based on bins
        aggregate_hist_counts = zeros(size(hist_counts));
        aggregate_hist_counts(1) = hist_counts(1);
        for i=1:length(hist_counts)-1
            aggregate_hist_counts(i+1) = aggregate_hist_counts(i) + hist_counts(i+1);
        end

        % find bin which contains 98% of aggregate counts (or closest to)
        hist_counts_sum_ninety_five = hist_counts_sum*0.98;
        aggregate_hist_counts_abs = abs(aggregate_hist_counts - round(hist_counts_sum_ninety_five));
        cut_off_idx = aggregate_hist_counts_abs == min(aggregate_hist_counts_abs);

        % find bin edges (in nm) for this 98% of aggregate count, and define the
        % colour limits in nm
        colour_lim_floor = histedges(1);
        colour_lim_ceil = histedges(cut_off_idx);
        clims = [colour_lim_floor colour_lim_ceil];

        figure();
        imagesc(Image, clims);

        colormap(copper);
        xlabel('Select centres and press enter when finished')
        set (gcf, 'WindowButtonMotionFcn', @mouseMove);
        [c, r] = ginput;
        r = round(r);
        c = round(c);
        Centres = [r c];

        end
        
        
        
        
        
        
        %%%=== Drift_Correction_Height_Data ===%%%
        
        % Drift-correct the kymographs.

        function [Height_data_shifted, x_axis_minus_10nm_idx, x_axis_plus_10nm_idx] = Drift_Correction_Height_Data(Height_data, x_axis_nm, image_width, drift_shift_window_nm)

        % find index of 10nm
        x_axis_10nm_abs = abs(x_axis_nm - 10);
        x_axis_10nm_idx = find(x_axis_10nm_abs == min(x_axis_10nm_abs),1);

        % take the data in the first 10nm, find the mean, and subtract from thedata
        % set
        z_data_in_first_10nm = Height_data(:,1:x_axis_10nm_idx);
        z_data_in_first_10nm_mean = mean(z_data_in_first_10nm(:));

        % Height_data = Height_data - z_data_in_first_10nm_mean;

        % also take average of final 10 nm. This is for later.

        % find index of end-10nm
        x_axis_end_10nm_abs = abs(x_axis_nm - (image_width-10));
        x_axis_end_10nm_idx = find(x_axis_end_10nm_abs == min(x_axis_end_10nm_abs),1);

        % take the data in the first 10nm, find the mean, and subtract from thedata
        % set
        z_data_in_end_10nm = Height_data(:,x_axis_end_10nm_idx:end);
        z_data_in_end_10nm_mean = mean(z_data_in_end_10nm(:));


        %% at this point, have the slope corrected concatonated data, and the background set to 0nm. 

        % Now want to correct XY drift

        % will do this using cross-correlation averaging in a defined space around
        % where the AFM tip hits the side of the NPC/DNA origami. For this will
        % need a template, which will be the average height. Therefore, first find
        % the average height.

        [row_t, col_z] = size(Height_data);
        z_mean         = zeros(1, col_z);

        for i = 1:col_z

            pixel_t_z = Height_data(:,i);
            z_mean(1,i) = mean(pixel_t_z);   

        end

        %% define a col, and find the columns (+/-) drift_window_half_nm

        centres = ImageFuncs.InputCentres_Copper(Height_data);

        col_def_idx_array = centres(:,2)';

        drift_window_half_nm = round(drift_shift_window_nm)/2;

        col_def_nm_array = zeros(size(col_def_idx_array));
        col_minus_10_nm_array = zeros(size(col_def_idx_array));
        col_plus_10_nm_array  = zeros(size(col_def_idx_array));
        x_axis_minus_10nm_idx_array = zeros(size(col_def_idx_array));
        x_axis_plus_10nm_idx_array = zeros(size(col_def_idx_array));

        for i = 1:length(centres(:,1))
            col_def_idx = col_def_idx_array(i);

            % define desired range in nm
            col_def_nm = x_axis_nm(col_def_idx);
        %     col_minus_10_nm = col_def_nm - 10;
        %     col_plus_10_nm  = col_def_nm + 10;
            col_minus_10_nm = col_def_nm - drift_window_half_nm;
            col_plus_10_nm  = col_def_nm + drift_window_half_nm;

            col_def_nm_array(i) = col_def_nm;
            col_minus_10_nm_array(i) = col_minus_10_nm;
            col_plus_10_nm_array(i)  = col_plus_10_nm;

            % now find their index values
            % for minus 10nm
            x_axis_minus_10nm_abs = abs(x_axis_nm - col_minus_10_nm);
            x_axis_minus_10nm_idx = find(x_axis_minus_10nm_abs == min(x_axis_minus_10nm_abs),1);
            x_axis_minus_10nm_idx_array(i) = x_axis_minus_10nm_idx;
            % and for plus 10 nm
            x_axis_plus_10nm_abs = abs(x_axis_nm - col_plus_10_nm);
            x_axis_plus_10nm_idx = find(x_axis_plus_10nm_abs == min(x_axis_plus_10nm_abs),1);
            x_axis_plus_10nm_idx_array(i) = x_axis_plus_10nm_idx;
        end

        %% now, have defined range over which to search. Find lowest SAD scores compared to average

        z_mean_def = z_mean(1, x_axis_minus_10nm_idx_array(1):x_axis_plus_10nm_idx_array(1));

        if length(col_def_idx_array)>1
            for i = 1:length(col_def_idx_array)-1

                z_mean_def =  horzcat(z_mean_def, z_mean(1, x_axis_minus_10nm_idx_array(i+1):x_axis_plus_10nm_idx_array(i+1)));


            end
        end

        %%
        % at this point, will work in idx values (rather than nm) until minimum SAD
        % scores found and reconstruction begins
        Numb_idx_half_range = round((x_axis_plus_10nm_idx_array(1) - x_axis_minus_10nm_idx_array(1))/2)+1;

        Distance_idx_to_begin_img = x_axis_minus_10nm_idx_array(1);
        Distance_idx_to_end_img = col_z - x_axis_plus_10nm_idx_array(end);

        % this if statement ensures that no error can occur by searching too close
        % to the edges of the image.
        while Numb_idx_half_range>=Distance_idx_to_begin_img || Numb_idx_half_range>=Distance_idx_to_end_img
            Numb_idx_half_range = Numb_idx_half_range-1;
        end

        %% Shift data and save into cell arrays

        % this section pre-allocates cell arrays for shifting the data to the left,
        % and to the right. It then shifts the data in both directions and saves
        % the new arrays. These can then be compared against the template for their
        % SAD scores, and subsequent drift correction.
        Shifted_left_cell_pages = cell(row_t, Numb_idx_half_range, length(col_def_idx_array));
        Shifted_right_cell_pages = cell(row_t, Numb_idx_half_range, length(col_def_idx_array));

        for m = 1:length(col_def_idx_array)

            for n = 1:row_t

                height_data_array = Height_data(n,:);

                for i = 1:Numb_idx_half_range

                    % first one stored is not shifted (original), last one stored (end) is
                    % most shifted. Therefore, SAD_idx 1 = no shift. SAD_idx 3 = shift 2.
                    shift_left_array = height_data_array(:, x_axis_minus_10nm_idx_array(m)+(i-1):x_axis_plus_10nm_idx_array(m)+(i-1));
                    Shifted_left_cell_pages{n,i,m} = shift_left_array;

                    % again, first one stored is not shifted (original), last one stored (end) is
                    % most shifted. Therefore, SAD_idx 1 = no shift. SAD_idx 3 = shift 2.
                    shift_right_array = height_data_array(:, x_axis_minus_10nm_idx_array(m)-(i-1):x_axis_plus_10nm_idx_array(m)-(i-1));
                    Shifted_right_cell_pages{n,i,m} = shift_right_array;

                end
            end
        end

        %% Now have created the shifted arrays. For each template for the SAD, there
        % is a page in the cell array. Now must go through the pages and horzcat
        % the arrays to make them into one longer array

        Shifted_left  = Shifted_left_cell_pages(:,:,1);
        Shifted_right = Shifted_right_cell_pages(:,:,1);

        for n = 1:row_t

            for i = 1:Numb_idx_half_range

                for m = 1:length(col_def_idx_array)-1

                    Shift_left_temp_array = Shifted_left_cell_pages{n,i,m+1};
                    Shifted_left{n,i} = horzcat(Shifted_left{n,i},Shift_left_temp_array);

                    Shift_right_temp_array = Shifted_right_cell_pages{n,i,m+1};
                    Shifted_right{n,i} = horzcat(Shifted_right{n,i},Shift_right_temp_array);


                end

            end

        end

        %% Find the SAD scores.

        Shifted_left_SADs = zeros(size(Shifted_left));
        Shifted_right_SADs = zeros(size(Shifted_right));

        for n = 1:row_t

            for i = 1:length(Shifted_left(n,:))

                % pull out the shift height array
                height_array_left_shifted = Shifted_left{n, i};
                % take the absolute difference from the averaged template
                height_array_left_shifted_abs_diff = abs(height_array_left_shifted - z_mean_def);
                % sum the score
                SAD_left_score = sum(height_array_left_shifted_abs_diff(:));
                % store score
                Shifted_left_SADs(n,i) = SAD_left_score;

                % repeat operation for right
                height_array_right_shifted = Shifted_right{n, i};
                height_array_right_shifted_abs_diff = abs(height_array_right_shifted - z_mean_def);
                SAD_right_score = sum(height_array_right_shifted_abs_diff(:));
                Shifted_right_SADs(n,i) = SAD_right_score;

            end


        end

        %% Find position of min SAD scores for each row, left- and right-shifted

        Shifted_left_SAD_minScore_array_idx = zeros(row_t, 1);
        Shifted_left_SAD_minScore_array_absolute_score = zeros(row_t, 1);

        Shifted_right_SAD_minScore_array_idx = zeros(row_t, 1);
        Shifted_right_SAD_minScore_array_absolute_score = zeros(row_t, 1);

        for n = 1:row_t

            % pull out SAD score array
            Shifted_left_SAD_array = Shifted_left_SADs(n,:);
            % find index and absolute values of minimum score
            Shifted_left_SAD_minScore_idx = find(Shifted_left_SAD_array == min(Shifted_left_SAD_array),1);
            Shifted_left_SAD_minScore = Shifted_left_SAD_array(1, Shifted_left_SAD_minScore_idx);
            % save into arrays
            Shifted_left_SAD_minScore_array_idx(n) = Shifted_left_SAD_minScore_idx;
            Shifted_left_SAD_minScore_array_absolute_score(n) = Shifted_left_SAD_minScore;

            % repeat for right
            Shifted_right_SAD_array = Shifted_right_SADs(n,:);
            Shifted_right_SAD_minScore_idx = find(Shifted_right_SAD_array == min(Shifted_right_SAD_array),1);
            Shifted_right_SAD_minScore = Shifted_right_SAD_array(1, Shifted_right_SAD_minScore_idx);
            Shifted_right_SAD_minScore_array_idx(n) = Shifted_right_SAD_minScore_idx;
            Shifted_right_SAD_minScore_array_absolute_score(n) = Shifted_right_SAD_minScore;

        end

        %% Need to find out if the min score is in the left- or right-shifted for each

        Shift_left = zeros(row_t,1);
        Shift_right = zeros(row_t,1);
        Shift_idx = zeros(row_t,1);

        for n = 1:row_t

            % pull out the left and right min idx values
           min_score_left_idx  = Shifted_left_SAD_minScore_array_idx(n);
           min_score_right_idx = Shifted_right_SAD_minScore_array_idx(n);

           % pull out the left and right min values
           min_score_left  = Shifted_left_SAD_minScore_array_absolute_score(n);
           min_score_right = Shifted_right_SAD_minScore_array_absolute_score(n);

           if min_score_left == min_score_right % if the scores are the same, must position 1 (no shift), so say Shift Left and take its idx 
               Shift_left(n) = 1;
               Shift_right(n) = 0;
               Shift_idx(n) = min_score_left_idx;
           elseif min_score_left<min_score_right % if left less than right, take its idx and note down a 1 for shift left
               Shift_left(n) = 1;
               Shift_right(n) = 0;
               Shift_idx(n) = min_score_left_idx;
           elseif min_score_left>min_score_right % if right less than left, take its idx and give shift right a 1
               Shift_left(n) = 0;
               Shift_right(n)= 1;
               Shift_idx(n) = min_score_right_idx;
           end



        end

        %% Now know which direction to shift and how much to shift by, so can now do the drift correction

        % Remember, need to shidt by idx-1, as the first idx is unshifted.

        Height_data_shifted = zeros(size(Height_data));

        for n = 1:row_t


            Height_data_array = Height_data(n,:);

            if Shift_idx(n) == 1 % if idx value is 1, that means want unshifted data, so take original and skip to next iteration
                Height_data_shifted(n,:) = Height_data_array;
                continue
            end

            Height_data_shifted_array = zeros(size(Height_data(n,:)));

            if Shift_left(n) == 1 % if shifting to the left

                idx = Shift_idx(n);

                % take the height data from idx to end, i.e., if idx=2, this means
                % shift 1 to left. Therefore take from original data, starting at
                % index position 2 and going to end.
                Height_data_shifted_array(1:end-(idx-1)) = Height_data_array(idx:end);

                % then, in for loop, add noise back into the remaining positions.
                % For left-shifted, the gaps will be at the end. This is a random
                % number between -1 and 1, added to the average nm value for the 
                % final 10nm of the data set
                for i = 1:idx-1

                    g = -1 + (1+1)*rand(1,1); % generate a random number between -1 and 1
                    r = g + z_data_in_end_10nm_mean;

                    Height_data_shifted_array(end-(idx-i)) = r;

                end

                % save shifted data into array
                Height_data_shifted(n,:) = Height_data_shifted_array;

            elseif Shift_right(n) == 1 % else, if shifting to the right

                idx = Shift_idx(n);

                % take the height data from 1 to end-(idx-1), i.e., if idx=2, this means
                % shift 1 to right. Therefore take from original data (ignoring the
                % final point) and put into new array, starting at position 2
                % (i.e., ignoring the first point).
                Height_data_shifted_array(idx:end) = Height_data_array(1:end-(idx-1));

                % then, in for loop, add noise back into the remaining positions - 
                % for right-shifted, the gaps will be at the beginning.
                % This is a random number between -1 and 1, added to the average nm
                % value for the final 10nm of the data set
                for i = 1:idx-1

                    g = -1 + (1+1)*rand(1,1); % generate a random number between -1 and 1
                    r = g + z_data_in_first_10nm_mean;

                    Height_data_shifted_array(i) = r;

                end

                % save shifted data into array
                Height_data_shifted(n,:) = Height_data_shifted_array;

            end

        end



        end
        
        
         
        
    end %<-- end methods
    
end %<-- end classdef

