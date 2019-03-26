classdef MyPlots
    % MyPlots: Contains my various plotting functions.
       
    methods(Static)
        
        function Height_Plane_Scatter(XYZ_remaining_data, XYZ_array_for_plane_fit, plane)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure();
        scatter3(XYZ_remaining_data(:,1),XYZ_remaining_data(:,2),XYZ_remaining_data(:,3),'b.');
        hold on
        scatter3(XYZ_array_for_plane_fit(:,1), XYZ_array_for_plane_fit(:,2), XYZ_array_for_plane_fit(:,3),'r.');
        colormap(gray)
        surf(plane, 'edgecolor', 'none','facecolor','green','facealpha',0.5)
        
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        title('Height data with plane fitting', 'FontSize', 13)
        
        end
        
        function SurfaceKymograph(clims_z, x_axis_nm, time, Height_data)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure();
        surf(x_axis_nm, time, Height_data, 'edgecolor', 'none')
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        caxis(clims_z)
        colorbar
        colormap(gray)

        xlabel('Distance (nm)', 'FontSize', 13)
        ylabel('Time (s)', 'FontSize', 13)
        
        end
        
        
        function KymographBeforeAfterWithStdSubplots(clims_z, x_axis_nm, time, Height_data_before_shift, z_unshifted_std, Height_data_analyse, z_shifted_std, Unique_save_name, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
            
        % make some nice colours for plotting
        C  = linspecer(4);
        c1 = C(1,:);
        
        std_mx = max(z_shifted_std);

        StdFigs = figure();

        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        subplot(411)
        imagesc(x_axis_nm, time, Height_data_before_shift)
        colormap(gray)
        caxis(clims_z)
        set(gca, 'FontSize', 13, 'LineWidth', 2, 'YDir','normal')
        xlabel('Distance (nm)', 'FontSize', 13)
        ylabel('Slow-scan axis (time)', 'FontSize', 13)
        title('Original Height Data', 'FontSize', 13)

        subplot(412)
        plot(x_axis_nm, z_unshifted_std, 'Color', c1, 'LineWidth', 2)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        xlabel('Distance (nm)', 'FontSize', 13) 
        ylabel('Standard deviation', 'FontSize', 13)
        title('SD of original height data', 'FontSize', 13)
        ylim([0 ceil(std_mx)+1])

        subplot(413)
        imagesc(x_axis_nm, time, Height_data_analyse)
        colormap(gray)
        caxis(clims_z)
        set(gca, 'FontSize', 13, 'LineWidth', 2, 'YDir','normal')
        xlabel('Distance (nm)', 'FontSize', 13), 
        ylabel('Slow-scan axis (time)', 'FontSize', 13)
        title('Drift-corrected height Data', 'FontSize', 13)

        subplot(414)
        plot(x_axis_nm, z_shifted_std, 'Color', c1, 'LineWidth', 2)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        xlabel('Distance (nm)', 'FontSize', 13),
        ylabel('Standard deviation', 'FontSize', 13)
        title('SD of drift-corrected height data', 'FontSize', 13)
        ylim([0 ceil(std_mx)+1])

        set(StdFigs, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);
        
        fig_name = strcat(Unique_save_name, ' - Standard_Deviations', ' ', num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'jpeg');
        end
        
        end
        
        
        function KymographPlot(clims_z, x_axis_nm, time, Kymograph)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figglywiggly = figure();
        imagesc(x_axis_nm, time, Kymograph)
        colormap(gray)
        caxis(clims_z)
        set(gca, 'FontSize', 13, 'LineWidth', 2, 'YDir','normal')
        xlabel('Distance (nm)', 'FontSize', 13), 
        ylabel('Slow-scan axis (time)', 'FontSize', 13)
        title('Kymograph', 'FontSize', 13)

        set(figglywiggly, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);
        
        end
        
        function KymographPlotColourbar(clims_z, x_axis_nm, time, Kymograph, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure()
        imagesc(x_axis_nm, time, Kymograph) 
        title('Kymograph', 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        colormap(gray)
        caxis(clims_z)
        gc = colorbar;
        set(gc, 'LineWidth', 2)
        ylabel(gc,'Height (nm)', 'FontSize', 13)
        set(gca,'ydir','normal', 'FontSize', 13, 'LineWidth', 2);
        pbaspect([1 2 1])
        fig_name = strcat(Unique_save_name, ' - ', FigureSaveName, ' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'pdf');
        end
        
        end
        
        
        function SimpleLinePlot_Ratio121(ylims, x_axis_nm, z_shifted_std, Unique_save_name, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
            
        % make some nice colours for plotting
        C  = linspecer(4);
        c1 = C(1,:);
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(),
        plot(x_axis_nm, z_shifted_std, 'Color', c1, 'LineWidth', 2)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        xlabel('Distance (nm)', 'FontSize', 13),
        ylabel('Standard deviation', 'FontSize', 13)
        ylim(ylims)
        xlim([0 max(x_axis_nm)])

        pbaspect([1 2 1])

        fig_name = strcat(Unique_save_name, ' - Standard_Deviation_DriftCorrectedOnly', ' ', num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'pdf');
        end
        
        end
        
        
        
        function KymographSubplots(clims_z, x_axis_nm, time, Height_data_before_shift, Height_data_analyse, Unique_save_name, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        heightplots_fig = figure();

        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        subplot(211)
        imagesc(x_axis_nm, time, Height_data_before_shift), 
        title('Original height data', 'FontSize', 13)
        xlabel('Distance (nm)', 'FontSize', 13) 
        ylabel('Time (s)',  'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2, 'YDir','normal')
        colormap(gray)
        caxis(clims_z)
        colorbar

        subplot(212),
        imagesc(x_axis_nm, time, Height_data_analyse), 
        title('Drift-corrected height data', 'FontSize', 13)
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)',  'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2, 'YDir','normal')
        colormap(gray)
        caxis(clims_z)
        colorbar

        set(heightplots_fig, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);

        fig_name = strcat(Unique_save_name, ' - Drift correction', ' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'jpeg');
        end
        
        end
        
        
        function AutoCorrelationHeatmapColourbar(clims_R, colourmap, x_axis_nm, l_range_t, R_matrix, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(), imagesc(R_matrix, 'XData', [min(x_axis_nm) max(x_axis_nm)], 'YData', [min(l_range_t) max(l_range_t)])
        title('Auto-correlation factor, R({\tau})', 'FontSize', 13)
        xlabel('Fast-scan axis (nm)', 'FontSize', 13)
        ylabel('\tau (s)',  'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2, 'ydir', 'norm')
        colormap(colourmap)
        caxis(clims_R)
        g = colorbar;
        set(g, 'LineWidth', 2)
        ylabel(g,'R({\tau})', 'FontSize', 13)
        pbaspect([1 2 1])
        fig_name = strcat(Unique_save_name, ' - ', FigureSaveName, ' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'pdf');
        end
        
        end
        
        
        function[l_range_t_plot_axis_scaled, R_normalised_matrix_scaled] = KymographAutoCorrelationHeatmapOverlay(clims_z, clims_R, colourmap, x_axis_nm, time, Height_data_analyse, l_range_t, R_not_normalised_matrix, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        R_normalised_matrix_scaled = imresize(R_not_normalised_matrix, size(Height_data_analyse));
        l_range_t_plot_axis_scaled = linspace(min(l_range_t), max(l_range_t), length(time));
            
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(), 
        f1 = imagesc(x_axis_nm, time, Height_data_analyse);
        f1.Parent.YAxisLocation = 'right';
        caxis(clims_z);
        colormap(gray);
        xlabel('Fast-scan axis (nm)', 'FontSize', 13);
        ylabel('Time (s)', 'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        set(gca,'ydir','normal');
        ax1 = gca;
        freezeColors
        hold on
        f2 = imagesc(x_axis_nm, l_range_t_plot_axis_scaled, R_normalised_matrix_scaled, 'AlphaData', 0.5);
        ax2 = gca;
        colormap(colourmap)
        caxis(clims_R)
        gc = colorbar;
        set(gc, 'LineWidth', 2)
        ylabel(gc,'R({\tau})', 'FontSize', 13)
        f1.Parent.YAxisLocation = 'right';
        f2.Parent.YAxisLocation = 'left';
        pbaspect([1 2 1])
        fig_name = strcat(Unique_save_name, ' - ', FigureSaveName, ' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'pdf');
        end

        end
        
        
        
        function SubplotKymographAutoCorrelationHeatmap(clims_z, clims_R, colourmap, x_axis_nm, time, Height_data_analyse, l_range_t, R_matrix, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(), f1 = subplot(211);
        imagesc(x_axis_nm, time, Height_data_analyse) 
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        title('Kymograph', 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        colormap(f1, 'gray')
        caxis(clims_z)

        gc = colorbar;
        set(gc, 'lineWidth', 2)
        ylabel(gc,'Height (nm)', 'FontSize', 13)
        pbaspect([1 2 1])

        subplot(212),
        imagesc(R_matrix, 'XData', [min(x_axis_nm) max(x_axis_nm)], 'YData', [min(l_range_t) max(l_range_t)])
        title('Autocorrelation factor, R(\tau)', 'FontSize', 13)
        xlabel('Fast-scan axis (nm)', 'FontSize', 13)
        ylabel('\tau (s)',  'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        colormap(colourmap)
        caxis(clims_R)
        g = colorbar;
        ylabel(g,'R(\tau)', 'FontSize', 13)
        set(gca,'ydir','normal', 'linewidth', 2);
        set(g, 'linewidth', 2);
        pbaspect([1 2 1])

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);


        fig_name = strcat(Unique_save_name, ' - ', FigureSaveName, ' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'pdf');
        end

        
        end
        
        
        
        function SubplotKymographAutoCorrelationHeatmap_simplesave(clims_z, clims_R, halfwindow_nm, dt_plot, R_tau_max, colourmap, x_axis_nm, time, Height_data_analyse, l_range_t, R_matrix, figure_save_name, save_outputs)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(),
        f1 = subplot(211);
        imagesc(x_axis_nm, time, Height_data_analyse) 
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 1)
        title('Cropped kymograph', 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        colormap(f1, 'gray')
        caxis(clims_z)
        xlim([-halfwindow_nm halfwindow_nm])

        set(gca,'XTick', [-15 -10 -5 0 5 10 15])
        set(gca,'XTickLabel',{'-15', '', '', '0', '', '', '15'})

        gc = colorbar;
        set(gc, 'lineWidth', 1)
        ylabel(gc,'Height (nm)', 'FontSize', 13)
        pbaspect([1 2 1])

        subplot(212);
        imagesc(x_axis_nm, l_range_t, R_matrix)
        colormap(colourmap)
        gc = colorbar;
        set(gc, 'LineWidth', 1)
        ylabel(gc,'R(\tau) (nm^{2})', 'FontSize', 13)
        caxis(clims_R)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 1)
        title(strcat('R(\tau) (nm^{2}); total t = ', num2str(max(time)), 's'), 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        xlim([-halfwindow_nm halfwindow_nm])
        ylim([dt_plot R_tau_max])

        set(gca,'XTick', [-15 -10 -5 0 5 10 15])
        set(gca,'XTickLabel',{'-15', '', '', '0', '', '', '15'})

        pbaspect([1 2 1])

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);

        if save_outputs == 1
            saveas(gca, figure_save_name, 'pdf');
        end

        
        end
        
        
        
        
        

        
        
        
        
        
        function[log10_l_range_t] = AutoCorrelationHeatmapLogarithmicPlot(clims_R, colourmap, dt_perLine, x_axis_nm, l_range_t, R_not_normalised_matrix, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        log10_l_range_t          = zeros(size(l_range_t));
        log10_l_range_t(1:end-1) = l_range_t(2:end);
        log10_l_range_t(end)     = l_range_t(end)+dt_perLine;
        log10_l_range_t          = log10(log10_l_range_t);

        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure(), mesh(x_axis_nm, log10_l_range_t, R_not_normalised_matrix, 'FaceColor', 'interp')
        axis([min(x_axis_nm) max(x_axis_nm) min(log10_l_range_t) max(log10_l_range_t)])
        title('Loagarithmic plot of R against time-lag', 'FontSize', 13)
        xlabel('Fast-scan axis (nm)', 'FontSize', 13)
        ylabel('log(Time-lag (s))',  'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        colormap(colourmap)
        caxis(clims_R)
        g = colorbar;
        ylabel(g,'R', 'FontSize', 13)
        set(gca,'ydir','normal');
        set(g, 'LineWidth', 2)
        pbaspect([1 2 1])
        view(0, 90);

        fig_name = strcat(Unique_save_name, ' - ', FigureSaveName, ' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'jpeg');
        end
        
        
        end
        
        
        
        
        
        
        
        function SubPlotKymographAutoCorrelationHeatmapLogarithmic(clims_z, clims_R, colourmap, dt_perLine, x_axis_nm, time, Height_data_analyse, l_range_t, R_not_normalised_matrix, Unique_save_name, FigureSaveName, Save_figure_directory, Line_Rate_Save, channel, save_outputs)
        
        log10_l_range_t          = zeros(size(l_range_t));
        log10_l_range_t(1:end-1) = l_range_t(2:end);
        log10_l_range_t(end)     = l_range_t(end)+dt_perLine;
        log10_l_range_t          = log10(log10_l_range_t);
            
        figure(), 

        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        f1 = subplot(211);
        imagesc(x_axis_nm, time, Height_data_analyse) 
        title('Kymograph', 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        colormap(f1, gray)
        caxis(clims_z)
        gc = colorbar;
        ylabel(gc,'Height (nm)', 'FontSize', 13)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        set(gc, 'FontSize', 13, 'LineWidth', 2)

        subplot(212),
        mesh(x_axis_nm, log10_l_range_t, R_not_normalised_matrix, 'FaceColor', 'interp')
        axis([min(x_axis_nm) max(x_axis_nm) min(log10_l_range_t) max(log10_l_range_t)])
        title('Loagarithmic plot of R against time-lag')
        xlabel('Fast-scan axis (nm)', 'FontSize', 13)
        ylabel('log(Time-lag (s))',  'FontSize', 13)
        set(gca, 'FontSize', 13, 'LineWidth', 2)
        colormap(colourmap)
        caxis(clims_R)
        g = colorbar;
        ylabel(g,'R', 'FontSize', 13)
        set(gca,'ydir','normal');
        set(g, 'LineWidth', 2)
        view(0, 90);

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);

        fig_name = strcat(Unique_save_name, ' - ', FigureSaveName,' - ',num2str(Line_Rate_Save), ' - ', num2str(channel));
        if save_outputs == 1
            saveas(gca, fullfile(Save_figure_directory, fig_name), 'jpeg');
        end
        
        
        end
        
        
        
        
        function SubPlotKymographAutoCorrelationHeatmapLogarithmic_simplesave(clims_z, clims_R, colourmap, halfwindow_nm, dt_plot, R_tau_max, x_axis_nm, time, Height_data_analyse, l_range_t, R_mat, figure_save_name, save_outputs)
        
        figure(),
        f3 = subplot(211);
        imagesc(x_axis_nm, time, Height_data_analyse) 
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 1)
        title('Cropped kymograph', 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        colormap(f3, 'gray')
        caxis(clims_z)
        xlim([-halfwindow_nm halfwindow_nm])

        set(gca,'XTick', [-15 -10 -5 0 5 10 15])
        set(gca,'XTickLabel',{'-15', '', '', '0', '', '', '15'})

        gc = colorbar;
        set(gc, 'lineWidth', 1)
        ylabel(gc,'Height (nm)', 'FontSize', 13)
        pbaspect([1 2 1])

        subplot(212);
    %     mesh(x_axis_nm_cropped, l_range_t, R_norm_mat_cropped, 'FaceColor', 'interp')
        surf(x_axis_nm, l_range_t, R_mat, 'FaceColor', 'interp', 'Edgecolor', 'none')
        colormap(colourmap)
        gc = colorbar;
        set(gc, 'LineWidth', 1)
        ylabel(gc,'R(\tau) (nm^{2})', 'FontSize', 13)
        set(gca, 'YScale', 'log')
        caxis(clims_R)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 1)
        title(strcat('Logarithmic plot of R(\tau) (nm^{2}); total t = ', num2str(max(time)), 's'), 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('\tau (s)', 'FontSize', 13)
        xlim([-halfwindow_nm halfwindow_nm])
        ylim([dt_plot R_tau_max])

        set(gca,'XTick', [-15 -10 -5 0 5 10 15])
        set(gca,'XTickLabel',{'-15', '', '', '0', '', '', '15'})

        view(0, 90);

        pbaspect([1 2 1])

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);

        if save_outputs == 1
            saveas(gca, figure_save_name, 'tif');
        end
        
        
        end
        
        
        
        function PlotAveragedAutoCorrelationHeatmap_simplesave(clims_R, halfwindow_nm, dt_plot, R_tau_max, colourmap, x_axis_nm, time, l_range_t, R_matrix, figure_save_name, save_outputs)
        
        set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

        figure();
        imagesc(x_axis_nm, l_range_t, R_matrix)
        colormap(colourmap)
        gc = colorbar;
        set(gc, 'LineWidth', 1)
        ylabel(gc,'R(\tau) (nm^{2})', 'FontSize', 13)
        caxis(clims_R)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 1)
        title(strcat('R(\tau) (nm^{2}); total t = ', num2str(max(time)), 's'), 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('Time (s)', 'FontSize', 13)
        xlim([0 halfwindow_nm])
        ylim([dt_plot R_tau_max])

        set(gca,'XTick', [0 5 10 15])
        set(gca,'XTickLabel',{ '0', '5', '10', '15'})

        pbaspect([1 2 1])

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);

        if save_outputs == 1
            saveas(gca, figure_save_name, 'pdf');
        end

        
        end
        
        
        
        function PlotAveragedAutoCorrelationHeatmapLogarithmic_simplesave(clims_R, colourmap, halfwindow_nm, dt_plot, R_tau_max, time, x_axis_nm, l_range_t, R_mat, figure_save_name, save_outputs)
        
        figure();
        surf(x_axis_nm, l_range_t, R_mat, 'FaceColor', 'interp', 'Edgecolor', 'none')
        view(0, 90);
        colormap(colourmap)
        gc = colorbar;
        set(gc, 'LineWidth', 1)
        ylabel(gc,'R(\tau) (nm^{2})', 'FontSize', 13)
        set(gca, 'YScale', 'log')
        caxis(clims_R)
        set(gca,'ydir','normal');
        set(gca, 'FontSize', 13, 'LineWidth', 1)
        title(strcat('Logarithmic plot of R(\tau) (nm^{2}); total t = ', num2str(max(time)), 's'), 'FontSize', 13),
        xlabel('Fast-scan axis (nm)', 'FontSize', 13), 
        ylabel('\tau (s)', 'FontSize', 13)
        xlim([0 halfwindow_nm])
        ylim([dt_plot R_tau_max])
        
        set(gca,'XTick', [0 5 10 15])
        set(gca,'XTickLabel',{ '0', '5', '10', '15'})

        view(0, 90);
        pbaspect([1 2 1])

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.1, 0.3, 0.9]);

        if save_outputs == 1
            saveas(gca, figure_save_name, 'tif');
        end
        
        
        end
        
        
          
        
    end
    
end

