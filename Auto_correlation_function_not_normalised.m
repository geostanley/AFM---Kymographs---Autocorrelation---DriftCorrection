%%%=== Auto_correlation_function ===%%%

% Calcualted the auto-correlation of each pixel as a function of time-lag.

% This is defined as R_(tau) = [(X_(t)-mu) * (X_(t+tau)-mu)]/V
% Where R is the autocorrelation factor, and is a function of time-lag,
% tau; X_(t) is the time-dependent discrete function, and in this case, is
% z: height of a pixel with time; mu is the average height of that pixel
% over that time-frame.

function [R_matrix, l_range_t] = Auto_correlation_function_not_normalised(Height_data_analyse, Samples_per_line, Line_Rate_Hz, tau_max_secs, autocorrelation_over_defined_secs, dt_perLine)

% Autocorrelation function
pix_nos = 1:Samples_per_line;
Auto_correlation_R = cell(size(pix_nos));

sec_per_line = 1/Line_Rate_Hz;
lines_per_defined_secs = round(tau_max_secs/sec_per_line);

%% calculate auto-correlation factor

for n = 1:length(pix_nos) %<-- for each line of height data

    z = Height_data_analyse(:,pix_nos(n)); %<-- take array of height data
    mu = mean(z); %<-- take mean

    if autocorrelation_over_defined_secs == 1
        l_mx     = lines_per_defined_secs; %<-- max calculated tau value (index value)
        l_range  = (1:l_mx); %<-- range of tau values to be calculated (index values)
    else
        l_mx     = round(length(z))/2;
        l_range  = (1:l_mx);
    end

    Rnl = cell(1,length(l_range)); %<-- pre-allocate

    for j = 1:length(l_range);
        for i = 1:length(z)-l_range(j);
            tau = l_range(j);
            Rnl{j}(i) = ((z(i)-mu)*(z(i+tau)-mu)); %<-- calculate all R values at a given time-lag and store in cell array
        end
    end

    R = zeros(1,length(Rnl)); %<-- pre-allocate
    
    for a = 1:length(Rnl);
        R_array = Rnl{a};
        R(a)    = sum(R_array)/length(Rnl{a}); %<-- calculate final R value at given time-lag
    end

    Auto_correlation_R{n} = R; %<-- store in cell array
    
end

%% Take Auto_correlation_R cell and turn into Matrix

R_matrix = zeros(length(Auto_correlation_R{1}), length(Auto_correlation_R));

for i = 1:length(Auto_correlation_R)
    R_vec         = Auto_correlation_R{i}(:);
    R_matrix(:,i) = R_vec;  
end

l_range_t = linspace(dt_perLine, dt_perLine*l_range(end), l_range(end));

end