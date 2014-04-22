function filtered=filter_data(data, sampling_interval)
    close all;
    Fs = (1000 / sampling_interval);  % Sampling Frequency


    % Calculate the coefficients using the FIRLS function.
    Hd = BartlettHanning;
    
   
    if (isprop(Hd, 'sosMatrix'))
        filtered = filtfilt(Hd.sosMatrix, 1, data);
    elseif (isprop(Hd, 'Numerator'))
        filtered = filtfilt(Hd.Numerator, 1, data);
    else
        display('error');
    end

end