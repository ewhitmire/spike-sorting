function spikesort(dataset_path, detect_mode, cluster_mode)
% spikesort  Carries out spikesorting on a Quiroga dataset
%   dataset_path: Path to Quiroga simulated dataset file
%   detect_mode:    1 - amplitude threshold
%                   2 - slope threshold
%                   3 - wavelet
%   cluster_mode:   1 - k-means
%                   2 - SPC
%

    dataset = load(dataset_path);
    assignin('caller', 'dataset', dataset);
    data = dataset.data;
    
    % FILTER DATA
    filtered = filter_data(data, dataset.samplingInterval);
    assignin('caller', 'filtered', filtered);
    
    % SPIKE DETECTION
    if (detect_mode==1)
        [spikes, spike_times] = spike_detect(filtered, ...
            dataset.spike_times{1}, dataset.spike_class{2});
    elseif (detect_mode==2)
        [spikes, spike_times] = spike_detect_slope(filtered, ...
            dataset.spike_times{1}, dataset.spike_class{2});
    else
        [spikes, spike_times] = spike_detect_wavelets(filtered, ...
            dataset.spike_times{1}, dataset.spike_class{2});
    end
    
    assignin('caller', 'spikes', spikes);
    assignin('caller', 'spike_times', spike_times);
    
    % REALIGN SPIKE
    [int_spikes] = spike_interpolate(spikes);
    assignin('caller', 'int_spikes', int_spikes);
    
    % FEATURE EXTRACTION
    spike_features = feature_extract(int_spikes);
    assignin('caller', 'spike_features', spike_features);
    
    % CLUSTERING
    if (cluster_mode==1)
        spike_classes = spike_cluster(spike_features);
    else
        spike_classes = spike_cluster_spc(spike_features);
    end  
    assignin('caller', 'spike_classes', spike_classes);
    
    % EVALUATE
    evaluate(spike_classes, spike_times, dataset)

end