function [spikes, spike_times] = spike_detect(data, spike_times_gt, spike_overlap_gt)
    sigma = median(abs(data)/.6745);
    threshold = 4 * sigma;
    
    spikes_indices = (data) > threshold;
    spike_pruned_indices = pruneSpikes(data, spikes_indices);
    spike_times = spike_pruned_indices;
    
    spikes = [];
    N = 64;
    data = [zeros(1, 20), data, zeros(1, 44)];
    for i = 1:N
        spikes = [spikes; data(spike_pruned_indices+i+20)];
    end

end

function indices = pruneSpikes(data, spikes)
    i_all = find(spikes);
    i_next_limit = 0;
    indices = [];
    for i = 1:length(i_all)
        index = i_all(i);
        if (index > i_next_limit)
            window = data(index:min(length(data), index+44));
            %if (max(window) > threshold) % this may be unnecessary
                max_index = find(window == max(window), 1)+index - 1;
                indices = [indices, max_index-20];
                i_next_limit = max_index + 20;%find(window < threshold/2, 1)+index - 1;
            %end
        end
    end
end