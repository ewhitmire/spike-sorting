% relies on waveletdetect.m

function [spikes, spike_max] = spike_detect_wavelets(data, spike_times_gt, ...
    spike_overlap_gt)
    
    SFr = 24; % sample frequency in kHz
    Wid = [.35 1.25]; % min and max durations of transients, in msec
    Ns = 8; % Number of scales to be used
    option = 'c'; % conservative--returns nothing if no spikes are detected
                  % using hard threshold
    L = -.1; % value that sets tendency for over- and under-detection of
              % spikes
    wname = 'bior1.5'; % name of wavelet family
    PltFlg = 0; % 1 --> generate figures, 0 --> do not
    CmtFlg = 1; % 1 --> display comments, 0 --> do not
    
    spikes_indices = detect_spikes_wavelet(...
    data, SFr, Wid, Ns, option, L, wname, PltFlg, CmtFlg);
    
    spike_max = [];
    for i=1:length(spikes_indices)
        index = spikes_indices(i);
        window = 20;
        max_index = find(data(index:index+window)==max(data(index:index+window))) - 1 + index;
        spike_max = [spike_max, max_index-20];
    end

    spikes = [];
    N = 64;
    data = [zeros(1, 20), data, zeros(1, 44)];
    for i = 1:N
        spikes = [spikes; data(spike_max+i+20)];

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
