function spikes_int = spike_interpolate(spikes)
    spikes_int = [];
    decimation = 4;
    x = 1:64;
    xx = 0:1/decimation:66;
    for i = 1:size(spikes, 2)
        spike = spikes(:, i);
        spike_int = spline(x, spike', xx);
        peak_range = (19 * decimation+1:21 * decimation-1);
        peak_index = find(spike_int(peak_range)==max(spike_int(peak_range))) - 1 + peak_range(1);
        new_spike = [spike_int(peak_index-decimation*19:decimation:peak_index-decimation), spike_int(peak_index:decimation:peak_index+decimation*44)];
        spikes_int = [spikes_int, new_spike'];
    end
end