function [spikes, spike_times] = spike_detect(data, spike_times_gt, spike_overlap_gt)
    slopes = diff(data);
%     figure; hold on; plot(slopes(1:length(slopes)));
    threshLine = ones(1,length(slopes));
    sigma = median(abs(slopes)/.6745);
    threshold = 4.8 * sigma;
    threshLine = threshLine * threshold;
%     plot(threshLine, 'Color', [1 0 0]);
    
%     figure; hold on;
%     
%     hData = plot([1:length(data)]/24, data, 'Color', [0 0 .8], 'LineWidth', .1);
%     
%     hTitle = title('Simulated Spike Data');
%     hXLabel = xlabel('Time (ms)');
%     hYLabel = ylabel('Voltage (uV)');
%     
%     
% 
%     set( gca                       , ...
%         'FontName'   , 'Helvetica' );
%     set([hTitle, hXLabel, hYLabel], ...
%         'FontName'   , 'AvantGarde');
%     set([hXLabel, hYLabel]  , ...
%         'FontSize'   , 10          );
%     set( hTitle                    , ...
%         'FontSize'   , 12          , ...
%         'FontWeight' , 'bold'      );
% 
%     set(gca, ...
%       'Box'         , 'off'     , ...
%       'TickDir'     , 'out'     , ...
%       'TickLength'  , [.02 .02] , ...
%       'XMinorTick'  , 'on'      , ...
%       'YMinorTick'  , 'on'      , ...
%       'YGrid'       , 'on'      , ...
%       'XColor'      , [.3 .3 .3], ...
%       'YColor'      , [.3 .3 .3], ...
%       'YTick'       , -8:2:8, ...
%       'LineWidth'   , 1         );
%     set(gcf, 'PaperPositionMode', 'auto');

            
            
    spikes_indices = abs(slopes) > threshold;
    spike_pruned_indices = pruneSpikes(data, spikes_indices, threshold);
    spike_times = spike_pruned_indices;
    
    spikes = [];
    N = 64;
    data = [zeros(1, 20), data, zeros(1, 44)];
    for i = 1:N
        spikes = [spikes; data(spike_pruned_indices+i+20)];
    end

    
%     figure; hold on;
%     plot(spike_times_gt, 'r');
%     plot(spike_times);
    
    % figure; hold on;
    % plot(spike_times_gt(find(spike_overlap_gt==0)));
    
    
end

function indices = pruneSpikes(data, spikes, threshold)
    i_all = find(spikes);
    i_next_limit = 0;
    indices = [];
    for i = 1:length(i_all)
        index = i_all(i)-5;
        if (index > i_next_limit)
            window = data(index:min(length(data), index+44));
            %if (max(window) > threshold) % this may be unnecessary
                max_index = find(window == max(window), 1)+index - 1;
                indices = [indices, max_index-20];
                i_next_limit = max_index + 37;%find(window < threshold/2, 1)+index - 1;
            %end
        end
    end
end