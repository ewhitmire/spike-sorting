function evaluate(spike_classes, spike_times, dataset)
    gt_times = dataset.spike_times{1};
    gt_classes = dataset.spike_class{1};
    gt_nonoverlapped_index = find(dataset.spike_class{2}==0);
    gt_class_non = gt_classes(gt_nonoverlapped_index);
    gt_times_non = gt_times(gt_nonoverlapped_index);
    gt_classes_adjusted = [];
    spike_classes_adjusted = [];
    
    gt_classes_adjusted_non = [];
    spike_classes_adjusted_non = [];
    % spike detector
    correct = 0;
    missed = 0;
    false_positive = 0;
    
    
    correct_non = 0;
    missed_non = 0;
    false_positive_non = 0;
    
    fps = [];
    
    found_gt_spike_indices = [];
    found_gt_spike_indices_non = [];
    
    TIME_THRESH = 15; % units
    
    
    for i = 1:length(spike_classes)
       % for each of our detected spikes
       start_time = spike_times(i);
       gt_start_time_index_non = find(abs(gt_times_non - start_time) == min(abs(gt_times_non - start_time)), 1);
       gt_start_time_index = find(abs(gt_times - start_time) == min(abs(gt_times - start_time)), 1);
       gt_nearest_start_time_non = gt_times_non(gt_start_time_index_non);
       gt_nearest_start_time = gt_times(gt_start_time_index);
       
       if (abs(gt_nearest_start_time_non - start_time) < TIME_THRESH && sum(ismember(found_gt_spike_indices_non, gt_start_time_index_non))==0)
           correct_non = correct_non + 1;
           gt_classes_adjusted_non = [gt_classes_adjusted_non, gt_class_non(gt_start_time_index_non)];
           spike_classes_adjusted_non = [spike_classes_adjusted_non, spike_classes(i)];
           found_gt_spike_indices_non = [found_gt_spike_indices_non, gt_start_time_index_non];
       else
           false_positive_non = false_positive_non + 1;
       end
       
       if (abs(gt_nearest_start_time - start_time) < TIME_THRESH && sum(ismember(found_gt_spike_indices, gt_start_time_index))==0)
           correct = correct + 1;
           gt_classes_adjusted = [gt_classes_adjusted, gt_classes(gt_start_time_index)];
           spike_classes_adjusted = [spike_classes_adjusted, spike_classes(i)];
           
           found_gt_spike_indices = [found_gt_spike_indices, gt_start_time_index];
       else
           false_positive = false_positive + 1;
           fps = [fps, gt_times(gt_start_time_index)];
       end
       


           
    end
    missed_non = length(gt_times_non) - correct_non;
    missed = length(gt_times) - correct;
    
    total = length(gt_times);
    display(length(gt_times_non));
    display(false_positive);
    
    display(correct_non);
    display(missed_non);
    C = confusionmat(spike_classes_adjusted_non, gt_classes_adjusted_non);
    display(C);
    disp 'done!';
end