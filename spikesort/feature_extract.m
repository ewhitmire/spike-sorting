function features=feature_extract(data)
    Cs = [];
    for i=1:size(data, 2)
        x = data(:,i);
        C = wavedec(x, 4, 'haar');
        Cs = [Cs, C];
    end
    
    
    Ks = []
    for i=1:size(Cs,1)
       feature = Cs(i, :);
       [h, p, ksstat, cv] = kstest(feature);
       Ks = [Ks; ksstat];
    end
    [Ks_sorted, index] = sort(Ks, 1, 'descend');
    N_FEATURES = 10;
    top_indices = index(1:N_FEATURES);
    
    features = Cs(top_indices, :);
    
end