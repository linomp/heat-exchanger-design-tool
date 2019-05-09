function top_candidates = getCandidatesFromPool(pool, n)
    
    ftns = zeros(length(pool), 1);
    for i = 1:length(pool)
        ftns(i) = pool{i}.fitness;
    end
    
    [~, idx] = sort(ftns, 'descend');
    orderedPool = pool(idx); 
    
    set = {orderedPool{1}};
    for i = 1:length(orderedPool) 
        if size(set, 2) >= n
            break 
        end
        dup = 0;
        for j = 1:length(set) 
            if prod(round(set{j}.chrom, 3) == round(orderedPool{i}.chrom, 3))
            %if eq(set{j}, orderedPool{i})
                dup = 1;
                break 
            end
        end
        if ~dup
            set{end+1} = orderedPool{i};
        end 
    end
    
    if length(set) > n
        top_candidates = set(:, 1:n);
    else
        top_candidates = set;
    end
    
end