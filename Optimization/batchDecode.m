function decodedCandidates = batchDecode(candidates, decoding_func)
    
    decodedCandidates = cell(length(candidates), 1);
    for i = 1:length(candidates)
        decodedCandidates{i} = decoding_func(candidates{i});
    end

end