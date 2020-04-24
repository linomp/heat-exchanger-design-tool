function decoded = decodeChromosomeDefault(ind)
    % function to decode a chromosome and return as set of parameters
    % with names n' all (this is problem dependent, hence kept as a 
    % separate function not belonging to any class)
    
    paramNames = 1:length(ind.chrom); %'a':'z';
    
    %fprintf('\nBest Individual Fitness: %f \n', ind.fitness);
    fprintf('\n-- BEST PARAMETERS --');
    for i = 1:length(ind.chrom)
        fprintf('\n * parameter %i: %f', paramNames(i), ind.chrom(i));
    end
    fprintf('\n\n');
    decoded = ind;
end