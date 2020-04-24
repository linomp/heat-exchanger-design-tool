function bestParameters = runOptimization(data)

problem = HBProblem();   
problem.objective_func = @HeatExchangerObjectiveFunc; 
% last 2 parameters: integers
problem.constraints_func = @(chrom) [chrom(1:3), round(chrom(4:end))];

% PARAMETERS UPPER AND LOWER BOUNDS
% Individual Structure: 
% Ds : Shell outer diameter
% do : tube outer diameter
% b : baffle spacing 
% n : number of passes
% Nt: number of tubes 
problem.LB = data.LB; 
problem.UB = data.UB;  
problem.bc = data;

problem.grid_dims = [2 2 2 2 2];

% Perform optimization
problem.n_candidates = 5;
[~, candidates, avgs, bests] = optimizationLauncher(problem);
[cost, extra] = HeatExchangerObjectiveFunc(candidates{1}.chrom, data);
          

raw = candidates{1}.chrom;
bestParameters.r_shell = raw(1)/2;
bestParameters.r_tube = raw(2)/2;
bestParameters.baffle = raw(3);
bestParameters.n_passes = raw(4);
bestParameters.n_tubes = raw(5);
bestParameters.U = extra.U*1000;
bestParameters.v_tube = extra.v_tube;
bestParameters.v_shell = extra.v_shell;
bestParameters.L = extra.L;
bestParameters.T_shell_out = extra.T_shell_out;
bestParameters.T_tube_out = extra.T_tube_out;
bestParameters.A = extra.A;
bestParameters.Q = extra.Q;
bestParameters.cost = cost;

% plot GA convergence
%{
figure(1)
hold
semilogy(1:length(bests(2:end)), bests(2:end), 'g-', 'linewidth', 3);
semilogy(1:length(avgs(2:end)), avgs(2:end), 'r--', 'linewidth', 2);
L = legend('Best Individual', 'Population Avg.');
set(L, 'location', 'best', 'fontSize', 12); 
set(gcf, 'color', 'white');
%}

end