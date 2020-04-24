
% -------------------------------------------------------------------------
%   Author: Lino Mediavilla
%   Contact: lino.mediavilla@estud.usfq.edu.ec
%   Theoretical & Computational Mechanics Group 
%   2017-2018
% -------------------------------------------------------------------------
    
function [population, best_candidates, avgs, bests] = optimizationLauncher(problem)  

    %------------------------- General Settings ---------------------------

    problem.tolerance = -1;         % Maximum error function value tolerance   
    problem.maxGens = 30;           % Maximum number of generations allowed     
    problem.S = problem.maxGens;    % Number of generations to stop GA if the best minimum hasn't improved.
    problem.delta_tol = 1e-3;       % Minimum Obj. Func. Val. improvement 
    problem.optimizationMode = 'max';
    problem.inherit_ftns = 0;

    % General plot options    
    problem.plotOpts.quiet = 1;     % value of 1 overrides everything; no plots 
    problem.plotOpts.ff_profiling = 0;     
    problem.plotOpts.plots_style = 'dark';
    problem.plotOpts.plots_textIntepreter = 'none'; 

    % Video options
    problem.vidsFrameRates = 8;
    problem.vidsQuality = 85;
    problem.GA_vidName = 'test.avi';    % GA video file name 
    %----------------------------------------------------------------------


    %--------------------------- GA Settings ------------------------------

    problem.int_search = 0; % make GA search only integers 
    problem.structured_initial_pop = 1;
    problem.inherit_ftns = 1; 

    problem.tSize = 3;              % tournament size
    problem.n_tournaments = 80;     % mating pool size = n_tournaments 
    problem.pc = 0.3;               % crossover probability
    problem.mutation_rate = 0.8;    % mutation method proposed in Haupt (2004)   
    problem.mutProb = 0.80;         % for gaussian mutation
    problem.mutSTD = 1e-3;          % std. deviation for gaussian mutation

    % 1ST GENERATION IS AN STRUCTURED GRID OF DESIRED DIMENSIONS   
    problem.n_vars = length(problem.grid_dims);     % # of parameters to estimate -> # of Search space dimensions  

    % GA Plot options 
    problem.plotOpts.plot_obj_func_evolution = 0;
    problem.plotOpts.plot_full_results = 1;  

    % Disable all visualization for 3+ vars.
    if problem.n_vars > 2 || problem.plotOpts.quiet
        problem.plotOpts.plot_full_results = 0; 
        problem.plotOpts.plot_obj_func_evolution = 0;
        problem.plotOpts.plot_gbm_path = 0;
        problem.plotOpts.plot_pca_results = 0;  
    end 

    %-------------------- Objective Function Definition -------------------

    if isempty(problem.decoding_func)
        problem.decoding_func = @decodeChromosomeDefault;
    end  
    if isempty(problem.constraints_func)
        problem.constraints_func = @(chrom) chrom;
    end  


   %----------------------------- GA EXECUTION ------------------------------

    % calculate number of individuals from the grid specified above
    problem.n_inds = prod(problem.grid_dims);

    % popInit returns a set of chromosomes to initialize the first
    % generation.
    firstGenChroms = problem.popInit(); 

    % Create a new instance of Population, with a set of chromosomes, upper
    % & lower variable bounds, and misc. GA parameters.
    population = Population(firstGenChroms, problem.LB, problem.UB, ...
                            problem.mutation_rate, problem.mutProb, problem.mutSTD, ... 
                            problem.tSize, problem.n_tournaments, problem.pc, problem.optimizationMode);          

    % GO!
    [population, best_candidates, avgs, bests] = theGARoutine(problem, population);   
 

end
    
    