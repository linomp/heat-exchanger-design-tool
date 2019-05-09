
% ---------------------------  CONTINUOUS G.A. ----------------------------
%   Author: Lino Mediavilla
%   Contact: lino.mediavilla@estud.usfq.edu.ec
%   Theoretical & Computational Mechanics Group 
%   Universidad San Francisco de Quito 
%   2017 - 2018
%--------------------------------------------------------------------------

function [population, candidates, avgs,bests] = theGARoutine(problem, population)  
     
    % Initialize other important variables
    problem.converged = 0;
    problem.plotOpts.terminated = 0;
    problem.current_generation = 1;
    problem.plotOpts.lastGen = 0;
    problem.constrain = 0; 
    problem.PCAconstrained = 0; 
    
    tic
    
    %--------------------- 1st generation fitness evaluation --------------
    if(problem.int_search)
       population.theIndividuals = population.roundAll(population.theIndividuals); 
    end 
    population.theIndividuals = population.checkBoundsAll(population.theIndividuals); 
    population.theIndividuals = population.applyUserConstraints(population.theIndividuals, problem.constraints_func);

    for i = 1 : population.popSize  
        population.theIndividuals{i}.fitness =  1/(problem.objective_func(population.theIndividuals{i}.chrom, problem.bc));
        problem.ff_evals = problem.ff_evals + 1;    
    end  
    
    % Identify the best individual to avoid replacing it
    % (elitism). DeSantos(pg.36) 
    population.updateBestIndividual(); 

    avg_pop_fitness = population.avgPopFitness();
    optimum_objective_func_value = population.bestIndividual.fitness;  
    %----------------------------------------------------------------------
    
    %------------------------- 1st generation stats -----------------------
    clc 
    disp('------------------------------------------------')
    disp('GA PROGRESS')
    fprintf('Elapsed time: %f\n', toc)
    fprintf('Current Generation: %d\n', problem.current_generation);
    fprintf('Obj. Func. Evaluations: %d\n', problem.ff_evals);
    fprintf('Best minimum found: %d\n', optimum_objective_func_value);
    fprintf('Avg. Population Fitness: %f\n', avg_pop_fitness);  
    problem.decoding_func(population.bestIndividual);
    
    problem.plotEvolution(population.theIndividuals, population.bestIndividual);
    problem.bestParameters = population.bestIndividual; 
    %---------------------------------------------------------------------- 
       
    
    for current_generation = 2:problem.maxGens   
        
       problem.current_generation = current_generation; 

       %------------------------ Convergence Criteria ---------------------
       % The following criteria will be considered:
       % 1. min error function value < tolerance  
       % 2. a number of generations has passed without achieving further
       %    improvement in the best individual fitness
       % 3. max number of generations has been reached

       if(optimum_objective_func_value < problem.tolerance)
         problem.converged = 1;
         break
       end  

       obj_func_vals(problem.current_generation) = optimum_objective_func_value; 
       if (problem.current_generation >= problem.S)  
           delta = [];
           for i = problem.current_generation-1:-1:(problem.current_generation-problem.S+1)
              delta = [delta; abs(obj_func_vals(problem.current_generation) - obj_func_vals(i))];
           end 
           if max(delta) < problem.delta_tol
               problem.converged = 1;
           end
       end 

       if(problem.converged || problem.current_generation == problem.maxGens)  
           break
       end
       %-------------------------------------------------------------------       
       
       % Didn't converge, go on..
       problem.current_generation = current_generation;        
      
      
       %-------------- Form pool from Crossover & Mutation offspring ------
       
       population.makeMatingPool();
       
       crossoverOffspring = population.doCrossover();
       
       % Haupt mutation method offspring
       mutationOffspring1 = population.doMutation();
       
       % Gaussian mutation offspring
       mutationOffspring2 = population.doGaussianMutation(crossoverOffspring);
       
       % Make a pool containing mutation offspring, current population and
       % crossover offspring.
       population.pool = { population.theIndividuals{:,:}, ...
                           mutationOffspring1{:,:}, ...
                           mutationOffspring2{:,:}, ...
                           crossoverOffspring{:,:} };  
       
       %-------------------------------------------------------------------    
       
       
       %------------------------ Fitness evaluation------------------------ 
       if(problem.int_search)
          population.pool = population.roundAll(population.pool); 
       end 
       population.pool = population.checkBoundsAll(population.pool);
       
       population.pool = population.applyUserConstraints(population.pool, problem.constraints_func); 
       
       for i = 1:length(population.pool)   
           population.pool{i}.fitness =  1/problem.objective_func(population.pool{i}.chrom,problem.bc); 
           problem.ff_evals = problem.ff_evals + 1; 
       end  
       
       % Update best individual 
       population.updateBestIndividual(); 
       problem.plotEvolution(population.pool, population.bestIndividual); 
       clc       
       %-------------------------------------------------------------------       
       
       
       %-------------------- Selection & New Population -------------------    
       
       % Get a cell array of readily-constructed Individual objects
       newInds = population.getNewInds(population.pool);  
       
       % Don't construct anything: just update the array of Individuals 
       % (data member of Population object).    
       population.updatePopulation(newInds);   
       %------------------------------------------------------------------- 
       
       
       %-------------- Plot and Display results in every iteration -------- 
       clc       
       population.updateBestIndividual();  
       avg_pop_fitness = population.avgPopFitness();
       optimum_objective_func_value = population.bestIndividual.fitness;  
       
       avgs(current_generation) = avg_pop_fitness;
       bests(current_generation) = population.bestIndividual.fitness;  
       
       disp('------------------------------------------------')
       disp('GA PROGRESS') 
       fprintf('\nMating pool size: %d', length(population.matingPool));
       fprintf('\nPool size: %d', length(population.pool));
       fprintf('\nhaupt mutants: %d', length(mutationOffspring1));
       fprintf('\ngaussian mutants: %d\n', length(mutationOffspring2));
       fprintf('Elapsed time: %f\n', toc);
       fprintf('Current Generation: %d\n', problem.current_generation);
       fprintf('Obj. Func. Evaluations: %d\n', problem.ff_evals);
       fprintf('Best objective function value: %d\n', optimum_objective_func_value);
       fprintf('Avg. Population Fitness: %d\n', avg_pop_fitness);  
       problem.decoding_func(population.bestIndividual);
       
       problem.plotOpts.lastGen = (current_generation == (problem.maxGens - 1));      
       %-------------------------------------------------------------------        
       
       %if problem.plotOpts.quiet
        avgs(current_generation) = 1/avg_pop_fitness;
        bests(current_generation) = 1/population.bestIndividual.fitness; 

           %{
           if(current_generation > 2)
               figure(2)
               
               semilogy(1:length(bests(2:end)), bests(2:end), 'g-', 'linewidth', 4);
               hold on
               semilogy(1:length(avgs(2:end)), avgs(2:end), 'r--', 'linewidth', 4);
               hold off
               L = legend('Best Individual Fitness', 'Avg. Population Fitness');
               set(L, 'location', 'best', 'fontSize', 12);
               drawnow
           end
            %}
       %end  
    end  
    
    candidates = getCandidatesFromPool(population.pool, problem.n_candidates);
   
end