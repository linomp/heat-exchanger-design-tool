
% -------------------------------------------------------------------------
%   Author: Lino Mediavilla
%   Contact: lino.mediavilla@estud.usfq.edu.ec
%   Theoretical & Computational Mechanics Group 
%   2017-2018
% -------------------------------------------------------------------------

classdef Population < handle
    
    % Population class
    % keeps an array of Individuals
    % hosts methods to perform selection, crossover and mutation
    
    properties
        theIndividuals;  
        pool;  
        matingPool;
        bestIndividual; 
        popSize; 
        nVars;
        nMutants; 
        tSize; 
        poolSize;
        probCrossover
        gMutProb = 0; 
        gMutSTD = 0;
        UB;
        LB;
    end
    
    methods
        
        % Constructor ( no overload )
        function obj = Population(inds, LB, UB, mutRate, gMutProb, gMutSTD, tSize, poolSize, pc, optimizationMode)
            obj.theIndividuals = cell (1, length(inds));
            for i = 1:length(inds) 
                obj.theIndividuals{i} = Individual(inds(i,:), LB, UB, optimizationMode);
            end
            obj.UB = obj.theIndividuals{1}.UB;
            obj.LB = obj.theIndividuals{1}.LB;
            obj.popSize = length(obj.theIndividuals);
            obj.nVars = obj.theIndividuals{1}.nVars;
            obj.nMutants = ceil((obj.popSize - 1)*mutRate*obj.nVars);
            obj.gMutProb = gMutProb;
            obj.gMutSTD = gMutSTD;
            obj.tSize = tSize;
            obj.poolSize = poolSize;
            obj.probCrossover = pc;
            obj.pool = {};
            obj.matingPool = {}; 
        end
        
        
        % --- Main GA methods ---
        
        % Invokes the mutation method proposed in Haupt (2004) on random 
        % individuals, which returns mutants based on the original 
        % chromosomes (new objects).
        function mutOffspring = doMutation(obj)
           mutOffspring = cell(1, obj.nMutants);
           for i = 1 : obj.nMutants  
               ind = obj.theIndividuals{ randi([1 obj.popSize]) }; 
               mutOffspring{i} = ind.mutate();
           end
        end
        
        
        % Invokes the Gaussian mutation method on random individuals, which returns
        % mutants based on the original chromosomes (object copies).
        function gMutOffspring = doGaussianMutation(obj, inds)
           gMutOffspring = {};
           c = 0;
           for i = 1 : length(inds)  
               %ind = obj.theIndividuals{ randi([1 obj.popSize]) }; 
               if(rand(1,1) < obj.gMutProb && length(inds) > 0)
                  c = c + 1;
                  gMutOffspring{c} = inds{ randi([1 length(inds)]) }.gMutate(obj.gMutProb, obj.gMutSTD);
               end
           end
        end
        
        
        % Selection routine, in each call builds a tournament and returns 
        % a winner 
        function winner = doSelection(obj, inds)
           idxs = randperm(length(inds), obj.tSize);  % get tSize random indexes
           tournament = { inds{1, idxs} }; % build tournament
           winner = tournament{1};
           for i = 2 : obj.tSize
                ind = tournament{i};
                if ind > winner  % operator '>' is overloaded for Individuals
                    winner = ind;
                end
           end            
        end
         
        
        % Pool consists of mutation and crossover offspring (parents 
        % selection is based on Tournament Selection Mode)
        function offspring = doCrossover(obj)
           offspring = {}; 
           maxMatings = round(obj.probCrossover * obj.popSize);
           for i = 1 : maxMatings      
               parent_1 = obj.doSelection(obj.matingPool);
               parent_2 = obj.doSelection(obj.matingPool);  
               [child_1, child_2] = parent_1.mateWith(parent_2);               
               offspring = {offspring{:,:},  Individual(child_1, obj.theIndividuals{1}.LB, obj.theIndividuals{1}.UB, parent_1.optimizationMode)}; 
               offspring = {offspring{:,:},  Individual(child_2, obj.theIndividuals{1}.LB, obj.theIndividuals{1}.UB, parent_1.optimizationMode)}; 
           end  
        end 
        
        % Pool consists of mutation and crossover offspring (parents 
        % selection is based on Tournament Selection Mode)
        function makeMatingPool(obj)
           obj.matingPool = {}; 
           while length(obj.matingPool) < obj.poolSize                                
               obj.matingPool = {obj.matingPool{:,:}, obj.doSelection(obj.theIndividuals)}; 
           end 
        end 
        
        function makePool(obj, inds)
            obj.pool = inds;
        end
        
        
        % Performs Tournament Selection on a previously formed pool and
        % returns a cell array of new Individuals
        function newIndividuals = getNewInds(obj, thePool) 
           newIndividuals = {obj.bestIndividual};            
           while size(newIndividuals, 2) < obj.popSize  
             % newIndividuals = {newIndividuals{:,:}, obj.doSelection(thePool)};
             newIndividuals{end+1} = obj.doSelection(thePool);
           end        
        end
        
        
        % To update the Individuals cell array
        function updatePopulation(obj, inds) 
           obj.theIndividuals = inds;            
        end    
        
        
                
        % --- Other utility methods ---   
        
        % Find best Individual routine. Useful in the elitist GA scheme
        function updateBestIndividual(obj)
            obj.bestIndividual = obj.theIndividuals{1}; 
            for i = 2 : obj.popSize
                % could be optimized as: "bestInd = max(ind{i}, bestInd)"
                ind = obj.theIndividuals{i};
                if ind > obj.bestIndividual 
                    obj.bestIndividual = ind; 
                end
            end
        end
        
        
        % Find average population fitness
        function avgFtns = avgPopFitness(obj)
            avgFtns = 0;
            for i = 1 : obj.popSize 
                avgFtns = avgFtns + obj.theIndividuals{i}.fitness;
            end
            avgFtns = avgFtns / obj.popSize;
        end
        
                
        % Round Individuals, used in Integer Search
        function rounded = roundAll(~, inds)
            rounded = cell(1,length(inds));
            for i = 1 : length(inds)
                rounded{i} = inds{i}.intChrom();
            end
        end
                
        % Checks and ajusts an input chromosome to the upper and lower 
        % bounds if neccessary. Returns an adjusted copy. Caller routine
        % should perform the overwrite.
        function adjusted = checkBounds(obj, chr)
            for i = 1:obj.nVars
%                 if( chr(i) > obj.theIndividuals{1}.UB(i) )
%                    chr(i) = obj.theIndividuals{1}.UB(i);
%                 elseif( chr(i) < obj.theIndividuals{1}.LB(i) )
%                    chr(i) = obj.theIndividuals{1}.LB(i);
%                 end
                if( chr(i) > obj.UB(i) )
                   chr(i) = obj.UB(i);
                elseif( chr(i) < obj.LB(i) )
                   chr(i) = obj.LB(i);
                end
            end
            adjusted = chr;
        end        
        function inds = checkBoundsAll(obj, inds)
            for i = 1:length(inds) 
                inds{i}.chrom = obj.checkBounds(inds{i}.chrom);
            end
        end
                       
        function inds = applyUserConstraints(obj, inds, userConstraintFunc)
            for i = 1:length(inds) 
                inds{i}.chrom = userConstraintFunc(inds{i}.chrom);
            end
        end
        
        % --- Getters & setters ---            
        
        % theIndividuals cell array
        function set.theIndividuals(obj, val)
            obj.theIndividuals = val;
        end 
        
        % Best individual
        function set.bestIndividual(obj, val)
            obj.bestIndividual = val;
        end
        
        % Mating pool
        function set.matingPool(obj, val)
            obj.matingPool = val;
        end
        
    end
    
end

