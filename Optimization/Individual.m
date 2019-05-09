
% -------------------------------------------------------------------------
%   Author: Lino Mediavilla
%   Contact: lino.mediavilla@estud.usfq.edu.ec
%   Theoretical & Computational Mechanics Group 
%   2017-2018
% -------------------------------------------------------------------------

classdef Individual < matlab.mixin.Copyable 
    
    % Individual Class
    % Keeps a chromosome, a fitness value, and the lower/upper bounds for
    % each 'gene'.
    % super.copy method is overriden to produce deep copies of Individuals
    
    % THESE TERMS ARE USED CONSISTENTLY THROUGHOUT THE CODE:
    % chromosome: vector with potential values of design variables:
    % Individual: a class from which each object holds a chromosome
    %             as property, as well as fitness and bounds.
    
    properties
        chrom;
        LB;
        UB;
        fitness;
        nVars;
        optimizationMode;
    end
    
    methods
        
        % Constructor (contains alternative no-parameter implementation)
        function obj = Individual(chr, lb, ub, optimizationMode)
            if nargin > 0  
                obj.chrom = chr; 
                obj.LB = lb;
                obj.UB = ub; 
                obj.nVars = length(lb); 
                obj.optimizationMode = optimizationMode;
            else
                obj.chrom = [];
                obj.LB = [];
                obj.UB = [];
                obj.fitness = 0;
                obj.nVars = 0;
                obj.optimizationMode = 'max';
            end
        end        
        
        
        % Mutation method proposed in Haupt (2004).
        % Returns a mutant (new object) based on *this 
        function mutant = mutate(obj) 
            rand_col = randi([1 length(obj.chrom)]);
            rand_gene = obj.LB(rand_col) + (obj.UB(rand_col)-obj.LB(rand_col)).*rand(1,1);
            mutant = Individual(obj.chrom, obj.LB, obj.UB, obj.optimizationMode);
            mutant.chrom(rand_col) = rand_gene;  
        end
        
        
        % Gaussian Mutation
        % Returns a mutant (new object) based on *this 
        function mutant = gMutate(obj, mutProb, mutSTD) 
            mutant = Individual(obj.chrom, obj.LB, obj.UB, obj.optimizationMode); 
            for i = 1:length(mutant.chrom)
                if(mutProb > rand(1,1))
                    mutant.chrom(i) = mutant.chrom(i) + mutSTD*randn(1,1);
                end
            end 
        end
        
       
        % Takes the chromosome of another individual as input and returns 
        % two new chromosomes, according to the method for continuous GAs 
        % proposed in Haupt (2004). 
        function [of_1, of_2] = mateWith(obj, parent_2)            
            %-------BLEND CROSSOVER------- 
            % Crossover method proposed in Haupt (2004). pg 116 
            new_material = zeros(1, obj.nVars); 
            for i = 1 : obj.nVars
                beta = 1*rand(1,1);
                new_material(1,i) = beta*(obj.chrom(i)-parent_2.chrom(i));
            end 
            % Build new chromosomes and make sure the offspring is within bounds 
            of_1 = obj.chrom - new_material;
            of_2 = parent_2.chrom + new_material;       
        end  
        
        
        % Concatenates the entire chromosome of the Individual into a string
        % and returns it. This is further used as 'key' to store individuals 
        % in a HashMap, where the associated value is the fitness. 
        function str = toString(obj)
            %str = num2str(obj.chrom,32);
            if(iscell(obj.chrom))
                str = num2str(obj.chrom{1},32);
            else
                str = num2str(obj.chrom,32);
            end
            %str='';
        end
        
        
        % Returns a new individual based on this* but with all components
        % of its chromosome rounded for integer search. 
        function rounded = intChrom(obj)
            rounded = Individual(round(obj.chrom), obj.LB, obj.UB, obj.optimizationMode);
        end
                       
        
        % -- Getters & Setters --
        
        % Fitness
        function set.fitness(obj, val)
            obj.fitness = val;
        end 
        
        % Chromosome
        function val = get.chrom(obj)
            %val = obj.chrom;
            %if(iscell(obj.chrom))
            %    val = obj.chrom{1};
            %else
                val = obj.chrom;
            %end
        end
        function set.chrom(obj, val)
            %if(iscell(val))
            %    obj.chrom = val{1};
            %else
                obj.chrom = val;
            %end
        end
        
        
        % -- Overloaded operators --
        
        % Greater than - by comparing fitness 
        function res = gt(obj1, obj2)  
            criteria = obj1.fitness > obj2.fitness; 
            
            if(criteria)
                res = true;
            else
                res = false;
            end
        end
        
        % Less than - base on the previous 'gt'
        function res = lt(obj1, obj2)
            res = ~(obj1 > obj2);
        end
          
        
        % Equality operator: compares chromosomes. Returns true if equal.
        function res = eq(obj, obj2) 
            if prod(obj.chrom == obj2.chrom) == 1
                res = true;
            else
                res = false;
            end
        end
        
    end    
    
        
    % -- Overriden copy method --
    methods(Access = protected)
        function cp = copyElement(obj)
           cp = Individual(obj.chrom, obj.LB, obj.UB, obj.optimizationMode);
           cp.fitness = obj.fitness;
        end
    end
    
end

