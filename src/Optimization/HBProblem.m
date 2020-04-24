
% -------------------------------------------------------------------------
%   Author: Lino Mediavilla
%   Contact: lino.mediavilla@estud.usfq.edu.ec
%   Theoretical & Computational Mechanics Group 
%   2017-2018
% -------------------------------------------------------------------------

classdef HBProblem < handle
    % HBProblem Class, hosts utility methods and fields that will be accessed
    % throughout the entire backanalysis procedure. 
    
    %  YEP, THERE'S A SHIIIITLOAD OF FIELDS, but this sort of interface is
    %  necessary to keep everything modular. Sorry , not sorry!!
    properties
        tolerance;
        maxGens;
        max_GBM_iter;
        S;
        delta_tol;
        ff_evals;
        plotOpts;
        int_search;
        structured_initial_pop;
        inherit_ftns;
        tSize;
        n_tournaments;
        poolSize;
        pc;
        mutation_rate;
        mutProb;
        mutSTD;
        grid_dims;
        n_vars;
        err_func_type;
        mod_func_type;
        GA_vidName;
        GBM_vidName;
        post_proc;
        hybrid;
        gbm_delta_tol;
        gbm_S;
        lambda;
        local_opt;
        UB;
        LB;
        global_min; 
        objective_func;
        mod_func;
        surfView;
        sig;
        X_me;
        C_x;
        C_x_cv2;
        AF;
        PCAconstrained;
        n_inds;
        V_GA;
        V_GBM;
        vidsFrameRates;
        vidsQuality;
        converged;
        current_generation;
        constrain;
        mem;
        bestParameters;
        best_obj_func_val;
        converged_at_generation;
        GA_FFE;
        seed;
        PCA_FFE;
        coeffs;
        latent;
        ellipse_center; 
        Local_Opt_FFE;
        optimizationMode;
        constraints_func;
        decoding_func;
        bc;
        n_candidates;
    end
    
    
    methods
        
        % Constructor 
        function obj = HBProblem()
            obj.n_candidates = 0;
            obj.ff_evals = 0;
            obj.tolerance = 0;
            obj.maxGens = 0;
            obj.max_GBM_iter = 0;
            obj.S = 0;
            obj.delta_tol = 0;
            obj.ff_evals = 0;
            obj.int_search = 0;
            obj.structured_initial_pop = 0;
            obj.inherit_ftns = 0;
            obj.tSize = 0;
            obj.n_tournaments = 0;
            obj.poolSize = 0;
            obj.pc = 0;
            obj.mutation_rate = 0;
            obj.mutProb = 0;
            obj.mutSTD = 0;
            obj.grid_dims = [];
            obj.n_vars = 0;
            obj.err_func_type = '';
            obj.mod_func_type = '';
            obj.GA_vidName = '';
            obj.GBM_vidName = '';
            obj.post_proc = 0;
            obj.hybrid = 0;
            obj.gbm_delta_tol = 0;
            obj.gbm_S = 0;
            obj.lambda = 0;
            obj.local_opt = '';
            obj.UB = [];
            obj.LB = [];
            obj.global_min = [];
            obj.objective_func = [];
            obj.mod_func = [];
            obj.surfView = [];
            obj.sig = [];
            obj.X_me = [];
            obj.C_x = [];
            obj.C_x_cv2 = [];
            obj.AF = 0;
            obj.PCAconstrained = 0;
            obj.n_inds = 0;
            obj.V_GA = [];
            obj.V_GBM = [];
            obj.vidsFrameRates = 0;
            obj.vidsQuality = 0;
            obj.converged = 0;
            obj.current_generation = 0;
            obj.constrain = 0;
            obj.mem = [];
            obj.bestParameters = Individual();
            obj.best_obj_func_val = 0;
            obj.converged_at_generation = 0;
            obj.GA_FFE = 0;
            obj.seed = [];
            obj.PCA_FFE = 0;
            obj.coeffs = [];
            obj.latent = [];
            obj.ellipse_center = []; 
            obj.Local_Opt_FFE = 0; 
            obj.optimizationMode = 'max';
            obj.constraints_func = [];
            obj.decoding_func = [];
            obj.bc = [];
            obj.plotOpts = struct(  'quiet',0, ...
                                    'ff_profiling', 0, ...
                                    'plots_style', '', ...
                                    'plots_textIntepreter', '', ...
                                    'plot_obj_func_evolution', 0, ...
                                    'plot_full_results', 0, ...
                                    'plot_pca_results', 0, ...
                                    'plot_gbm_path', 0, ...
                                    'surfFig', 0, ...
                                    'bounds', [], ...
                                    'lwidths', 0, ...
                                    'terminated', 0, ...
                                    'lastGen', 0 ...
                                );          
        end
                
        
        % ---- Utility methods ----        
                
        % Generate equally spaced 1st generation chromosomes
        function firstGenChroms = popInit(obj) 
            
            if(obj.structured_initial_pop == 0)
                %Random initial population instead of grid
                firstGenChroms = zeros(obj.n_inds,obj.n_vars);  
                for i = 1:obj.n_vars 
                  rand_genes = obj.LB(i) + (obj.UB(i)-obj.LB(i)).*rand(size(firstGenChroms,1),1);
                  firstGenChroms(:,i) = rand_genes;
                end

            else
                %Structured initial population
                %In this GA Scheme, an structured initial population is generated as a grid
                %of user-specified dimensions. This could help sample the search space in a
                %better way and offersthe possiblity of applying some extra knowledge about 
                %the problem. First, the step for each dimension is calculated, since
                %variables may lie in domains of different sizes. Then an orthogonal base is
                %generated for the search space, and then the complete set of points is
                %obtained by generating all the possible combinations of the base vectors.
                step = zeros(obj.n_vars,1);   
                for i = 1:obj.n_vars
                    step(i) = abs(obj.UB(i)-obj.LB(i))/(obj.grid_dims(i)-1);
                end 

                grid_vector_1 = (obj.LB(1):step(1):obj.UB(1))';
                for i = 2:obj.n_vars   
                    grid_vector_2 = (obj.LB(i):step(i):obj.UB(i))';
                    if i == 2 
                        grid_pop = combvec(grid_vector_1',grid_vector_2')';
                    else
                        grid_pop = combvec(grid_pop',grid_vector_2')';
                    end
                end

                %The population is a matrix in which each row corresponds to the entire
                %chromosome of an individual (# of rows = # of individuals). And each column
                %contains a gene (# of columns = # of target parameters).
                firstGenChroms = grid_pop;
            end 

        end         
        
        
        % Plots objective function surface, minima, & does setup for video
        % writing (both GA and GBM stages)
        function visualize( obj )
           obj.plotOpts.surfFig = 1;
           obj.plotOpts.bounds = [min(obj.LB) max(obj.UB)];
           obj.plotOpts.lwidths = 2.5;  
           
           if(obj.plotOpts.plot_full_results)  % TODO - ADD GUARD AGAINST LARGE SEARCH SPACE
               
               obj.V_GA = VideoWriter(obj.GA_vidName);      % GA video
               obj.V_GA.FrameRate = obj.vidsFrameRates;
               obj.V_GA.Quality = obj.vidsQuality;
           
               fig1 = figure(obj.plotOpts.surfFig); 
               %set(fig1,'Position',[-10 15 683 768]);        
               hold on; grid on

               % Plot Objective Function Surface
               mshgridpts = 15e1;
               [par1,par2] = meshgrid(linspace(obj.plotOpts.bounds(1), obj.plotOpts.bounds(2), mshgridpts));
               obj_func_surf = zeros(size(par1,1),size(par1,2));
               for i = 1:size(par1,1)
                   for j = 1:size(par1,2)
                       obj_func_surf(i,j) = obj.objective_func([par1(i,j) par2(i,j)], obj); 
                   end
               end       
               surf(par1,par2,obj_func_surf,'linestyle','none')

               % Plot global minima
               for i = 1 : size(obj.global_min, 1)
                   if( obj.global_min(i,1) > obj.plotOpts.bounds(1) ) && ( obj.global_min(i,1) < obj.plotOpts.bounds(2) ) ...
                           && ( obj.global_min(i,2) > obj.plotOpts.bounds(1) ) && ( obj.global_min(i,2) < obj.plotOpts.bounds(2) )
                       z_min = obj.objective_func(obj.global_min(i,:), obj);
                       plot3(obj.global_min(i,1), obj.global_min(i,2), z_min, 'mp', 'markersize', 15, 'linewidth', 2); %, 'markerfacecolor', 'm', );
                   end
               end

               tcmgPlotFormat('', 'Parameter 1', 'Parameter 2', 'Objective function',  obj.plotOpts.plots_style, '', obj.plotOpts.plots_textIntepreter);
               %axis([obj.plotOpts.bounds obj.plotOpts.bounds]); axis manual
               axis tight; axis manual
               view(obj.surfView)
               alpha(0.75) 
               camlight

               % Open new file to write GA frames (if GA plot is turned
               % off, existing file won't get overwritten)
               open(obj.V_GA);
           end
           
           if (obj.plotOpts.plot_gbm_path) 
               obj.V_GBM = VideoWriter(obj.GBM_vidName);    
               obj.V_GBM.FrameRate = floor(obj.vidsFrameRates*.75);
               obj.V_GBM.Quality = obj.vidsQuality;
               % Open new file to write GBM frames (if GBM path plot is turned
               % off, existing file won't get overwritten)
              open(obj.V_GBM);
           end           
        end               
        
                
        % Designed to plot GA individuals as the evolve, over an objective
        % function surface
        function plotEvolution(obj, inds, best_ind) 
            if(obj.plotOpts.plot_full_results)     % TODO - FIX PLOTTING W.O. SURFACE   
                h = []; 
                obj.plotOpts.surfFig = 1;
                figure(obj.plotOpts.surfFig) 
                for i = 1 : length(inds)
                    h(i) = plot3(inds{i}.chrom(1), inds{i}.chrom(2), 1/inds{i}.fitness, 'ko', ...
                          'linewidth', 2.35,'markerfacecolor', 'r', 'markersize', 9.5);    
                end  

                h(end+1) = plot3(best_ind.chrom(1), best_ind.chrom(2), 1/best_ind.fitness, 'gd', ...
                          'linewidth', 2.5, 'markersize', 14); %'markerfacecolor', 'k');
                padding = 1;
                axis([min(inds{1}.LB)-padding max(inds{1}.UB)+padding min(inds{1}.LB)-padding max(inds{1}.UB)+padding]) 
                drawnow
                writeVideo(obj.V_GA, getframe); 

                if(~obj.plotOpts.lastGen && ~obj.plotOpts.terminated)
                   delete(h)
                end
            end
            
        end
        
        
        % Set objective function (to test), surface angle, variable bounds
        % and test function global minima
        function defineObjectiveFunc(obj)
    
            % Defaults to be overriden 
            obj.surfView = [35 -20];  
            local_UB = [10 10]; local_LB = [-10 -10]; 

            if (strcmp(obj.mod_func_type,'f'))
                obj.global_min = zeros(1,obj.n_vars);
                obj.surfView = [-25 50];
                obj.mod_func = @f_adaptive_test_OF; % test function f

            elseif(strcmp(obj.mod_func_type,'g'))
                obj.global_min = -1.*ones(1,obj.n_vars);
                obj.mod_func = @g_adaptive_test_OF; % test function g

            elseif(strcmp(obj.mod_func_type,'h'))
                obj.global_min = ones(1,obj.n_vars);
                obj.mod_func = @h_adaptive_test_OF; % test function h

            elseif (strcmp(obj.mod_func_type,'haupt'))
                obj.global_min = 10.*[0.9039, 0.8668; -0.9039, 0.8668; 0.9039, -0.8668; -0.9039, -0.8668];
                obj.surfView = [40 -10];
                local_UB = [10 10]; local_LB = [0 0];
                obj.mod_func = @haupt;

            elseif (strcmp(obj.mod_func_type,'rosenbrock'))
                obj.global_min = ones(1,obj.n_vars);
                obj.mod_func = @rosenbrock;

            elseif (strcmp(obj.mod_func_type,'cross')) 
                obj.global_min = [1.34941,1.34941;-1.34941,1.34941; 1.34941,-1.34941; -1.34941,-1.34941];
                obj.surfView = [35 -15]; 
                local_UB = [15 15]; local_LB = -1*local_UB;
                obj.mod_func = @cross;

            elseif (strcmp(obj.mod_func_type,'kursawe'))
                obj.global_min = [-1.151541601899755  -1.152781473642893];
                obj.surfView = [105 5];
                local_UB = [10 10]; local_LB = -1*local_UB; 
                obj.mod_func = @kursawe;

            elseif (strcmp(obj.mod_func_type,'h1'))
                obj.global_min = [8.6998, 6.7665];
                obj.surfView = [35 25];  
                obj.mod_func = @h1;

            elseif (strcmp(obj.mod_func_type,'schaffer'))
                obj.global_min = [0,0];
                obj.surfView = [35 -10];
                obj.mod_func = @schaffer;

            elseif (strcmp(obj.mod_func_type,'table')) 
                obj.global_min = [8.05502,9.66459 ;-8.05502,9.66459 ; 8.05502,-9.66459; -8.05502,-9.66459];
                obj.surfView = [65 25];        
                obj.mod_func = @table_fun;

            elseif (strcmp(obj.mod_func_type,'modschaffer1'))
                obj.global_min = [0,0];
                obj.surfView = [35 -45];        
                local_UB = [3 3]; local_LB = -1*local_UB;
                obj.mod_func = @modschaffer1;

            elseif (strcmp(obj.mod_func_type,'schweffel'))
                obj.global_min = [0,0];
                obj.mod_func = @schweffel;

            elseif (strcmp(obj.mod_func_type,'peaks'))
                obj.global_min = [0.227744735326231  -1.628074888356166];
                obj.surfView = [-45 20];
                local_UB = [3 3]; local_LB = -1*local_UB;
                obj.mod_func = @zpeaks;

            end        

            if (strcmp(obj.err_func_type,'RE'))
                obj.objective_func = @error_func_RE;
            elseif(strcmp(obj.err_func_type,'LS'))
                obj.objective_func = @error_func_LS;
            elseif(strcmp(obj.err_func_type,'ML'))
                obj.objective_func = @error_func_ML;
            end  

            if(obj.n_vars == 2)
                obj.UB = local_UB;
                obj.LB = local_LB;
            end
            
        end
        
                        
        
        % ---- Getters & Setters ----
        % A SHITLOAD OF THESE AS WELL. All rutinary one-liners.
        
        function set.tolerance(obj, val) 
            obj.tolerance = val; 
        end
        function set.maxGens(obj, val) 
            obj.maxGens = val; 
        end
        function set.max_GBM_iter(obj, val) 
            obj.max_GBM_iter = val; 
        end
        function set.S(obj, val) 
            obj.S = val; 
        end
        function set.delta_tol(obj, val) 
            obj.delta_tol = val;
        end
        function set.ff_evals(obj, val) 
            obj.ff_evals = val; 
        end
        function set.plotOpts(obj, val) 
            obj.plotOpts = val; 
        end
        function set.int_search(obj, val) 
            obj.int_search = val; 
        end
        function set.structured_initial_pop(obj, val) 
            obj.structured_initial_pop = val; 
        end
        function set.inherit_ftns(obj, val) 
            obj.inherit_ftns = val; 
        end
        function set.tSize(obj, val) 
            obj.tSize = val; 
        end
        function set.n_tournaments(obj, val) 
            obj.n_tournaments = val; 
        end
        function set.poolSize(obj, val) 
            obj.poolSize = val; 
        end
        function set.pc(obj, val)
            obj.pc = val;
        end
        function set.mutation_rate(obj, val) 
            obj.mutation_rate = val; 
        end
        function set.mutProb(obj, val) 
            obj.mutProb = val; 
        end
        function set.mutSTD(obj, val) 
            obj.mutSTD = val; 
        end
        function set.grid_dims(obj, val) 
            obj.grid_dims = val; 
        end
        function set.n_vars(obj, val) 
            obj.n_vars = val; 
        end
        function set.err_func_type(obj, val) 
            obj.err_func_type = val; 
        end
        function set.mod_func_type(obj, val) 
            obj.mod_func_type = val; 
        end
        function set.GA_vidName(obj, val) 
            obj.GA_vidName = val; 
        end
        function set.GBM_vidName(obj, val) 
            obj.GBM_vidName = val; 
        end
        function set.post_proc(obj, val) 
            obj.post_proc = val; 
        end
        function set.hybrid(obj, val) 
            obj.hybrid = val; 
        end
        function set.gbm_delta_tol(obj, val) 
            obj.gbm_delta_tol = val; 
        end
        function set.gbm_S(obj, val) 
            obj.gbm_S = val; 
        end
        function set.lambda(obj, val) 
            obj.lambda = val; 
        end
        function set.local_opt(obj, val) 
            obj.local_opt = val; 
        end
        function set.UB(obj, val) 
            obj.UB = val; 
        end
        function set.LB(obj, val) 
            obj.LB = val; 
        end
        function set.global_min(obj, val) 
            obj.global_min = val; 
        end 
        function set.objective_func(obj, val) 
            obj.objective_func = val; 
        end
        function set.mod_func(obj, val) 
            obj.mod_func = val; 
        end
        function set.surfView(obj, val) 
            obj.surfView = val; 
        end
        function set.sig(obj, val) 
            obj.sig = val; 
        end
        function set.X_me(obj, val) 
            obj.X_me = val; 
        end
        function set.C_x(obj, val) 
            obj.C_x = val; 
        end
        function set.C_x_cv2(obj, val) 
            obj.C_x_cv2 = val; 
        end
        function set.AF(obj, val) 
            obj.AF = val; 
        end
        function set.PCAconstrained(obj, val) 
            obj.PCAconstrained = val; 
        end
        function set.n_inds(obj, val) 
            obj.n_inds = val; 
        end
        function set.V_GA(obj, val) 
            obj.V_GA = val; 
        end
        function set.V_GBM(obj, val) 
            obj.V_GBM = val; 
        end
        function set.vidsFrameRates(obj, val)
            obj.vidsFrameRates = val;
        end
        function set.vidsQuality(obj, val)
            obj.vidsQuality = val;
        end
        function set.converged(obj, val) 
            obj.converged = val; 
        end
        function set.current_generation(obj, val) 
            obj.current_generation = val; 
        end
        function set.constrain(obj, val) 
            obj.constrain = val; 
        end
        function set.mem(obj, val) 
            obj.mem = val; 
        end
        function set.bestParameters(obj, val) 
            obj.bestParameters = val; 
        end
        function set.best_obj_func_val(obj, val) 
            obj.best_obj_func_val = val; 
        end
        function set.converged_at_generation(obj, val) 
            obj.converged_at_generation = val; 
        end
        function set.GA_FFE(obj, val) 
            obj.GA_FFE = val; 
        end
        function set.seed(obj, val) 
            obj.seed = val; 
        end
        function set.PCA_FFE(obj, val) 
            obj.PCA_FFE = val; 
        end
        function set.coeffs(obj, val) 
            obj.coeffs = val; 
        end
        function set.latent(obj, val) 
            obj.latent = val; 
        end
        function set.ellipse_center(obj, val) 
            obj.ellipse_center = val; 
        end 
        function set.Local_Opt_FFE(obj, val) 
            obj.Local_Opt_FFE = val; 
        end
        function set.optimizationMode(obj, val) 
            obj.optimizationMode = val; 
        end
        function set.constraints_func(obj, val) 
            obj.constraints_func = val; 
        end
        function set.decoding_func(obj, val) 
            obj.decoding_func = val; 
        end
        function set.bc(obj, val) 
            obj.bc = val; 
        end
        function set.n_candidates(obj, val) 
            obj.n_candidates = val; 
        end
    end
    
end

