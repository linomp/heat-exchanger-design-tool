
% -------------------------------------------------------------------------
%   Author: Lino Mediavilla
%   Contact: lino.mediavilla@estud.usfq.edu.ec
%   Theoretical & Computational Mechanics Group 
%   2017-2018
% -------------------------------------------------------------------------

classdef GAPostProc < handle
    % GAPostProc class, hosts methods for PCA and Gradient-Based Search 
    
    % TO-DO: WRITE CONJUGATE GRADIENTS METHOD TO TEST AGAINST MARQUARDT/NEWTON
    %        ALGORITHMS
    
    properties
        pcaInds; chroms;  
        pData; pcaData; gbmData;
        bestGAIndividual;
    end
    
    methods
        
        % Constructor
        function obj = GAPostProc(inds, best, problem)
            obj.pData = problem;
            obj.pcaInds = inds;
            obj.bestGAIndividual = best;            
            for i = 1 : length(obj.pcaInds)
                obj.chroms(i,:) = obj.pcaInds{i}.chrom;
            end
        end
        
        
        % ------------------ MAIN METHODS --------------------
        % Carry out PCA using Matlab function
        function [seed, FFE, coeffs, latent, e_c] = doPCA(obj)            
            
            FFE = 0; % fitness func. evals count
            obj.pcaData.ellipse_center = mean(obj.chroms);
            [obj.pcaData.coeffs, ~, obj.pcaData.latent] = pca(obj.chroms);   
                
            % Standard deviation amplification factor for the pca ellipse axes
            obj.pcaData.latent = obj.pData.AF .* obj.pcaData.latent;
            
            seed = Individual(obj.pcaData.ellipse_center, obj.pData.LB, obj.pData.UB);
            
            % Evaluate fitness of the individual initialized to the
            % coordinates of the ellipse center.
            if ~eq(seed, obj.bestGAIndividual)
                seed.fitness = 1/obj.pData.objective_func(seed.chrom, obj.pData);
                FFE = FFE + 1; 
            end
            
            % Remains as seed only iff it has a higher fitness value than 
            % the best GA individual.
            if seed < obj.bestGAIndividual 
                fprintf('\n%s\n\n%s', ...
                    'Ellipse Center yields a larger O.F. value than the best GA individual', ...
                    'Seed has been set to the best individual found by GA')        
                seed = obj.bestGAIndividual;
                obj.pcaData.ellipse_center = obj.bestGAIndividual.chrom;
            end            
            fprintf('\nSeed OF Value: %d \n', 1/seed.fitness);  
            
            % Output data assignment
            coeffs = obj.pcaData.coeffs;
            latent = obj.pcaData.latent;
            e_c = obj.pcaData.ellipse_center;
            
        end
        
        
        % --- Carry out Gradient-Based Search --- 
        
        % Manage the requested GBM
        function [best, FFE] = GBM(obj, pData)                    
            if(strcmp(pData.local_opt,'marquardt'))
                [best, FFE] = obj.marquardt(pData); 
            elseif(strcmp(pData.local_opt,'gauss-newton'))
                [best, FFE] = obj.gaussNewton(pData); 
            elseif(strcmp(pData.local_opt,'conjugate gradients'))
                [best, FFE] = obj.conjGrad(pData);
            end
        end
        
        % Marquardt Method        
        function [best, FFE] = marquardt(obj, pData)
            
            FFE = 0;
            % Deep copies ...
            best = copy(pData.seed);
            potential = copy(best);
            
            obj.gbmData.f_min = 1/best.fitness; 
            obj.gbmData.points = inf*ones(pData.max_GBM_iter, best.nVars);
            
            % For convergence check;
            %obj_func_vals = inf.*ones(pData.max_GBM_iter,1);
            gbm_delta = inf.*ones(pData.gbm_S-1,1);
            
            % Matrix to store partial derivatives
            A = zeros(length(pData.X_me), best.nVars);
            % Step along each dimension to compute partial derivative
            % TODO - generate different ffd_h according domain size ??
            ffd_h = (1e-3)*eye(best.nVars);
            
            % Auxiliary data for RE Error Function
            if(strcmp(pData.err_func_type,'RE')) 
                B = zeros(length(pData.X_me));
                for i = 1:length(pData.X_me)
                   B(i,i) = 1/pData.X_me(i);
                end
            end
            
            % Marquardt Algorithm Parameters (De Santos pg. 108)
            if(strcmp(pData.err_func_type,'LS')) 
                mu = 1e-6; %1e-1 %1
                rho = 10; %.5; % 1.5;
            else
                mu = 1e-6; %1e-1
                rho = 15; 
            end
            
            % Smaller objective function contour for gbm path
            if(obj.pData.plotOpts.plot_gbm_path)
                figure(2) 
                point = pData.seed.chrom;
                offset = (obj.pData.plotOpts.bounds(2) - obj.pData.plotOpts.bounds(1))/500;
                reduced_bounds = [ ... 
                    min( [min([obj.pData.global_min(1,1), point(1)])-offset, min([obj.pData.global_min(1,2), point(2)])-offset] ), ...
                    max( [max([obj.pData.global_min(1,1), point(1)])+offset, max([obj.pData.global_min(1,2), point(2)])+offset] ) ...
                    ];                   
                [par1,par2] = meshgrid(linspace(reduced_bounds(1), reduced_bounds(2), 8e1));
                obj_func_surf = zeros(size(par1,1),size(par1,2));
                for i = 1:size(par1,1)
                    for j = 1:size(par1,2)
                        obj_func_surf(i,j) = pData.objective_func([par1(i,j) par2(i,j)], pData); 
                    end
                end     
                contour(par1,par2,obj_func_surf, 100);
            end
            
            for n_iter = 1 : pData.max_GBM_iter  
                
                % Save gbm path to plot it
                obj.gbmData.points(n_iter,:) = best.chrom; 
                if(obj.pData.plotOpts.plot_gbm_path)
                    obj.plotGBMPath(pData.seed.chrom, 0);
                    %obj.plotGBMPath(best.chrom, 0);
                end
                
                % Evaluate Objective Function
                [~, X_calc] = pData.objective_func(best.chrom,pData);
                FFE = FFE + 1;
                
                % Calculation of FFD matrix"A"
                
                for j = 1 : best.nVars                    
                    %{
                     Each component in matrix A stores the FFD approximation of the 
                     partial derivative of the objective function with respect to a small
                     increment ffd_h(j) in parameter 'j'. 

                     Strategy:
                     try with a diagonal matrix ffd_h; add row 'j'  to the whole old_p
                     vector. Only the parameter 'j' will increase and the others will 
                     remain constant.  
                    %}                                         
                    temp = best.chrom + ffd_h(j,:); 

                    % Evaluate Objective Function with partial change in 
                    % parameters vector "p".
                    [~, X_calc_ffd] = pData.objective_func(temp ,pData); 
                    FFE = FFE + 1; 
                    for i = 1 : length(pData.X_me)         
                        num = X_calc_ffd(i) - X_calc(i);       
                        den = max(ffd_h(j,:));
                        A(i,j) = num/den;
                    end   
                end 
                                 
                % Create and initialize an individual that could potentially 
                % yield a better objective function value than the best one 
                % found so far. Its chromosomewill be modified by an
                % advance vector. 
                first_flag = 1; 
                ctr = 0;
                
                % Generate a new advance vector until a descent direction
                % is found.
                while( potential < best || first_flag) && ctr < 15;  
                    
                    % Update marquardt parameters if necessary
                    if( ~first_flag )   
                        mu = mu * rho;
                    else
                        first_flag = 0;
                    end                    
                    
                    % Calculation of advance vector
                    marquardt_matrix = mu.*eye(best.nVars); 
                    if(strcmp(pData.err_func_type,'ML'))
                        % For ML Error Function
                        delta_p = inv(((A') * inv(pData.C_x) * A) + marquardt_matrix) * (A') * inv(pData.C_x) * (pData.X_me - X_calc)';
                        
                    elseif(strcmp(pData.err_func_type,'RE')) 
                        % For RE Error Function
                        a_x = zeros(1,length(pData.sig));
                        for i = 1:length(pData.sig) 
                            a_x(i) = (pData.X_me(i) - X_calc(i)) / pData.X_me(i); 
                        end 
                        delta_p = inv(((A') * (B') * inv(pData.C_x_cv2) * B * A) + marquardt_matrix) * (A') * (B') * inv(pData.C_x_cv2) * a_x';
                        
                    elseif(strcmp(pData.err_func_type,'LS'))
                        % For LS/Markov Error Function
                        % Weighted diagonal matrix (W = I is the basic quadratic error approach).
                        W = eye(length(pData.X_me)); 
                        M1 = ((A')* (W\A) ) + marquardt_matrix; 
                        M2 = M1\(A'); 
                        delta_p = M2 * (W\(pData.X_me - X_calc)');  
                        
                    end      
                    potential.chrom = best.chrom + reshape(delta_p, size(best.chrom,1), size(best.chrom,2)); 
                    potential.fitness = 1/pData.objective_func(potential.chrom,pData);
                    if(isnan(sum(potential.chrom)))
                        break
                    end
                    FFE = FFE + 1; ctr = ctr+1; 
                end
                
                % If a descent direction was found; update best individual
                % and Marquardt parameter
                if(~isnan(sum(potential.chrom)))
                    best = copy(potential);
                end
                obj.gbmData.f_min = 1/best.fitness;
                mu = mu/rho;
                
                % Check for convergence
                obj_func_vals(n_iter) = obj.gbmData.f_min;
                if (n_iter >= pData.gbm_S)  
                    gbm_delta = [];
                    for i = n_iter-1:-1:(n_iter-pData.gbm_S+1)
                       gbm_delta = [gbm_delta; abs(obj_func_vals(n_iter) - obj_func_vals(i))];
                    end 
                end  
                if(obj.gbmData.f_min <= pData.tolerance || max(gbm_delta) <= pData.gbm_delta_tol)  
                    break
                end                 
            end 
            
        end  
        
        % Gauss-Newton Method        
        function [best, FFE] = gaussNewton(obj, pData)
            best = pData.seed;
            FFE = 0;
            % TODO ...
        end    
        
        % Conjugate Gradients Method
        function [best, FFE] = conjGrad(obj, pData)
            best = pData.seed;
            FFE = 0;
            % TODO ...
        end 
                
        
        % -----------------------------------------------------
        
        
        % ----------------- METHODS FOR PLOTS -----------------
        % --- PCA ---        
        % Plot PCA results
        function pca_fig = plotPCA(obj, pData)     
             
            pc = [];
            for k = 1:length(obj.pcaData.latent)
                q = [-obj.pcaData.latent(k) obj.pcaData.latent(k)];
                for i = 1:2
                    pc = [pc ;obj.pcaData.ellipse_center + (q(i).*obj.pcaData.coeffs(:,k))'];
                end
            end 
            a = obj.pcaData.latent(1); b = obj.pcaData.latent(2);

            %Extract info from the principal component directions
            pc1 = obj.pcaData.coeffs(:,1);  pc2 = obj.pcaData.coeffs(:,2);  
            m_pc1 = pc1(2)/pc1(1);  m_pc2 = pc2(2)/pc2(1); 

            %Traslation to ellipse center
            pc1 = obj.pcaData.coeffs(:,1) + obj.pcaData.ellipse_center';
            pc2 = obj.pcaData.coeffs(:,2) + obj.pcaData.ellipse_center';  
            x = (-10:1:10);
            l_pc1 = m_pc1.*(x - pc1(1)) + pc1(2);
            l_pc2 = m_pc2.*(x - pc2(1)) + pc2(2);

            %Constrains
            pc1_angle = atan(1/m_pc1);
            b_prime = b/sin(pc1_angle);
            
            pc2_angle = atan(m_pc2);
            a_prime = a/sin(90 - pc2_angle); 

            l_1 = l_pc1 - b_prime; l_2 = l_pc1 + b_prime;
            l_3 = l_pc2 + a_prime; l_4 = l_pc2 - a_prime; 
            
            % Begin plotting
            pca_fig = figure(2);
            tcmgPlotFormat('', 'Parameter 1', 'Parameter 2', '',  pData.plotOpts.plots_style, '', pData.plotOpts.plots_textIntepreter);
            
            hold on; grid on;
            set(pca_fig,'Position',[690 20 683 768])
            ang = atan( (pc(2,2) - pc(1,2)) / (pc(2,1) - pc(1,1)) );
            [ex, ey] = obj.plotEllipse([obj.pcaData.ellipse_center(1), obj.pcaData.ellipse_center(2), a, b, -ang]);
            
            % Ellipse
            p1 = plot(ex, ey ,'g','linewidth', pData.plotOpts.lwidths);
            
            % PCA individuals
            p2 = plot(obj.chroms(:, 1), obj.chroms(:, 2), 'ko', 'markerfacecolor', 'r', 'markersize',8,'linewidth',1.3);
            
            % Objective Function contour plot
            reduced_bounds = [obj.pcaData.ellipse_center(1)-a-2 ...
                              obj.pcaData.ellipse_center(2)+a+2];
            %[par1,par2] = meshgrid(reduced_bounds(1) : .25 : reduced_bounds(2)); 
            [par1,par2] = meshgrid(linspace(reduced_bounds(1), reduced_bounds(2), 8e1));
            obj_func_surf = zeros(size(par1,1),size(par1,2));
            for i = 1:size(par1,1)
                for j = 1:size(par1,2)
                    obj_func_surf(i,j) = pData.objective_func([par1(i,j) par2(i,j)], pData); 
                end
            end            
            p3 = contour(par1,par2,obj_func_surf, 150);
            axis([reduced_bounds reduced_bounds]); 
            
            % Ellipse center
            p4 = plot(obj.pcaData.ellipse_center(1), obj.pcaData.ellipse_center(2),'kd','MarkerSize',12,'MarkerFaceColor','g','linewidth',1.5);
            
            % Global Min
            for i = 1 : size(pData.global_min, 1)
                p5(i) = plot(pData.global_min(i,1), pData.global_min(i,2),'mp','MarkerSize',15,'linewidth',1.5);
            end 
            
        end                
        % Generate data to plot the PCA ellipse in 2-D
        function [rot_x, rot_y] = plotEllipse(obj, args)
            x_center = args(1); y_center = args(2);
            a = args(3); % major axis
            b = args(4); % minor axis
            
            % Ellipse CCW rotation angle (major axis w.r.t horizontal axis)
            alpha = args(5); 

            s = 0:pi/200:2*pi;
            ellipse_x = a*cos(s);  ellipse_y = b*sin(s);

            rot_ellipse_x = ellipse_x*cos(-alpha) - ellipse_y*sin(-alpha);
            rot_ellipse_y = ellipse_x*sin(-alpha) + ellipse_y*cos(-alpha);

            rot_x = x_center + rot_ellipse_x; 
            rot_y = y_center + rot_ellipse_y;
        end
        
        % --- GBM ---
        function plotGBMPath(obj, point, status)
            figure(2)
            hold on
            plot(obj.gbmData.points(:,1), obj.gbmData.points(:,2), '-r*', 'MarkerSize',8,'linewidth',5);
            if (status)
                plot(obj.gbmData.points(end,1), obj.gbmData.points(end,2), 'go','MarkerSize',13,'MarkerFaceColor','b','linewidth',2);
            end
            %offset = 10;%1e-2;
            offset = (obj.pData.plotOpts.bounds(2) - obj.pData.plotOpts.bounds(1))/500;
            axis([min([obj.pData.global_min(1,1), point(1)])-offset, ... 
                  max([obj.pData.global_min(1,1), point(1)])+offset, ... 
                  min([obj.pData.global_min(1,2), point(2)])-offset, ...
                  max([obj.pData.global_min(1,2), point(2)])+offset]);
            drawnow
            writeVideo(obj.pData.V_GBM, getframe); 
        end
        % ----------------------------------------------------
        
        % Getter & Setters
        function set.pcaInds(obj, val)
            obj.pcaInds = val;
        end        
        function set.chroms(obj, val)
            obj.chroms = val;
        end
        function set.pData(obj, val)
            obj.pData = val;
        end
        function set.pcaData(obj, val)
            obj.pcaData = val;
        end
        function set.gbmData(obj, val)
            obj.gbmData = val;
        end
        function set.bestGAIndividual(obj, val)
            obj.bestGAIndividual = val;
        end
        
    end
        
end

