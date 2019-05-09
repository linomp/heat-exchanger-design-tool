
function [bestParameters]  = mainFunc(data, parent)
    % Optimization routine (Our Method) 

    bestParameters = runOptimization(data);
    data.L = bestParameters.L; %[m]
    data.r_tube = bestParameters.r_tube; %[m]
    data.r_shell = bestParameters.r_shell; %[m]
    data.U = bestParameters.U; %[W/m^2*K]
    data.baffle = bestParameters.baffle;  %[??]
    data.v_tube = bestParameters.v_tube;
    data.v_shell = bestParameters.v_shell; 
    data.T_tube_out = bestParameters.T_tube_out + 273;
    data.T_shell_out = bestParameters.T_shell_out + 273;
    data.Q = bestParameters.Q;


    fprintf('\n\n****Our method****\n')
    bestParameters
    
    fprintf('Total heat transfer: %.2f [kW]\n', data.Q/1000);
    fprintf('Annual Cost: %.2f [USD]\n', bestParameters.cost);
    

    % Kumar et. al Method
    fprintf('\n\n****Kumar et. al Method****\n')
    [kumarCost, kumarResults] = HeatExchangerObjectiveFunc([1.05,0.0459,0.5,1,60, 3.83], data)
    fprintf('Total heat transfer: %.2f [kW]\n', kumarResults.Q/1000);
    fprintf('Annual Cost: %.2f [USD]\n', kumarCost);
    

    % PDE SOLUTION
     
    % Solver parameters
    data.To = data.T_shell_in; %[K] 
    data.D =  1E-3; 
    data.tsim = 200; %[sg]
    t = linspace(0, data.tsim, 150);
    x = linspace(0, data.L, 100);
    dx = x(2)-x(1); 
    Temp_tube = tube (data, t, x, 0);
    Temp_shell = shell (data, t, x, 0); 
    for i = 1:size(Temp_tube, 1)
        hold off     
        plot(parent, x, Temp_tube(i, :), 'linewidth', 2.5)

        hold on     
        plot(parent, x, Temp_shell(i, :), 'linewidth', 2.5)

        xlabel('X [m]');
        ylabel('T [°C]'); 
        text(0, 515, sprintf('t =%8.2f [s]', t(i)), 'fontSize', 12);    
        axis([0 data.L 250 500]);
        axis manual
        drawnow 
    end 
    legend({'Tube', 'Shell'}, 'location', 'best')   
    

end
