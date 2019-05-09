function Temp_tube = tube(data, t, x, doPlots) 
 
T_tube_in = data.T_tube_in; 
T_tube_out = data.T_tube_out;  
T_shell_in = data.T_shell_in;
To = T_shell_in; 

m = 0; % (cartesian coordinates)

options = []; 
tube_sol = pdepe(m, @(x,t,u,DuDx) pde_tube(x, t, u, DuDx, data), ...
                    @(x) ic_tube(x, To), ...
                    @(xl,ul,xr,ur,t) bc_tube(xl,ul,xr,ur,t,T_tube_in,T_tube_out), ... 
                x, t, options);
                
Temp_tube = tube_sol(:,:,1); 
if(doPlots)
    for i = 1:size(Temp_tube, 1) 
        idx = find(Temp_tube(i, :) < T_tube_in);        
        T = Temp_tube(i,idx);
        x2 = x(idx);  
        
        plot(x2, T)
        drawnow
    end 
end


function [c,f,s] = pde_tube(x, t, u, DuDx, data) 
    
    U = data.U;
    r_tube = data.r_tube;
    T_tube_in = data.T_tube_in; 
    T_tube_out = data.T_tube_out;
    rho_tube = data.rho_tube; 
    cp_tube = data.cp_tube*1000; % [J/kg]  
    v_tube = data.v_tube;      
    D = data.D;  
    
    c = 1; % constant coefficient for du/dt
    f = (-v_tube*u) + (D*DuDx); % because df_dx = -v*du_dx
    s = -U*2*pi*(T_tube_in-T_tube_out)*1/(rho_tube*cp_tube*r_tube);
    
    
function uo = ic_tube(x, To)    
    uo = To;
 
    
function [pl,ql,pr,qr] = bc_tube(xl, ul, xr, ur, t, T_tube_in, T_tube_out)  
    ql = 0;
    pl = ul-T_tube_in;
    qr = 0; 
    pr = ur-T_tube_out;  
    
    



