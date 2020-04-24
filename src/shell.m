function Temp_shell = shell(data, t, x, doPlots) 
 
T_shell_in = data.T_shell_in; 
T_shell_out = data.T_shell_out; 
To = T_shell_in;  

m = 0;


options = []; 
shell_sol = pdepe(m, @(x,t,u,DuDx) pde_shell(x, t, u, DuDx, data), ...
                    @(x) ic_shell(x, To), ...
                    @(xl,ul,xr,ur,t) bc_shell(xl,ul,xr,ur,t,T_shell_in,T_shell_out), ... 
                 x, t, options);
                
Temp_shell = shell_sol(:,:,1); 

if(doPlots)
    for i = 1:size(Temp_shell, 1)         
        idx = find(Temp_shell(i, :) < 400);        
        T = Temp_shell(i,idx);
        x2 = x(idx);  
        
        plot(x2, T)
        drawnow
    end 
end

function [c,f,s] = pde_shell(x, t, u, DuDx, data)   
     
    U = data.U;
    r_tube = data.r_tube;
    r_shell = data.r_shell; 
    T_shell_in = data.T_shell_in; 
    T_shell_out = data.T_shell_out;
    rho_shell = data.rho_shell; 
    cp_shell = data.cp_shell*1000; % [J/kg]  
    v_shell = data.v_shell; 
    D = data.D;

    c = 1; % constant coefficient for du/dt 
    f = (-v_shell*u) + (D*DuDx); % because df_dx = -v*du_dx
    
    s = U*2*pi*r_tube*(abs(T_shell_in-T_shell_out))*1/(rho_shell*cp_shell*((r_shell^2)-(r_tube^2))); % no sources
    
    
function uo = ic_shell(x, To)    
    uo = To;
 
    
function [pl,ql,pr,qr] = bc_shell(xl, ul, xr, ur, t, T_shell_in, T_shell_out)  
    ql = 0;       
    pl = ul-T_shell_out;
    qr = 0;
    pr = ur-T_shell_in;



