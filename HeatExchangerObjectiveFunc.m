function [cost, extra] = HeatExchangerObjectiveFunc(x, data)
        
    Ds = x(1);
    do = x(2);
    B = x(3);
    n = x(4);%passes
    Nt = x(5);
    
    a1 = 8000;    % $
    a2 = 259.2;   % $/m^2 (of heat exchanger surface) 
    a3 = 0.93;    % dimensionless 
   
    cost = inf;
    extra = [];
    
    % tube side parameters    
    mt = data.mass_flow_tube;    
    Tci = data.T_tube_in - 273;    
    Tco = data.T_tube_out - 273; % WHAT THE FUCK   
    rhot = data.rho_tube;    
    mut = data.tube_viscosity; %viscosity
    muwt = .75*mut; 
    cpt = data.cp_tube;
    kt = data.k_tube;
    Rft = 0.00061; % fouling factor, m^2*K/W
    st = 1.25*x(2); % tube pitch
    vt = (mt/(rhot*0.8*x(2).^2*pi/4))*n/Nt; %velocity of fluid on tube side
    di = 0.8*x(2); % assuming internal diam. is 20% smaller than outer
    Ret = rhot*vt*0.8*x(2)/mut; %Reynold's number tube side
    ft = 0.079/(Ret^0.25); %Darcy friction factor
    Prt = mut*cpt/kt; %Prandtl number tube side

    % shell side parameters 
    ms =  data.mass_flow_shell;
    Thi = data.T_shell_in - 273;
    Tho = data.T_shell_out - 273;  
    rhos = data.rho_shell;
    cps = data.cp_shell;
    mus = data.shell_viscosity;
    muws = .75*mus; 
    ks = data.k_shell;
    Rfs = 0.00061; % fouling factor, m^2*K/W
    
    % energy
    ce = data.energyCost; % $/kWh energy cost
    H = data.annualOperatingTimeInHours;  %annual operating time in hours
    eff = data.pumpEfficiency; % Pump effficiency    
        
    de = 4*(0.43*st.^2 - (0.5*pi*x(2).^2/4))/(0.5*pi*x(2)); %equivalent dia
    As = x(1)*x(3)*(1-x(2))/(1.25*x(2));  %shell side cross section area
    vs = ms/(rhos*As);   %velocity of fluid on shell side
    Res = ms*de/(As*mus);  %shell side Reynold's number
    Prs = mus*cps/ks;  %shell side Prandtl's number
    L = 4.5;
    
    %Thermal Calculations   
    hs = 0.36*(kt/de)*(Res^0.55)*(Prs^(1/3))*(mus/muws)^0.14; %shell side heat transfer coefficient

    if Ret < 2300       
        ht = kt/di*(3.657 + (0.0677*(Ret*Prt*((di/L).^1.33).^(1/3)))/(1+0.1*Prt*(Ret*(di/L)).^0.3));       
    elseif   Ret > 2300 && Ret < 10000       
        ht = kt/di*((1+di/L).^0.67*(ft/8)*(Ret -1000)*Prt/(1+12.7*((ft/8).^0.5)*((Prt.^(2/3))-1)));       
    elseif   Ret > 10000          
        ht = 0.027*kt/do*(Ret.^0.8)*(Prt.^(1/3)).*(mut/muwt).^0.14;   
    end

    U = 1/((1/hs)+Rfs+(x(2)/di)*(Rft+(1/ht))); 
     
    R = (Thi - Tho)/(Tco - Tci);  %correction coefficient    
    P = (Tco - Tci)/(Thi - Tci);  %efficiency
    F = sqrt((R.^2+1)/(R.^2-1));  %correction factor    
    LMTD = ((Thi - Tco) - (Tho - Tci))/(log((Thi-Tco)/(Tho-Tci)));    
    Q = ms*cps*(Thi - Tho);   
    A = Q/(U*F*LMTD);        
    L = A/(pi*x(2)*Nt);
    
    %Pressure drops
    
    %tube side pressure drop
    p = 4;  
    pt = 0.5*rhot*vt.^2*(L*ft/(0.8*x(2))+p)*n;
    
    %shell side pressure drop
    b0 = 0.72;                                            
    fs = 2*b0*Res;    
    ps = fs*(rhos*vs.^2/2)*(L/x(2))*(x(1)/de);    

    % FINAL COST FUNCTION CALCULATION
    cost = (a1 + a2*A^a3) + (ce*H*(mt*pt/rhot + ms*ps/rhos)/eff);  
    
    extra.U = U;
    extra.L = L;
    extra.v_tube = vt;
    extra.v_shell = vs;
    extra.Q = Q;
    extra.T_shell_out = Tho;
    extra.T_tube_out = Tco;
    extra.A = A;
    
end