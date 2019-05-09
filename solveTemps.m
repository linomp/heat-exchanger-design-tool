% A priori steady state temps calculation (epsilon-NTU method) (hacky but correct)
% credits: https://www.mathworks.com/matlabcentral/fileexchange/46303-heat-exchanger-solver
 
% program to test Epsilon - NTU method.
% assuming a counterflow heat exhanger with known properties, known inlet
% temperatures and mass flowrates.


function [T_hot_out, T_cold_out] = solveTemps(data, U, tube_diam, L) 
    
    A = pi*tube_diam*L;  % m^2

    if (data.T_tube_in > data.T_shell_in)         
        c_p_hot = data.cp_tube; % KJ/kg-K
        m_dot_hot = data.mass_flow_tube;
        T_hot_in = data.T_tube_in; % K

        c_p_cold = data.cp_shell; % KJ/kg-K 
        m_dot_cold = data.mass_flow_shell;
        T_cold_in = data.T_shell_in; 
    else 
        c_p_hot = data.cp_shell; % KJ/kg-K 
        m_dot_hot = data.mass_flow_shell;
        T_hot_in = data.T_shell_in;  
        
        c_p_cold = data.cp_tube; % KJ/kg-K
        m_dot_cold = data.mass_flow_tube;
        T_cold_in = data.T_tube_in; % K   
    end
          
    [T_hot_out, T_cold_out] = HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,'Counter Flow');
        
    
    %Q_hot=m_dot_hot*c_p_hot*(T_hot_in-T_hot_out)
    %Q_cold=m_dot_cold*c_p_cold*(T_cold_in-T_cold_out)

    %LMTD = ((T_hot_in-T_cold_out)-(T_hot_out-T_cold_in))/log((T_hot_in-T_cold_out)/(T_hot_out-T_cold_in));
    %Q_Log_Mean_Temp_Diff = U*A*LMTD


function [T_hot_out,T_cold_out]=HeatExchanger(m_dot_hot,c_p_hot,T_hot_in,m_dot_cold,c_p_cold,T_cold_in,U,A,HE_Type);
% [T_hot_out,T_cold_out]=HeatExchanger(c_p_hot,m_dot_hot,T_hot_in,c_p_cold,m_dot_cold,T_cold_in,U,A,HE_Type);
% This function calculates the outlet temperatures of a heat exchanger
% using Epsilon-NTU method. This function uses effectiveness.m as a
% function and should have access to that function.
% 
% The inputs are as follows:
% Hot Flow: c_p_hot, m_dot_hot, T_hot_in.
% Cold Flow: c_p_cold, m_dot_cold, T_cold_in.
% Heat exchanger design parameters: U,A, HE_Type.
% 
% HE_Type defines the type of heat exchanger: (see reference)
%   'Parallel Flow'
%   'Counter Flow'
%   'One Shell Pass'
%   'N Shell Pass'
%   'Cross Both Unmixed'
%   'Cross Cmax Mixed'
%   'Cross Cmin Mixed'
%
% Reference:
% Frank P. Incropera, Introduction to heat transfer. New York:Wiley, 1985, Section 11.4.
% Programmer: Seyyed Ali Hedayat Mofidi (seyyed4li@yahoo.com)

C_hot = m_dot_hot*c_p_hot;
C_cold = m_dot_cold*c_p_cold;
C_min = min(C_hot,C_cold); % finds the flow with lower heat capacity and higher temperature change.
C_max = max(C_hot,C_cold); % finds the flow with higher heat capacity and lower temperature change. 
C_r=C_min/C_max;
NTU = U*A/C_min;
epsilon = effectiveness (NTU,C_r,HE_Type);
Q_max = C_min*(T_hot_in-T_cold_in);
Q = epsilon * Q_max ;
T_hot_out = T_hot_in - Q/C_hot ;
T_cold_out = T_cold_in + Q/C_cold ;


function epsilon=effectiveness(NTU,C_r,HE_Type,N)
% This function calculates the effectiveness of a heat exchanger.
% NTU is the number of transfer units of the Heat Exchanger:
% NTU = U*A / C_min
% epsilon = f (NTU,C_r)
% where C_min is the total heat capacity of the hot or cold flow, whichever
% is smaller and C_max is total heat capacity of other flow:
% C_min = min(C_hot,C_cold); 
% C_max = max(C_hot,C_cold); 
% C_r is the ratio of C_min and C_max:
% C_r=C_min/C_max;
%
% Regardless of heat exchanger type, if C_r=0, either hot flow is
% condensing ( means no change in T_hot) or cold flow is evaporating ( no
% change in T_cold), therefore C_max =inf its temperature does not change. 
%
% HE_Type defines the type of heat exchanger: (see reference)
%   'Parallel Flow'
%   'Counter Flow'
%   'One Shell Pass'
%   'N Shell Pass'
%   'Cross Both Unmixed'
%   'Cross Cmax Mixed'
%   'Cross Cmin Mixed'
%
% Reference:
% Frank P. Incropera, Introduction to heat transfer. New York:Wiley, 1985, Section 11.4. 
% Programmer: Seyyed Ali Hedayat Mofidi (seyyed4li@yahoo.com)

if nargin == 3
    N=1;
end
%% ===== Calculating effectiveness =====
% Special case of boiling or condensing:
if C_r == 0
    epsilon = 1-exp(-NTU);
    return;
end

switch HE_Type
    case 'Parallel Flow'
        epsilon = (1-exp(-NTU*(1+C_r)))/(1+C_r);
    case 'Counter Flow'
        if C_r==1
            epsilon = NTU/(1+NTU);
        else
            epsilon = (1-exp(-NTU*(1-C_r)))/(1-C_r*exp(-NTU*(1-C_r)));
        end
    case 'One Shell Pass'
        epsilon = 2/(1+C_r+sqrt(1+C_r^2)*(1+exp(-NTU*sqrt(1+C_r^2)))/(1+exp(-NTU*sqrt(1+C_r^2))));
    case 'N Shell Pass'
        NTUN = NTU/N;
        epsilon1 = 2/(1+C_r+sqrt(1+C_r^2)*(1+exp(-NTUN*sqrt(1+C_r^2)))/(1+exp(-NTUN*sqrt(1+C_r^2))));
        epsilon = (((1-epsilon1*C_r)/(1-epsilon1))^N-1) / (((1-epsilon1*C_r)/(1-epsilon1))^N-C_r);
    case 'Cross Both Unmixed'
        epsilon = 1-exp(1/C_r * NTU^0.22 * (exp(-C_r*NTU^0.78)-1));
    case 'Cross Cmax Mixed'
        epsilon = 1/C_r*(1-exp(-C_r*(1-exp(-NTU))));
    case 'Cross Cmin Mixed'
        epsilon = 1 - epx(-1/C_r*(1-exp(-C_r*NTU)));
    otherwise % the type is not in the list, therefore we assume there's no heat exchanger.
        epsilon = 0; 
end



    

    
    
