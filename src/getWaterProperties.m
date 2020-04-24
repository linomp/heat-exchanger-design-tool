
function [rho, cp, visc, k] = getWaterProperties(workingTemp)

% load data tables & results from previous interpolation
load('coeffs_water_props.mat') 

% Truncar segun datos disponibles
if (workingTemp > C(end, 1))
    workingTemp = C(end, 1);
elseif (workingTemp < C(1, 1))
    workingTemp = C(1, 1);
end

%Densidad
rho = polyval(coeffs_rho, workingTemp);

% Viscosidad
visc = polyval(coeffs_visc, workingTemp);

% Thermal Conductivity
k = polyval(coeffs_k, workingTemp);

% Cp
cp = polyval(coeffs_cp, workingTemp);

end
