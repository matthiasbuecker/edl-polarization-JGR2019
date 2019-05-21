function zeta = sigma2zeta(c0,epsilon_a,sigma,mod,T)
% sigma2zeta computes the zeta potential from equivalent surface charge
% density of the diffuse layer.
%
% Input paramters
% c0:               Bulk ion concentration [mol/m^3]
% epsilon_a:     	Relative permitivity of electrolyte [-]
% sigma:            Surface charge densities in diffuse layer [C/m^2]
% mod:              'd' - sigma is surface charge density in diffuse layer
%                   'd+' - sigma is surface charge density of counterions 
%                   in diffuse layer
% T:                Absolute temperature [K]

% Matthias Buecker, May 2019

% Fundamental constants
epsilon_0 = 8.854e-12;      % Vaccum permitivity [F/m]
kB = 8.617e-5;           	% Boltzmann's constant [eV/K]
F = 96485.3;                % Faraday's constant (C/mol)

% Calculated parameters
% Inverse Debye's screening length [1/m] [eq. (2)]
kappa = sqrt(2*c0*F/(epsilon_0*epsilon_a*kB*T));
switch mod
    case 'd'
        % Zeta potential [V]
        zeta = 2*kB*T*asinh(-sigma*kappa/(4*F*c0));
    case 'd+'
        % Zeta potential [V]
        zeta = -2*kB*T*log(sigma*kappa/(2*F*c0)+1);
end