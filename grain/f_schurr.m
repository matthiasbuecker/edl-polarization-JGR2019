function f = f_schurr(omega,c0,mu,muS,epsilon_a,epsilon_i,R,zeta,sigmaS,T)
% F_SCHURR computes the complex frequency-dependent reflection coefficient
% of the Schurr model (Schurr, 1964).
%
% Input paramters
% omega:            Angular frequency [rad/s]
% c0:               Bulk ion concentration [mol/m^3]
% mu:               Ion mobility in bulk electrolyte [m^2/V/s]
% muS:              Ion mobility in Stern layer [m^2/V/s]
% epsilon_a:     	Relative permitivity of electrolyte [-]
% epsilon_i:    	Relative permitivity of particle [-]
% R:                Particle radius [m]
% zeta:             Zeta potential [V]
% sigmaS:           Surface charge density of Stern layer [C/m^2]
% T:                Absolute temperature [K]

% Matthias Buecker, May 2019

% Fundamental constants 
epsilon_0 = 8.854e-12;      % Vaccum permitivity [F/m]
kB = 8.617e-5;           	% Boltzmann's constant [eV/K]
F = 96485.3;                % Faraday's constant (C/mol)

% Calculated parameters
% Uniform conductivity of electrolyte [S/m]
Ka = 2*c0*mu*F+1i*omega*epsilon_0*epsilon_a;
% Inverse Debye's screening length [1/m] [eq. (2)]
kappa = sqrt(2*c0*F/(epsilon_0*epsilon_a*kB*T));
% Surface charge densities in the diffuse layer [C/m^2]
sigmadp = 2*c0*F/kappa*(exp(-zeta/(2*kB*T))-1);
sigmadm = -2*c0*F/kappa*(exp(zeta/(2*kB*T))-1);
% Surface conductivity of diffuse charges [S]
lambda = abs(sigmadp-sigmadm)*mu;
% Surface conductivity of bound charges (Stern layer) [S] (eq. 9)
lambda0 = abs(sigmaS)*muS; 
% Diffusion coefficient [m^2/s]
DS = muS*kB*T;
% Time constant [s] (eq. 7)
tau = R.^2./(2*DS);
% Surface conductivity
% Dielectric conductivity of the bound charge diffusion process (eq. 11)
epsilon_s = tau.*lambda0./(1+1i*omega.*tau)./epsilon_0;

% Equivalent conductivity of the sphere (eq. 15)
kip = 2*lambda/R;
% Equivalent permittivity of the sphere (eq. 16)
epsilon_ip = epsilon_i+2*epsilon_s/R;
% Effective complex conductivity of the sphere (eq. 14)
Kip = kip + 1i*omega*epsilon_0.*epsilon_ip;

% Dimensionless frequency-dependent reflection coefficient (from eq. 13)
f = (Kip-Ka)./(2*Ka+Kip);