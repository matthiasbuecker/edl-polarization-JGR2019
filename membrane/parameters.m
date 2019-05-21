% Prepare and save parameter structure
% Matthias Buecker, May 2019
clear all; clc;

%% Constants and model paramters
% Fundamental constants
m.eps0 = 8.854e-12;     % Vaccum permitivity (F/m)
m.F = 96485.3;          % Faraday's constant (C/mol)
m.kB = 8.617e-5;        % Boltzmann's constant (eV/K)
m.NA = 6.02214086e23;   % Avogadro's constant (1/mol)
m.e = 1.6e-19;          % Elementary charge (As) 

% Model parameters
m.T = 293;              % Absolute temperature [K]
m.mu = 5e-8;            % Ion mobility in bulk electrolyte (m^2/(Vs))
m.muS = m.mu;           % Ion mobility in Stern layer (m^2/(Vs))
m.epsA = 80;            % Relative permitivity of electrolyte (-)
m.epsI = 4.5;           % Relative permitivity of solid phase (-)
m.c0 = 1;               % Bulk conc. (mol/m^3) anions
m.eta = 1e9;            % Electrolyte viscosity (kg*m/s)
m.R1 = 2e-6;            % Radius of wide pore (m)
m.R2 = 0.2e-6;          % Radius of narrow pore (m)
m.L1 = 90e-6;           % Length of wide pore (m)
m.L2 = 10e-6;           % Length of narrow pore (m)
m.E0 = 1;               % External electrical field (V/m)
m.sigma = 0.01;         % Total surface charge density (C/m^2)
m.p = [0,0.05,0.2:0.2:1];% Partitioning factor (-) 

% Fixed model parameters and derived quantities
m.D = m.mu*m.kB*m.T;        % Ion diffusivity in bulk electrolyte (m^2/s) 
m.DS = m.muS*m.kB*m.T;      % Ion diffusivity in Stern layer (m^2/s) 
m.ka = 2*m.c0*m.mu*m.F;     % Electrolyte conductivity (S/m)
m.ld = sqrt((m.eps0*m.epsA*...      % Debye length (m)
    m.kB*m.T)/(2*m.c0*m.F));

% Model parameters are collected in the structure "m" and saved in order to
% make them available for all post-processing routines
save('./parameters','m')