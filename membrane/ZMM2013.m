function s = ZMM2013(T,c0,mu,epsA,w,L1,R1,L2,R2,f0,fQ,wide)
    % ZMM2013 computes the effective conductivity of a sequence of two
    % cylindrical pores based on the 1D Marshall-Madden impedance (Marshall
    % and Madden, 1959) and described in Buecker and Hoerdt (2013)
    %
    % Matthias Buecker, May 2019
    %
    % Input parameters:
    % T:        Temperature (K)
    % c0:       Bulk ion concentration (mol/m^3)
    % mu:       Ion mobility in bulk electrolyte (m^2/(Vs))
    % epsA:     Relative permitivity of electrolyte (-)
    % w:        Angular frequency (1/s)
    % L1:       Length of pore 1 (m)
    % R1:       Radius of pore 1 (m)
    % L2:       Length of pore 2 (m)
    % R2:       Radius of pore 2 (m)
    % f0:       Zeta potential (V) 
    % fQ:       Partitioning factor (-) 
    % wide:     0 - no surface charge on walls of wide pore
    %           1 - surface charge on walls of wide pore
    
    % CONSTANTS
    % Faraday constant (C/mol)
    F = 96485.3;        

    % MODEL PARAMETERS
    % Bulk conductivity (S/m)
    ka = 2*mu*c0*F;
    % Mean ion concentrations zone 1 (p - cations, n - anions)
    [cp1,cn1]=meanc(R1,c0,epsA,T,f0*wide,fQ*wide);
    % Mean ion concentrations zone 2 (p - cations, n - anions)
    [cp2,cn2]=meanc(R2,c0,epsA,T,f0,fQ);
    % Total cross-sectional areas of cylindrical pores (m^2/pi)
    A1 = R1.^2;
    A2 = R2.^2;
    
    % Normalized mean ion concentrations and mobilities. Normalized to 
    % largest cross-sectional area (passive zone)
    Anorm = max(A1,A2);
    % (Effective) cation mobility in zone 1 (m^2/V/s)
    mp1=mu*cp1/c0.*A1./Anorm;
    % (Effective) anion mobility in zone 1 (m^2/V/s)
    mn1=mu*cn1/c0.*A1./Anorm;
    % (Effective) cation mobility in zone 2 (m^2/V/s)
    mp2=mu*cp2/c0.*A2./Anorm;
    % (Effective) anion mobility in zone 2 (m^2/V/s)
    mn2=mu*cn2/c0.*A2./Anorm;

    % Magnitude and phase of the Marshall-Madden impedance (Marshall and 
    % Madden, 1959)
    [absZ,phiZ]=ZMM(T,c0,w,L1,mp1,mn1,L2,mp2,mn2);
    
    % Convert 1D impedance to normalized effective conductivity (S/m)
    Zm = absZ.*(cos(phiZ)+1i*sin(phiZ));
    s = (L1+L2)./(Zm.*ka);