function [rho,cp,mu,nu,k,alpha,Pr,beta] = vaporfun(T)
%===============================BEGIN-HEADER============================
% FILE: vaporfun.m
% AUTHOR: Spencer Talbot
% DATE: 10/29/19
% 
% PURPOSE: Return water vapor (steam) properties (e.g. rho,cp,mu,nu,k,
% alpha,Pr,beta) for any given a temperature within the bounds of Table 
% A.6. This function with give interpolated properties as well. I hope you 
% can find this code as useful as I have in saving time looking up table 
% values.
%
% INPUT: Temperature [Kelvin](usually T_m, T_f, or T_mbar--depends on the
% situation)
%
%
% OUTPUT: rho (density [kg/m^3]), cp (specific heat [kJ/(kg*K)]), mu
% (dynamic viscosity [N*s/m^2]), nu (kinematic viscosity [m^2/s]), k
% (fluid thermal conductivity [W/(m*K)]), alpha (thermal diffusivity
% [m^2/s]), Pr (Prandtl Number [Nondimensional]), beta (coefficient of
% volumetric thermal expansion [K^-1])
%
%
% NOTES: You need to have the intrp() function in the same folder as this
% function in order for this function to run. This function only includes
% the fluid values of the water table. Also note that some of the table
% values for higher temperatures do not include beta. Also, this table
% technically assumes that the water is a saturated vapor, but this 
% function does not return the corresponding pressure to the temperature value.
%
%
% VERSION HISTORY
% V1 - Vapor1
% V2 - 
% V3 - 
% 
%==========================================END-HEADER======================

T_vapor = [380	400	450	500	550	600	650	700	750	800	850]';
rho = 1./([0.58630	0.55420	0.49020	0.44050	0.40050	0.36520	0.33800	0.31400	0.29310	0.27390	0.25790]');
cp = [2.06000	2.01400	1.98000	1.98500	1.99700	2.02600	2.05600	2.08500	2.11900	2.15200	2.18600]'*10^3;
mu = [127.10000	134.40000	152.50000	170.40000	188.40000	206.70000	224.70000	242.60000	260.40000	278.60000	296.90000]'*10^(-7);
nu = [21.68000	24.25000	31.11000	38.68000	47.04000	56.60000	66.48000	77.26000	88.84000	101.70000	115.10000]'*10^(-6);
k = [24.60000	26.10000	29.90000	33.90000	37.90000	42.20000	46.40000	50.50000	54.90000	59.20000	63.70000]'*10^(-3);
alpha = [20.40000	23.40000	30.80000	38.80000	47.40000	57.00000	66.80000	77.10000	88.40000	100.00000	113.00000]'*10^(-6);
Pr = [1.06000	1.04000	1.01000	0.99800	0.99300	0.99300	0.99600	1.00000	1.00000	1.01000	1.02000]';
beta = 1/T;

% prop = [T_vapor rho cp mu nu k alpha Pr beta]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_vapor)
    if T == T_vapor(i)
        rho = rho(i);
        cp = cp(i);
        mu = mu(i);
        nu = nu(i);
        k = k(i);
        alpha = alpha(i);
        Pr = Pr(i);
        beta = beta(i);
    end
end

% Interpolation occurs here if needed to find properties
for i = 2:(length(T_vapor))
    if T > T_vapor(i-1) && T < T_vapor(i)
        T1 = T_vapor(i-1);
        T2 = T_vapor(i);

        rho1 = rho(i-1,:);
        rho2 = rho(i,:);    
        rho = intrp(T1,T,T2,rho1,rho2);

        cp1 = cp(i-1,:);
        cp2 = cp(i,:);
        cp = intrp(T1,T,T2,cp1,cp2);

        mu1 = mu(i-1,:);
        mu2 = mu(i,:);
        mu = intrp(T1,T,T2,mu1,mu2);

        nu1 = nu(i-1,:);
        nu2 = nu(i,:);
        nu = intrp(T1,T,T2,nu1,nu2);

        k1 = k(i-1,:);
        k2 = k(i,:);
        k = intrp(T1,T,T2,k1,k2);

        alpha1 = alpha(i-1,:);
        alpha2 = alpha(i,:);
        alpha = intrp(T1,T,T2,alpha1,alpha2);

        Pr1 = Pr(i-1,:);
        Pr2 = Pr(i,:);
        Pr = intrp(T1,T,T2,Pr1,Pr2);
        
        beta1 = beta(i-1,:);
        beta2 = beta(i,:);
        beta = intrp(T1,T,T2,beta1,beta2);
    end
end

% Error will appear if temperature exceeds table temperature bounds
if T < T_vapor(1) || T > T_vapor(end)
    error('Input temperature exceeds table range')
end
end