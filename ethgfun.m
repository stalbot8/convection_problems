function [rho,cp,mu,nu,k,alpha,Pr,beta] = ethgfun(T)
%===============================BEGIN-HEADER============================
% FILE: ethgfun.m
% AUTHOR: Spencer Talbot
% DATE: 10/29/19
% 
% PURPOSE: Return ethylene glycol properties (e.g. rho,cp,mu,nu,k,
% alpha,Pr,beta) for any given a temperature within the bounds of Table 
% A.5. This function with give interpolated properties as well. I hope you 
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
% function in order for this function to run. 
%
%
% VERSION HISTORY
% V1 - Ethg1
% V2 - 
% V3 - 
% 
%==========================================END-HEADER======================

% Inputting table values into MatLab
T_ethg = [273 280:10:370 373]';
rho = [1130.8 1125.8 1118.8 1114.4 1103.7 1096.2 1089.5 1083.8 1079 1074 1066 1058.5]';
cp = [2.294 2.323 2.368 2.415 2.460 2.505 2.549 2.592 2.637 2.682 2.728 2.742]'*10^3;
mu = [6.51 4.2 2.47 1.57 1.07 .757 .561 .431 .342 .278 .228 .215]'*10^(-2);
nu = [57.6 37.3 22.1 14.1 9.65 6.91 5.15 3.98 3.17 2.59 2.14 2.03]'*10^(-6);
k = [242 244 248 252 255 258 260 261 261 261 262 263]'*10^(-3);
alpha = [.933 .933 .936 .939 .939 .94 .936 .929 .917 .906 .9 .906]'*10^(-7);
Pr = [617 400 263 151 103 73.5 55 42.8 34.6 28.6 23.7 22.4]';
beta = ones(17,1)*.65;

% prop = [T_ethg rho cp mu nu k alpha Pr beta]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_ethg)
    if T == T_ethg(i)
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
for i = 2:(length(T_ethg))
    if T > T_ethg(i-1) && T < T_ethg(i)
        T1 = T_ethg(i-1);
        T2 = T_ethg(i);

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
if T < T_ethg(1) || T > T_ethg(end)
    error('Input temperature exceeds table range')
end
end