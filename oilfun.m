function [rho,cp,mu,nu,k,alpha,Pr,beta] = oilfun(T)
%===============================BEGIN-HEADER============================
% FILE: oilfun.m
% AUTHOR: Spencer Talbot
% DATE: 10/29/19
% 
% PURPOSE: Return engine oil (unused) properties (e.g. rho,cp,mu,nu,k,
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
% V1 - Oil1
% V2 - 
% V3 - 
% 
%==========================================END-HEADER======================

% Inputting table values into MatLab
T_oil = [273 280:10:430]';
rho = [899.1 895.3 890 884.1 877.9 871.8 865.8 859.9 853.9 847.8 841.8 836 830.6 825.1 818.9 812.1 806.5]';
cp = [1.796 1.827 1.868 1.909 1.951 1.993 2.035 2.076 2.118 2.161 2.206 2.250 2.294 2.337 2.381 2.427 2.471]'*10^3; 
mu = [385 217 99.9 48.6 25.3 14.1 8.36 5.31 3.56 2.52 1.86 1.41 1.10 .874 .698 .564 .470]'*10^(-2);
nu = [4280 2430 1120 550 288 161 96.6 61.7 41.7 29.7 22 16.9 13.3 10.6 8.52 6.94 5.83]'*10^(-6);
k = [147 144 145 145 145 143 141 139 138 138 137 136 135 134 133 133 132]'*10^(-3);
alpha = [.910 .88 .872 .859 .847 .823 .8 .779 .763 .753 .738 .723 .709 .695 .682 .675 .662]'*10^(-7);
Pr = [47000 27500 12900 6400 3400 1965 1205 793 546 395 300 233 187 152 125 103 88]';
beta = ones(17,1)*.7;

% prop = [T_oil rho cp mu nu k alpha Pr beta]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_oil)
    if T == T_oil(i)
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
for i = 2:(length(T_oil))
    if T > T_oil(i-1) && T < T_oil(i)
        T1 = T_oil(i-1);
        T2 = T_oil(i);

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
if T < T_oil(1) || T > T_oil(end)
    error('Input temperature exceeds table range')
end
end