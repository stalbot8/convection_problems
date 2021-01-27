function [rho,cp,mu,nu,k,alpha,Pr,beta] = vaporfun(T)
%===============================BEGIN-HEADER============================
% FILE: vaporfun.m
% AUTHOR: Spencer Talbot
% DATE: 10/29/19
% 
% PURPOSE: Return water properties (e.g. rho,cp,mu,nu,k,
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
% V1 - Waterf1
% V2 - 
% V3 - 
% 
%==========================================END-HEADER======================

T_waterg = [273.15 275:5:370 373.15 375:5:390 400:10:620 625:5:645 647.3]';
rho = 1./([206.3 181.7 130.4 99.4 69.7 51.94 39.13 29.74 22.93 17.82 13.98 11.06 8.82 7.09 5.74 4.683 3.846 3.180 2.645 2.212 1.861 1.679 1.574 1.337 1.142 .98 .731 .553 .425 .331 .261 .208 .167 .136 .111 .0922 .0766 .0631 .0525 .0445 .0375 .0317 .0269 .0228 .0193 .0163 .0137 .0115 .0094 .0085 .0075 .0066 .0057 .0045 .0032]');
cp = [1.854 1.855 1.858 1.861 1.864 1.868 1.872 1.877 1.882 1.888 1.895 1.903 1.911 1.920 1.93 1.941 1.954 1.968 1.983 1.999 2.017 2.029 2.036 2.057 2.080 2.104 2.158 2.2221 2.291 2.369 2.46 2.56 2.68 2.79 2.94 3.1 3.27 3.47 3.7 3.96 4.27 4.64 5.09 5.67 6.4 7.35 8.75 11.1 15.4 18.3 22.1 27.6 42 NaN NaN]'*10^3;
mu = [8.02 8.09 8.29 8.49 8.69 8.89 9.09 9.29 9.49 9.69 9.89 10.09 10.29 10.49 10.69 10.89 11.09 11.29 11.49 11.69 11.89 12.02 12.09 12.29 12.49 12.69 13.05 13.42 13.79 14.14 14.5 14.85 15.19 15.54 15.88 16.23 16.59 16.95 17.33 17.72 18.1 18.6 19.1 19.7 20.4 21.5 22.7 24.1 25.9 27 28 30 32 37 45]'*10^(-6);
nu = mu./rho;
k = [18.2 18.3 18.6 18.9 19.3 19.5 19.6 20.1 20.4 20.7 21 21.3 21.7 22 22.3 22.6 23 23.3 23.7 24.1 24.5 24.8 24.9 25.4 25.8 26.3 27.2 28.2 29.8 30.4 31.7 33.1 34.6 36.3 38.1 40.1 42.3 44.7 47.5 50.6 54 58.3 63.7 76.7 76.7 84.1 92.9 103 114 121 130 141 155 178 238]'*10^(-3);
alpha = k./(rho.*cp);
Pr = [.815 .817 .825 .833 .841 .849 .857 .865 .873 .883 .894 .901 .908 .916 .925 .933 .942 .951 .96 .969 .978 .984 .987 .999 1.004 1.013 1.033 1.054 1.054 1.075 1.1 1.12 1.14 1.17 1.2 1.23 1.25 1.28 1.31 1.35 1.39 1.43 1.47 1.52 1.59 1.68 1.84 2.15 2.6 3.46 4.2 4.8 6 9.6 26 NaN]';
beta = [-68.05 -32.74 46.04 114.1 174 227.5 276.1 320.6 361.9 400.4 436.7 471.2 504 535.5 566 595.4 624.2 652.3 697.9 707.1 728.7 750.1 761 788 814 841 896 952 1010 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]'*10^(-6); 

prop = [T_waterg rho cp mu nu k alpha Pr beta]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_waterg)
    if T == T_waterg(i)
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
for i = 2:(length(T_waterg))
    if T > T_waterg(i-1) && T < T_waterg(i)
        T1 = T_waterg(i-1);
        T2 = T_waterg(i);

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
if T < T_waterg(1) || T > T_waterg(end)
    error('Input temperature exceeds table range')
end
end