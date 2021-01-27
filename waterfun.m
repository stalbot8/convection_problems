function [rho,cp,mu,nu,k,alpha,Pr,beta] = waterfun(T)
%===============================BEGIN-HEADER============================
% FILE: waterfun.m
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
% technically assumes that the water is a saturated liquid, but this 
% function does not return the corresponding pressure to the temperature value.
%
%
% VERSION HISTORY
% V1 - Waterf1
% V2 - 
% V3 - 
% 
%==========================================END-HEADER======================

T_water = [273.15 275:5:370 373.15 375:5:390 400:10:620 625:5:645 647.3]';
rho = 1./([1 1 1 1 1.001 1.002 1.003 1.005 1.007 1.009 1.011 1.013 1.016 1.018 1.021 1.024 1.027 1.03 1.034 1.038 1.041 1.044 1.045 1.049 1.053 1.058 1.067 1.077 1.088 1.099 1.110 1.123 1.137 1.152 1.167 1.184 1.203 1.222 1.244 1.268 1.294 1.323 1.355 1.392 1.433 1.482 1.541 1.612 1.705 1.778 1.856 1.935 2.075 2.351 3.170]'*10^(-3));
cp = [4.217 4.211 4.198 4.189 4.184 4.181 4.179 4.178 4.178 4.179 4.180 4.182 4.184 4.186 4.188 4.191 4.195 4.199 4.203 4.209 4.214 4.217 4.220 4.226 4.232 4.239 4.256 4.278 4.302 4.331 4.36 4.40 4.44 4.48 4.53 4.59 4.66 4.74 4.84 4.95 5.08 5.24 5.43 5.68 6 6.41 7 7.85 9.35 10.6 12.6 16.4 26 90 NaN]'*10^3;
mu = [1750 1652 1422 1225 1080 959 855 769 695 631 577 528 489 453 420 389 365 343 324 306 289 279 274 260 248 237 217 200 185 173 162 152 143 136 129 124 118 113 108 104 101 97 94 91 88 84 81 77 72 70 67 64 59 54 45]'*10^(-6);
nu = mu./rho;
k = [569 574 582 590 598 606 613 620 628 634 640 645 650 656 660 664 668 671 674 677 679 680 681 683 685 686 688 688 688 685 682 678 673 667 660 651 642 631 621 608 594 580 563 548 5528 513 497 467 444 430 412 392 367 331 238]'*10^(-3);
alpha = k./(rho.*cp);
Pr = [12.99 12.22 10.26 8.81 7.56 6.62 5.83 5.20 4.62 4.16 3.77 3.42 3.15 2.88 2.66 2.45 2.29 2.14 2.02 1.91 1.8 1.76 1.7 1.61 1.53 1.47 1.34 1.24 1.16 1.09 1.04 0.99 .95 .92 .89 .87 .86 .85 .84 .85 .86 .87 .9 .94 .99 1.05 1.14 1.30 1.52 1.65 2 2.7 4.2 12 NaN]';
beta = [-68.05 -32.74 46.04 114.1 174 227.5 276.1 320.6 361.9 400.4 436.7 471.2 504 535.5 566 595.4 624.2 652.3 697.9 707.1 728.7 750.1 761 788 814 841 896 952 1010 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]'*10^(-6); 

% prop = [T_water rho cp mu nu k alpha Pr beta]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_water)
    if T == T_water(i)
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
for i = 2:(length(T_water))
    if T > T_water(i-1) && T < T_water(i)
        T1 = T_water(i-1);
        T2 = T_water(i);

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
if T < T_water(1) || T > T_water(end)
    error('Input temperature exceeds table range')
end
end