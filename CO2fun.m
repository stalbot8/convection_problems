function [rho,cp,mu,nu,k,alpha,Pr] = CO2fun(T)
%===============================BEGIN-HEADER============================
% FILE: CO2fun.m
% AUTHOR: Spencer Talbot
% DATE: 10/29/19
% 
% PURPOSE: Return CO2 properties (e.g. rho,cp,mu,nu,k,alpha,Pr) for any
% given a temperature within the bounds of Table A.5. This function with 
% give interpolated properties as well. I hope you can find this code
% as useful as I have in saving time looking up table values.
%
% INPUT: Temperature [Kelvin](usually T_m, T_f, or T_mbar--depends on the
% situation)
%
%
% OUTPUT: rho (density [kg/m^3]), cp (specific heat [kJ/(kg*K)]), mu
% (dynamic viscosity [N*s/m^2]), nu (kinematic viscosity [m^2/s]), k
% (fluid thermal conductivity [W/(m*K)]), alpha (thermal diffusivity
% [m^2/s]), Pr (Prandtl Number [Nondimensional])
%
%
% NOTES: You need to have the intrp() function in the same folder as this
% function in order for this function to run. 
%
%
% VERSION HISTORY
% V1 - Air1
% V2 - 
% V3 - 
% 
%==========================================END-HEADER======================

% Inputting table values into MatLab
T_CO2 = [280:20:400 450:50:800]'; % [Kelvin]
rho = [1.9022 1.7730 1.6609 1.5618 1.4743 1.3961 1.3257 1.1782 1.0594 .9625 .8826 .8143 .7564 .7057 .6614]';
cp = [.83 .851 .872 .891 .908 .926 .942 .981 1.02 1.05 1.08 1.1 1.13 1.15 1.17]'*10^3;
mu = [140 149 156 165 173 181 190 210 231 251 270 288 305 321 337]'*10^(-7);
nu = [7.36 8.4 9.39 10.6 11.7 13 14.3 17.8 21.8 26.1 30.6 35.4 40.3 45.5 51]'*10^(-6);
k = [15.2 16.55 18.05 19.7 21.2 22.75 24.3 28.3 32.5 36.6 40.7 44.5 48.1 51.7 55.1]'*10^(-3);
alpha = [9.63 11 12.5 14.2 15.8 17.6 19.5 24.5 30.1 36.2 42.7 49.7 56.3 63.7 71.2]'*10^(-6);
Pr = [.765 .766 .754 .746 .741 .737 .737 .728 .725 .721 .717 .712 .717 .714 .716]';

% prop = [T_CO2 rho cp mu nu k alpha Pr]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_CO2)
    if T == T_CO2(i)
        rho = rho(i);
        cp = cp(i);
        mu = mu(i);
        nu = nu(i);
        k = k(i);
        alpha = alpha(i);
        Pr = Pr(i);
    end
end

% Interpolation occurs here if needed to find properties
for i = 2:(length(T_CO2))
    if T > T_CO2(i-1) && T < T_CO2(i)
        T1 = T_CO2(i-1);
        T2 = T_CO2(i);

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
    end
end

% Error will appear if temperature exceeds table temperature bounds
if T < T_CO2(1) || T > T_CO2(end)
    error('Input temperature exceeds table range')
end
end
