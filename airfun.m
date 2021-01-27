function [rho,cp,mu,nu,k,alpha,Pr] = airfun(T)
%===============================BEGIN-HEADER============================
% FILE: airfun.m
% AUTHOR: Spencer Talbot
% DATE: 10/29/19
% 
% PURPOSE: Return air properties (e.g. rho,cp,mu,nu,k,alpha,Pr) for any
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
T_air = [100:50:1000 1100:100:2500 3000]'; % [Kelvin]
rho = [3.5562 2.3364 1.7458 1.3947 1.1614 .9950 .8711 .7740 .6964 .6329 .5804 .5356 .4975 .4643 .4354 .4097 .3868 .3666 .3482 .3166 .2902 .2679 .2488 .2322 .2177 .2049 .1935 .1833 .1741 .1658 .1582 .1513 .1448 .1389 .1135]';
cp = [1.032 1.012 1.007 1.006 1.007 1.009 1.014 1.021 1.030 1.040 1.051 1.063 1.075 1.087 1.099 1.110 1.121 1.131 1.141 1.159 1.175 1.189 1.207 1.230 1.248 1.267 1.286 1.307 1.337 1.372 1.417 1.478 1.558 1.665 2.726]'*10^3;
mu = [71.1 103.4 132.5 159.6 184.6 208.2 230.1 250.7 270.1 288.4 305.8 322.5 338.8 354.6 369.8 384.3 398.1 411.3 424.4 449.0 473.0 496.0 530 557 584 611 637 663 689 715 740 766 792 818 955]'*10^(-7);
nu = [2 4.426 7.590 11.44 15.89 20.92 26.41 32.39 38.79 45.57 52.69 60.21 68.1 76.37 84.93 93.8 102.9 112.2 121.9 141.8 162.9 185.1 213 240 268 298 329 362 396 431 468 506 547 589 841]'*10^(-6);
k = [9.34 13.8 18.1 22.3 26.3 30 33.8 37.3 40.7 43.9 46.9 49.7 52.4 54.9 57.3 59.6 62.0 64.3 66.7 71.5 76.3 82 91 100 106 113 120 128 137 147 160 175 196 222 486]'*10^(-3);
alpha = [2.54 5.84 10.3 15.9 22.5 29.9 38.3 47.2 56.7 66.7 76.9 87.3 98 109 120 131 143 155 168 195 224 257 303 350 390 435 492 534 589 646 714 783 869 960 1570]'*10^(-6);
Pr = [.786 .758 .737 .720 .707 .7 .69 .686 .684 .683 .685 .690 .695 .702 .709 .716 .72 .723 .726 .728 .728 .719 .703 .685 .688 .685 .683 .677 .672 .667 .655 .647 .63 .613 .536]';

% prop = [T_air rho cp mu nu k alpha Pr]; % This isn't used below, but I thought it may be useful to have everything in one matrix

% If the temperature requires no interpolation, this section of the code
% will find it
for i = 1:length(T_air)
    if T == T_air(i)
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
for i = 2:(length(T_air))
    if T > T_air(i-1) && T < T_air(i)
        T1 = T_air(i-1);
        T2 = T_air(i);

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
if T < T_air(1) || T > T_air(end)
    error('Input temperature exceeds table range')
end
end

