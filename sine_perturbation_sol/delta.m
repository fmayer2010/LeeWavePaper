function [ delta ] = delta( K,N,U,h_hat,xx,zz,order )
%returns the dimensionless vertical perturbation of a streamline from 
%its upstream location using the specified order accurate in J=Nh/U perturbation 
%expansion of the bottom boundary condition eta(xx,0)=hcos(kx). This
%function is written such that it can be called in a fourier synthesis to
%produce streamfunctions for arbitrary topography

% the input variables are as follows:
% K is the wavelength of the sine function
% N is the constant upstream buoyancy frequency
% U is the constant upstream horizontal velocity
% h_hat is the Kth component of the fourier tranform of the bathymetry, in
% units of L^2
% xx and zz are matrices of the grid, produced via [xx,zz]=meshgrid(x,z)
% order is the order of accuracy of the solution. Can be 0, 1, or 2.

H = h_hat*K;
Kc = N/U;
J = N*H/U;

if K<=Kc %propagating solution

    delta0 = cos(K.*xx + Kc.*zz);

    delta1 = 1/2.*sin(2*K.*xx + Kc.*zz);

    delta2 = 1/2.*cos(K.*xx + Kc.*zz);
else %evanescent solution (currently only solved to 0th order
    delta0 = cos(K.*xx).*exp(-Kc.*zz);
    delta1 = 0;
    delta2 = 0;
end

if order==0
    delta = delta0;
elseif order==1
    delta = delta0 + J.*delta1;
elseif order==2
    delta = delta0 + J.*delta1 + J^2.*delta2;
else
    error('order must be 0, 1, or 2')
end





end
