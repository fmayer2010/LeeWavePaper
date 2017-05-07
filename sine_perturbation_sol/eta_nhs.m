function [ eta ] = eta_nhs( K,N,U,h_hat,xx,zz,order )
%returns the dimensionless vertical perturbation of a non-hydrostatic lee-wave
%streamline (solving Long's non-hydrostatic model) from 
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

    M = sqrt(Kc^2-K^2); 
    eta0 = cos(K.*xx + M.*zz);

    eta1 = 1/2.*sin(2*K.*xx + M.*zz) - 1/2.*sin(2*K.*xx + 2*M.*zz);

    eta2 = 1/2.*cos(K.*xx + M.*zz)+...
        (1/2.* cos(K.*xx + M.*zz).*cos(2.*K.*xx + M.*zz) ...
        - 1/2.*sin(K.*xx + M.*zz).*sin(2.*K.*xx + M.*zz) ...
        + 1/2.*sin(K.*xx + M.*zz).*sin(2.*K.*xx + 2.*M.*zz) ...
        -cos(K.*xx + M.*zz).^3);
else %evanescent solution
    M = sqrt(K^2-Kc^2);
    eta0 = cos(K.*xx).*exp(-M.*zz);
    
    eta1 = 1/2.*sin(2*K.*xx).*exp(-M.*zz) - 1/2.*sin(2*K.*xx).*exp(-2*M.*zz);

    eta2 = 1/2.*cos(K.*xx).*exp(-M.*zz)+...
        (1/2.* cos(K.*xx).*exp(-M.*zz).*cos(2.*K.*xx).*exp(-M.*zz) ...
        - 1/2.*sin(K.*xx).*exp(-M.*zz).*sin(2.*K.*xx).*exp(-M.*zz) ...
        + 1/2.*sin(K.*xx).*exp(-M.*zz).*sin(2.*K.*xx).*exp(-2.*M.*zz) ...
        -cos(K.*xx).*exp(-M.*zz).^3);
end

if order==0
    eta = eta0;
elseif order==1
    eta = eta0 + J.*eta1;
elseif order==2
    eta = eta0 + J.*eta1 + J^2.*eta2;
else
    error('order must be 0, 1, or 2')
end





end
