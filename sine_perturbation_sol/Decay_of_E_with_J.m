% This short script plots the nondimensional energy flux over a sine wave
% bottom as computed with a second order accurate in J=Nh/U assymptotic 
% solution to Long's Equation for the perturbation streamlines of a 
% hydrostatic leewave. The computation of E was performed in MAPLE by
% solving for eta to order J^2 accuracey, and then u and w (each to order 
% J^3 accuracy, since derivatives of eta), and then P (to order J^3 
% accuracy, since involves the integral of eta^2), and then multiplying p 
% and w to get E(x,z) (to order J^3 accuracy), and then integrating along 
% the bottom to get E0, which is the energery fluxed from the bottom (per 
% unit length spanwise per wavelength, units of Watts per length) (to J^3
% accuracy).
% E0 is what is plotted here, to 1st, 3rd, and 5th order accuracy. The 5th
% order result is probably not totally correct. But the 3rd order result
% should be.

J = 0:0.01:1;


E5=-(1/8).*(J.^4+2*J.^2-8);
E3=-(1/8).*(2*J.^2-8);
E1=-(1/8).*(-8).*ones(size(J));

figure(1)
plot(J,E1,J,E3,J,E5)
axis([0 1 0 1.1])
xlabel('Juice','Interpreter','latex')
ylabel('$E/E_{linear}$','Interpreter','latex')
leg=legend('O($J$)','O($J^3$)','O($J^5$)','Location','sw');
set(leg,'Interpreter','latex')
title({'Energy flux from h(x)=Asin(kx) as a function of Juice=NA/U',...
    'Solutions via perturbation expansions of Longs equation'},...
    'Interpreter','latex')

print('decay_of_E_with_J', '-depsc');

Ri=J.^-2 - 3/2;
Ri_high=(-6.*J.^2+4)./(J.^2.*(169.*J.^4-52.*J.^2+4));
figure(2)
plot(J,Ri,J,Ri_high,J,0.25.*ones(size(J)),'--k',J,zeros(size(J)),'k')
ylabel('Ri','Interpreter','latex')
xlabel('J','Interpreter','latex')
title({'Maximum gradient Richardson number in the lee wave above a sine bottom, h(x)=Asin(kx)',...
    'as a function of Juice=NA/U, Solutions via perturnation expansion of Longs equation' },...
    'Interpreter','latex')
leg=legend('O($J^2$)','O($J^6$)','Location','sw');
set(leg,'Interpreter','latex')
axis([0.5,1,-0.5,2])



