% gdc_demo16: coherent beam combiner
%
% gdc_demo16 covers the following topic:
%	- multiple incident orders
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% bal(t) is a measure of the balance between the grating's 0 and 1 orders,
% as a function of grating thickness t. bal(t) is zero when the orders are
% in balance.
bal = @(t) diff(abs(calc(0,[0,1],t)));
% Find t for which bal(t)==0. First bracket the zero.
disp(['[bal(0.5),bal(0.6)]: ' ...
    num2str([bal(0.5),bal(0.6)])]) % -0.069204     0.22194
% Optimize t for balanced 0- and 1st-order efficiencies:
t = fzero(bal,[0.5,0.6]);
disp(['[t,bal(t)]: ' num2str([t,bal(t)])]) % 0.5248 -1.0325e-14

% Diffraction efficiencies, orders 0 and 1 (cf. GD-Calc.pdf, Eq's A.29,
% case A = 1, B = 1; with Tps = 0):
[Tss,f1i,f1t] = calc(0,[0,1],t);
eff = abs(Tss).^2.*f1t./f1i;
disp(['Balanced order 0 and 1 efficiencies: ' num2str(eff)]) ...
    % 0.48317     0.48317

% Get the grating's transmission scattering amplitude matrix Tss (TE
% polarization) from incident orders 0 and 1 (Tss rows) to transmitted
% orders 0 and 1 (Tss columns). The incident f1 frequencies in orders 0 and
% 1 are f1i(1) and f1i(2), and the transmitted frequencies in orders 0 and
% 1 are f1t(1) and f1t(2).
[Tss,f1i,f1t] = calc([0;1],[0,1],t);
disp('Tss:')
disp(num2str(Tss))
% 0.88826+0.0017056i     0.0089276-0.88822i  
% 0.0089274-0.88822i     0.88826+0.0017314i

% Define the incident field amplitudes A0i in order 0 and A1i in order 1,
% with phase modulation in order 1.
A0i = 1;
phase = 0:360; % transmitted order 1 phase (deg)
A1i = complex(cosd(phase),sind(phase));
% Calculate the transmitted field amplitudes A0t in order 0 and A1t in
% order 1, via linear superposition.
A0t = A0i*Tss(1,1)+A1i*Tss(2,1);
A1t = A0i*Tss(1,2)+A1i*Tss(2,2);
% Get the incident power Pi, summed over orders 0 and 1 (GD-Calc.pdf, Eq
% A.14 with lambda factor omitted).
Pi = -f1i(1).*abs(A0i).^2-f1i(2).*abs(A1i).^2;
% Get the transmitted power P0t in order 0 and P1t in order 1 (GD-Calc.pdf,
% Eq A.18 with lambda factor omitted).
P0t = -f1t(1).*abs(A0t).^2;
P1t = -f1t(2).*abs(A1t).^2;
% diffraction efficiency eff0 in order 0 and eff1 in order 1:
eff0 = P0t./Pi;
eff1 = P1t./Pi;
figure, plot(phase,eff0)
hold on, plot(phase,eff1)
title('Phase-modulated transmission efficiencies')
xlabel('Phase (deg) of incident order 1')
ylabel('Efficiency')
legend('Order 0','Order 1')

return

function [Tss,f1i,f1t] = calc(m1_inc,m1,t)
% Calculate E-field (TE polarization) transmission scattering matrix Tss
% through the grating from incident orders m1_inc (Tss rows) to transmitted
% orders m1 (Tss columns). m1_inc is a column vector and m1 is a row
% vectors. Also Get f1 frequencies for incident orders (denoted as f1i,
% column vector corresponding to m1_inc) and for tranmitted orders (f1t,
% row vector corresponding to m1). t is the grating thickness.

% Define parameters for grating structure and incident field.
pmt = 2.25; % superstrate and grating permittivity
lambda = 0.633; % wavelength, micron
period = lambda; % grating period
c = 0.5; % line width/period
theta = 30.0; % exit angle (deg)
m_max = 10; % maximum diffraction order index

clear grating
grating.pmt = {1,pmt};
grating.pmt_sub_index = 1;
grating.pmt_sup_index = 2;
grating.d21 = period;
grating.d31 = 0;
grating.d22 = 0;
grating.d32 = period;
clear stratum
stratum.type = 1;
stratum.thick = t;
stratum.h11 = 1;
stratum.h12 = 0;
clear stripe
stripe.c1 = -c/2;
stripe.pmt_index = 1;
stratum.stripe{1} = stripe;
stripe.c1 = c/2;
stripe.pmt_index = 2;
stratum.stripe{2} = stripe;
grating.stratum{1} = stratum;
clear stratum stripe

clear inc_field
inc_field.wavelength = lambda;
inc_field.f2 = -sind(theta)./lambda;
inc_field.f3 = 0;

clear order
order(1).m2 = 0;
order(1).m1 = -m_max:m_max;

[~,scat_field,inc_field] = ...
    gdc(grating,inc_field,order,false,m1_inc.*[1,0]);

m1_inc_ = [scat_field.m1_inc];
m1_ = [scat_field.m1];
Tss = [];
for j = length(m1_inc):-1:1
    for k = length(m1):-1:1
        Tss(j,k) = scat_field(m1_inc_==m1_inc(j) & m1_==m1(k)).Tss;
    end
end
if nargout<2
    return
end
f1t = zeros(size(m1));
for k = 1:length(m1)
    f1t(k) = scat_field(m1_inc_==0 & m1_==m1(k)).f1t;
end
m1_inc_ = [inc_field.m1_inc];
f1i = zeros(size(m1_inc));
for j = 1:length(m1_inc)
    f1i(j) = inc_field(m1_inc_==m1_inc(j)).f1;
end

end % calc

