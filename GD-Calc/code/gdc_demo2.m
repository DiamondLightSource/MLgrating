% gdc_demo2: Biperiodic grating comprising rectangular pyramids
%
% Test case from https://doi.org/10.1007/BF00935902. See also
% https://doi.org/10.1364/OL.17.001740;
% https://doi.org/10.1016/0030-4018(93)90546-H.
%
% This demo covers the following topics:
%   - GDC interface and setup for a simple biperiodic grating
%   - Compare computation results with published numeric data.
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo2.m')

% Define parameters for grating structure and incident field:
d1 = 1.50; % grating period in x2 direction
d2 = 1.00; % grating period in x3 direction
h = 0.25; % grating height
wavelength = 1.533;
phi = 30; % incidence polar angle from x1 axis, deg
psi = 45; % incidence azimuthal angle around x1 axis (zero on x2 axis), deg
grating_pmt = 1.5^2; % grating permittivity
L1 = 16; % number of grating strata
m_max = 5; % maximum diffraction order index

% Construct grating.
clear grating
grating.pmt = {1,grating_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
grating.d21 = d1; % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = d2; % second grating period: x3 projection
grating.stratum = {};
clear stratum
stratum.type = 2; % biperiodic stratum
stratum.thick = h/L1; % stratum thickness
% The following h11 ... h22 spec indicates that the stratum's period
% vectors match the grating periods (GD-Calc.pdf, Eq. 3.20).
stratum.h11 = 1;
stratum.h12 = 0;
stratum.h21 = 0;
stratum.h22 = 1;
clear stripe
stripe.type = 0; % first stripe is homogeneous
stripe.pmt_index = 1; % first stripe's permittivity index
stratum.stripe{1} = stripe;
clear stripe
stripe.type = 1; % second stripe is inhomogeneous
stripe.block{1}.pmt_index = 1; % first block's permittivity index
stripe.block{2}.pmt_index = 2; % second block's permittivity index
stratum.stripe{2} = stripe;
clear stripe
for l1 = 1:L1
    % Set first and second stripes' boundaries (c1) on positive side:
    stratum.stripe{1}.c1 = -((L1-l1+0.5)/L1)/2;
    stratum.stripe{2}.c1 = -stratum.stripe{1}.c1;
    % Set first and second block' boundaries (c2) on positive side:
    stratum.stripe{2}.block{1}.c2 = -((L1-l1+0.5)/L1)/2;
    stratum.stripe{2}.block{2}.c2 = -stratum.stripe{2}.block{1}.c2;
    grating.stratum{end+1} = stratum;
end
clear stratum

% Define the indicent field.
clear inc_field
inc_field.wavelength = wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector, Eq's 4.1-2 in
% GD-Calc.pdf. The grating-normal (x1) projection is implicitly
%   f1 = -cosd(phi)/wavelength
inc_field.wavelength = wavelength;
inc_field.f2 = sind(phi)*cosd(psi)/wavelength;
inc_field.f3 = sind(phi)*sind(psi)/wavelength;

% Specify which diffracted orders are to be retained in the calculations,
% Eq's 4.12-13 in GD-Calc.pdf.
order = [];
for m2 = -m_max:m_max
    order(end+1).m2 = m2; %#ok<SAGROW> 
    order(end).m1 = -m_max:m_max;
end

% Run the diffraction calculations.
tic
[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
toc

% Compute the diffraction efficiencies.
[R,T] = gdc_eff(scat_field,inc_field);
% Discard diffracted waves that decay exponentially with distance from the
% grating. These include evanescent waves and, if the substrate's
% permittivity is not real-valued, all transmitted waves.
R = R(imag([scat_field.f1r])==0);
T = T(imag([scat_field.f1t])==0);
% Tabulate the diffraction order indices and diffraction efficiencies for
% an incident field polarized parallel to the incident plane. Also tabulate
% the fractional energy loss in the grating.
disp(' ');
disp('Diffraction efficiencies (m1, m2, eff2)');
disp('R:');
disp(num2str([[R.m1].' [R.m2].' [R.eff2].']));
disp('T:');
disp(num2str([[T.m1].' [T.m2].' [T.eff2].']));
disp('Energy loss:');
disp(num2str(1-sum([[R.eff2].'; [T.eff2].'])));

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75];
pmt_display(2).alpha = 1;
x_limit = [-0.5*h,-1.5*d1,-1.5*d2;1.5*h,1.5*d1,1.5*d2];
h_plot = gdc_plot(grating,1,pmt_display,x_limit);
view(-25,15)
