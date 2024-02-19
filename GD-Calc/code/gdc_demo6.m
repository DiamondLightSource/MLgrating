% gdc_demo6: Biperiodic grating - skewed metal grid, unit cell A
%
% Test case from https://doi.org/10.1364/JOSAA.14.002758, Example 3. The
% metal permittivity (stratum_pmt) for this case is stated in JOSA as 1+5i,
% but it is apparently supposed to be (1+5i)^2; see Example 4 and Figure 6
% in https://doi.org/10.1088/1464-4258/5/4/307.
%
% Note: The number of stripes used in the JOSA paper cannot be plotted; see
% N assignment.
%
% gdc_demo6 and gdc_demo7 cover the following topics:
%   - biperiodic grating modeling with different unit cells and stripe
%     orientations
%   - biperiodic grating with non-orthogonal period vectors
%	- metallic grating
%   - diffraction order selection
%   - stripe-induced symmetry error
%   - Compare computation results with published numeric data.
%     (Note: The results do not agree with the published data. The
%     published reflectance efficiencies in the JOSA paper appear
%     unexpectedly high in comparison to a bare metal substrate
%     reflectivity.)
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo6.m')

alt_stripe = 1; % selector for alternate stripe orientation (1 or 2)

% Define parameters for grating structure and incident field:
substrate_pmt = 2.25; % substrate permittivity
stratum_pmt = (1+5i)^2; % metal grid permittivity (See header comment.)
wavelength = 1;
h = wavelength; % grating height
d = 2*wavelength; % grating period
theta = 30; % incidence polar angle from x1 axis, deg
phi = 30; % incidence azimuthal angle around x1 axis (zero on x2 axis), deg
if 0
    % This N setting matches the JOSA reference, but cannot be plotted.
    N = 1000; % number of partition stripes covering each hole
else
    disp('The number of stripes (N) is reduced to allow plotting.')
    N = 25;
end
m_max = 12; % maximum diffraction order index

% Construct grating.
clear grating
grating.pmt = ...
    {1,substrate_pmt,stratum_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
% Define the x2 and x3 projections of the first grating period
% (d21,d31) and second grating period (d22,d32), corresponding to unit
% cell A. The second period is parallel to the x2 axis, and the first
% period is at an angle zeta to the x3 axis. (Note: The second period
% defines the stripe orientation - see Figure 3 in GD-Calc.pdf.)
zeta = 30; % deg
grating.d21 = d*sind(zeta);
grating.d31 = d*cosd(zeta);
grating.d22 = d;
grating.d32 = 0;
% Construct the stratum. (Refer to Figure 3 in GD-Calc.pdf and view the
% grating plot with a small N value to follow the construction logic.)
% The x2, x3 coordinate origin is at the center of a grid hole.
clear stratum
stratum.type = 2;
stratum.thick = h;
if alt_stripe==1
    % Stripes are parallel to [grating.d22,grating.d32].
    stratum.h11 = 1;
    stratum.h12 = 0;
    stratum.h21 = 0;
    stratum.h22 = 1;
else
    % Stripes are parallel to [grating.d21,grating.d31].
    stratum.h11 = 0;
    stratum.h12 = 1;
    stratum.h21 = 1;
    stratum.h22 = 0;
end
stratum.stripe = cell(1,N+1);
clear stripe block
stripe.type = 1;
block.c2 = 0.25;
block.pmt_index = 1;
stripe.block{1} = block;
block.c2 = 0.75;
block.pmt_index = 3;
stripe.block{2} = block;
clear block
for n = 1:N
    stripe.c1 = 0.5*n/N-0.25;
    stratum.stripe{n} = stripe;
end
clear stripe
stripe.type = 0;
stripe.c1 = 0.75;
stripe.pmt_index = 3;
stratum.stripe{N+1} = stripe;
clear stripe
grating.stratum = {stratum};
clear stratum

% Define the indicent field.
clear inc_field
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1 = -cosd(theta)/wavelength.
inc_field.wavelength = wavelength;
inc_field.f2 = sind(theta)*cosd(phi)/wavelength;
inc_field.f3 = sind(theta)*sind(phi)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
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
% grating. (These include evanescent waves and, if the substrate's
% permittivity is not real-valued, all transmitted waves.)
R = R(imag([scat_field.f1r])==0);
T = T(imag([scat_field.f1t])==0);
% Tabulate the diffraction order indices and diffraction efficiencies for
% an incident field polarized parallel to the incident plane.
disp(' ');
disp('Diffraction efficiencies (m1, m2, eff2)');
disp('R:');
disp(num2str([[R.m1].' [R.m2].' [R.eff2].']));
disp('T:');
disp(num2str([[T.m1].' [T.m2].' [T.eff2].']));

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.5,.5,.5];
pmt_display(2).alpha = 1;
pmt_display(3).name = '';
pmt_display(3).color = [.75,.75,.75];
pmt_display(3).alpha = 1;
x_limit = [-0.5*h,-1.5*d,-1.5*d;1.5*h,1.5*d,1.5*d];
gdc_plot(grating,1,pmt_display,x_limit);
view(-45,45)
