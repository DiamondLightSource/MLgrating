% gdc_demo5: Biperiodic grating comprising cylindrical pillars
%
% Test case from https://doi.org/10.1364/JOSAA.14.002758, Example 2.
%
% Note: The number of blocks used in the JOSA paper cannot be plotted; see
% N assignment.
%
% gdc_demo5 covers the following topics:
%   - biperiodic grating with non-orthogonal period vectors
%   - approximation of cylindrical structures by block partitioning
%   - choice of unit cell and stripe orientations
%   - harmonic indices
%   - diffraction order selection
%   - stripe-induced symmetry error
%   - Compare computation results with published numeric data.
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo5.m')

alt_stripe = 1; % selector for alternate stripe orientation (1, 2, or 3)

% Define parameters for grating structure and incident field:
grating_pmt = 2.56; % dielectric permittivity
wavelength = 1;
h = wavelength; % grating height
d = 2*wavelength; % grating period
r = 0.5*wavelength; % pillar radius
theta = 30; % incidence polar angle from x1 axis, deg
phi = 30; % incidence azimuthal angle around x1 axis (zero on x2 axis), deg
if 0
    % This N setting matches the JOSA reference, but cannot be plotted.
    N = 1000; % number of partition blocks per semicircle (for pillars)
else
    disp('The number of blocks (N) is reduced to allow plotting.')
    N = 25;
end
m_max = 9; % maximum diffraction order index

% Construct grating.
clear grating
grating.pmt = {1,grating_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
% Define the x2 and x3 projections of the first grating period
% (d21,d31) and second grating period (d22,d32). The second period is
% parallel to the x2 axis, and the first period is at an angle zeta to
% the x3 axis.
zeta = 30; % deg
grating.d21 = d*sind(zeta);
grating.d31 = d*cosd(zeta);
grating.d22 = d;
grating.d32 = 0;
clear stratum
stratum.type = 2; % biperiodic stratum
stratum.thick = h; % stratum thickness
% The following h11 ... h22 spec defines the stratum's period vectors
% (GD-Calc.pdf, Eq 3.20) and determines the stripe orientation.
if alt_stripe==1
    % Stripes are parallel to [grating.d22,grating.d32].
    stratum.h11 = 1;
    stratum.h12 = 0;
    stratum.h21 = 0;
    stratum.h22 = 1;
elseif alt_stripe==2
    % Stripes are parallel to [grating.d21,grating.d31].
    stratum.h11 = 0;
    stratum.h12 = 1;
    stratum.h21 = 1;
    stratum.h22 = 0;
else % alt_stripe==3
    % Stripes are parallel to
    % [grating.d21,grating.d31]-[grating.d22,grating.d32].
    stratum.h11 = 1;
    stratum.h12 = 1;
    stratum.h21 = 0;
    stratum.h22 = -1;
end
stratum.stripe = cell(1,2*N);
% Define vertex coordinates for block-partitioned unit circle (see
% GD-Calc.pdf, Figure 4). The j-th block vertex in the first quadrant
% has coordinates [x(j),x(N+1-j)], j = 1...N. x(j) is monotonic
% decreasing with j.
x = circle_partition(N);
clear stripe block
stripe.type = 1; % inhomogeneous stripe
% Construct the stratum. (Refer to Figures 3 and 4 in GD-Calc.pdf and
% view the grating plot with a small N value to follow the construction
% logic.)
% The x2, x3 coordinate origin is at the center of a pillar.
for n = 1:N
    if n<N
        % The next stripe intercepts a row of pillars between the x3
        % coordinate limits -x(n)*r and -x(n+1)*r.
        stripe.c1 = -x(n+1)*r/grating.d31;
    else
        % The N-th stripe is centered on the pillar axes, and its x3
        % limits are -x(N)*r and +x(N)*r.
        stripe.c1 = x(N)*r/grating.d31;
    end
    % The first block defines the open space between adjacent pillars.
    % The block's x2 coordinate limits are x(N+1-n)*r-d and
    % -x(N+1-n)*r.
    block.c2 = (-x(N+1-n)*r-stripe.c1*grating.d21)/d;
    block.pmt_index = 1;
    stripe.block{1} = block;
    % The second block traverses the pillar interior. Its x2 coordinate
    % limits are -x(N+1-n)*r and +x(N+1-n)*r.
    block.c2 = (x(N+1-n)*r-stripe.c1*grating.d21)/d;
    block.pmt_index = 2;
    stripe.block{2} = block;
    stratum.stripe{n} = stripe;
end
for n = 2:N
    % The next stripe intercepts a row of pillars between the x3
    % coordinate limits x(N+2-n)*r and x(N+1-n)*r.
    stripe.c1 = x(N+1-n)*r/grating.d31;
    % The first block defines the open space between adjacent pillars.
    % The block's x2 coordinate limits are x(n)*r-d and -x(n)*r.
    block.c2 = (-x(n)*r-stripe.c1*grating.d21)/d;
    block.pmt_index = 1;
    stripe.block{1} = block;
    % The second block traverses the pillar interior. Its x2 coordinate
    % limits are -x(n)*r and +x(n)*r.
    block.c2 = (x(n)*r-stripe.c1*grating.d21)/d;
    block.pmt_index = 2;
    stripe.block{2} = block;
    stratum.stripe{N+n-1} = stripe;
end
clear stripe
% The next stripe defines the open space between adjacent rows of
% pillars. Its x3 coordinate limits are x(1)*r and grating.d31-x(1)*r.
stripe.type = 0; % homogeneous stripe
stripe.c1 = 1-x(1)*r/grating.d31;
stripe.pmt_index = 1;
stratum.stripe{2*N} = stripe;
clear stripe
grating.stratum = {stratum};
clear stratum

% Define the indicent field.
clear inc_field
inc_field.wavelength = wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1 = -cosd(theta)/wavelength.
inc_field.f2 = sind(theta)*cosd(phi)/wavelength;
inc_field.f3 = sind(theta)*sind(phi)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
% The orders are specified so that the retained diffracted orders' spatial
% frequencies form a hexagonal pattern such as that illustrated in
% GD-Calc.pdf, Eq's 4.14 and 4.15 and Figure 5.
order = [];
for m2 = -m_max:m_max
    order(end+1).m2 = m2; %#ok<SAGROW> 
    m1 = -m_max:m_max;
    % The following line implicitly applies rectangular index truncation
    % with a different unit cell.
    m1 = m1(abs(m1-m2)<=m_max);
    order(end).m1 = m1;
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
disp(num2str(1-sum([[[R.eff2].']; [[T.eff2].']])));

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75];
pmt_display(2).alpha = 1;
x_limit = [-0.5*h,-1.5*d,-1.3*d;1.5*h,1.5*d,1.3*d];
gdc_plot(grating,1,pmt_display,x_limit);
