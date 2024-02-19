% gdc_demo4: Biperiodic checkerboard grating
%
% Test case from https://doi.org/10.1364/JOSAA.14.002758, Example 1 with
% unit cell A.
%
% gc_demo4 covers the following topics:
%   - GDC interface and setup for a simple biperiodic grating
%   - choice of unit cell and stripe orientations
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
disp('gdc_demo4.m')

alt_order = false; % Toggle for alternate order truncation.

% Define parameters for grating structure and incident field:
grating_pmt = 2.25; % grating permittivity
wavelength = 1;
h = wavelength; % grating height
w = 1.25*wavelength; % width of squares (half-period along x2, x3 axes)
m_max = 10; % maximum diffraction order index

% Construct grating.
clear grating
grating.pmt = {1,grating_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
% Define the x2 and x3 projections of the first grating period
% (d21,d31) and second grating period (d22,d32). The x2 and x3 axes are
% aligned to unit cell A, and the period vectors define the unit cell's
% edges.
grating.d21 = 2*w;
grating.d31 = 0;
grating.d22 = 0;
grating.d32 = 2*w;
% Construct the stratum. (View the grating plot to follow the construction
% logic.)
clear stratum stripe block
stratum.type = 2;
stratum.thick = h;
stratum.h11 = 1;
stratum.h12 = 0;
stratum.h21 = 0;
stratum.h22 = 1;
stripe.c1 = 0.25;
stripe.type = 1;
block.c2 = 0.25;
block.pmt_index = 2;
stripe.block{1} = block;
block.c2 = 0.75;
block.pmt_index = 1;
stripe.block{2} = block;
stratum.stripe{1} = stripe;
stripe.c1 = 0.75;
block.c2 = 0.25;
stripe.block{1} = block;
block.c2 = 0.75;
block.pmt_index = 2;
stripe.block{2} = block;
stratum.stripe{2} = stripe;
grating.stratum{1} = stratum;
clear stratum stripe block

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength = wavelength;
inc_field.f2 = 0;
inc_field.f3 = 0;

% Specify which diffracted orders are to be retained in the calculations.
% First define order indices relative to unit cell B.
m1 = repmat(-m_max:m_max,[2*m_max+1,1]);
m2 = m1.';
m1 = m1(:).';
m2 = m2(:).';
if alt_order
    % Alternate order truncation: |m1|+|m2|<=m_max
    find_ = find(abs(m1)+abs(m2)<=m_max);
    m1 = m1(find_);
    m2 = m2(find_);
end
% Redefine order indices relative to unit cell A (GD-Calc_Demo.pdf, Eq's 18
% and 19).
[m1,m2] = deal(m1-m2,m1+m2);
% Construct the order struct.
order = [];
for m2_ = unique(m2)
    order(end+1).m2 = m2_; %#ok<SAGROW> 
    order(end).m1 = m1(m2==m2_);
end
clear m1 m2

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
% Extract the diffraction order indices for the reflected and transmitted
% waves.
R_m1 = [R.m1].';
R_m2 = [R.m2].';
T_m1 = [T.m1].';
T_m2 = [T.m2].';
% Redefine order indices relative to unit cell B (GD-Calc_Demo.pdf, Eq's 22
% and 23).
[R_m1,R_m2] = deal((R_m1+R_m2)/2,(R_m2-R_m1)/2);
[T_m1,T_m2] = deal((T_m1+T_m2)/2,(T_m2-T_m1)/2);
% Sort the orders by order index
[~,k] = sortrows([R_m1,R_m2],[2,1]);
R_m1 = R_m1(k);
R_m2 = R_m2(k);
R = R(k);
[~,k] = sortrows([T_m1,T_m2],[2,1]);
T_m1 = T_m1(k);
T_m2 = T_m2(k);
T = T(k);
% Extract the reflection efficiencies (R_e2,R_e3) and transmission
% efficiencies (T_e2,T_e3), where R_e2 and T_e2 correspond to an incident
% field polarized in the e2 direction, and R_e3 and T_e3 correspond to an
% incident field polarized in the e3 direction. For the current diffraction
% geometry, e2 and e3 are equal to the polarization basis vectors p and s,
% respectively (GD-Calc_Demo.pdf, Eq's 13-14, and GD-Calc.pdf, Eq's 4.2,
% 4.9, 4.19-20), so the e2 and e3 polarization efficiencies correspond
% respectively to eff2 and eff1.
R_e2 = [R.eff2].';
R_e3 = [R.eff1].';
T_e2 = [T.eff2].';
T_e3 = [T.eff1].';

% Tabulate the diffraction order indices and diffraction efficiencies for
% an incident field polarized parallel to the x2 or x3 axis (unit basis
% vector e2 or e3). Also tabulate the fractional energy loss in the
% grating.
disp(' ');
disp('Diffraction efficiencies (with incident E parallel to e2 or e3)');
disp('R:');
disp(num2str([R_m1 R_m2 R_e2 R_e3]));
disp('T:');
disp(num2str([T_m1 T_m2 T_e2 T_e3]));
disp('Energy loss:');
disp(num2str(1-sum([[R_e3 R_e2]; [T_e3 T_e2]])));

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75];
pmt_display(2).alpha = 1;
x_limit = [-0.5*h,-2.5*w,-2.5*w;1.5*h,2.5*w,2.5*w];
h_plot = gdc_plot(grating,1,pmt_display,x_limit);
view(-30,45)
