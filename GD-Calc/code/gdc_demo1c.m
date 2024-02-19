% gdc_demo1c: Uniperiodic, sinusoidal grating
%
% Test case from https://doi.org/10.1364/JOSAA.10.002581, Table 4 with
% h/d  =  100. (This script is adapted from gdc_demo1a.) Toggle the
% "branch = ..." statement for the TE or TM case. The TE data in Table 4
% corresponds to eff1, and TM corresponds to eff2.
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo1c.m')

branch = true; % true for TE case, false for TM case

% If ctr_sect==true, define stratum widths by center-sectioning; otherwise
% use (slightly) more accurate averaging.
ctr_sect = true;

% Define parameters for grating structure and incident field:
grating_pmt = (0.3+7.0i)^2; % grating permittivity
d = 1.0; % grating period
wavelength = d/1.7;
h = 100*d; % grating height
if branch
    disp('TE case in JOSA publication')
    L1 = 50; % number of grating strata (for "staircase" approximation)
    m_max = 25; % max diffraction order index
else
    disp('TM case in JOSA publication')
    L1 = 10; % number of grating strata (for "staircase" approximation)
    m_max = 52; % max diffraction order index
end
% theta1 = angle between the grating normal (x1 axis) and the incident wave
% vector
theta1 = 30; % deg
% phi1 = angle between the grating-tangential direction normal to the
% grating lines (x2 axis) and the incident wave vector's grating-tangential
% projection (in the x2, x3 plane)
phi1 = 0; % deg

if ctr_sect
    % Define stratum half-widths (in units of d) by center-sectioning.
    c1 = acos(-1+2*((1:L1)-0.5)/L1)/(2*pi);
else
    % Define stratum half-widths using slightly more accurate averaging
    % method.
    c1 = 2*(1:L1-1)/L1-1;
    c1 = diff([-0.25, (c1.*acos(c1)-sqrt(1-c1.^2))/(4*pi), 0])*L1;
end

% Construct grating. (Note: Two grating periods must be specified
% because the grating can generally be biperiodic, although for this
% example the second period (d22,d32) is irrelevant.)
clear grating
grating.pmt = {1,grating_pmt}; % grating material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
grating.d21 = d; % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = d; % second grating period: x3 projection
grating.stratum = {};
clear stratum stripe
stratum.type = 1; % uniperiodic stratum
stratum.thick = h/L1; % stratum thickness
% The following h11, h12 spec indicates that the stratum's period
% vector matches the first grating period (GD-Calc.pdf, Eq's 3.23-25).
stratum.h11 = 1;
stratum.h12 = 0;
for l1 = 1:L1
    stripe.pmt_index = 1; % first stripe's permittivity index
    stripe.c1 = -c1(l1); % first stripe's boundary on positive side
    stratum.stripe{1} = stripe;
    stripe.pmt_index = 2; % second stripe's permittivity index
    stripe.c1 = c1(l1); % second stripe's boundary on positive side
    stratum.stripe{2} = stripe;
    grating.stratum{end+1} = stratum;
end
clear c1 stratum stripe

% Define the indicent field.
clear inc_field
inc_field.wavelength = wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector, Eq's 4.1-2 in
% GD-Calc.pdf. The grating-normal (x1) projection is implicitly
%   f1 = -cosd(theta1)/wavelength
inc_field.f2 = sind(theta1)*cosd(phi1)/wavelength;
inc_field.f3 = sind(theta1)*sind(phi1)/wavelength;

% Specify which diffracted orders are to be retained in the calculations,
% Eq's 4.12-13 in GD-Calc.pdf. (m2 is zero because this is a uniperiodic
% grating - all diffraction orders for m2~=0 have zero amplitude.)
clear order
order(1).m2 = 0;
order(1).m1 = -m_max:m_max;

% Run the diffraction calculations.
tic
[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
toc

% Compute the diffraction efficiencies. (Only the reflected waves are
% retained.)
R = gdc_eff(scat_field,inc_field);
% Discard diffracted waves that decay exponentially with distance from the
% grating. These include evanescent waves.
% (Note: "[scat_field.f1r]" is MATLAB syntax for
% "[scat_field(1).f1r, scat_field(2).f1r, ...]".)
R = R(imag([scat_field.f1r])==0);
% Extract the diffraction order indices for the reflected waves. (Only the
% m1 indices matter; the m2's are all zero.)
R_m1 = [R.m1].';
% Extract the reflection efficiencies (R1...R4) corresponding to the four
% incident polarization states defined in gdc_eff.m.
R1 = [R.eff1].';
R2 = [R.eff2].';
R3 = [R.eff3].';
R4 = [R.eff4].';
% Tabulate the diffraction order indices and diffraction efficiencies.
disp(' ');
disp('Diffraction efficiencies (m1, eff1, eff2, eff3, eff4)');
disp('R:');
disp(num2str([R_m1 R1 R2 R3 R4]));
