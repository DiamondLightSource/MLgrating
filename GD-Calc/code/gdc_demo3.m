% gdc_demo3: Biperiodic checkerboard grating
%
% Test case from https://doi.org/10.1364/JOSAA.14.002758, Example 1 with
% unit cell B.
%
% Note: The number of stripes used in the JOSA paper cannot be plotted; see
% L2 assignment.
%
% gc_demo3 covers the following topics:
%   - GDC interface and setup for a simple biperiodic grating
%   - choice of unit cell and stripe orientations
%   - diffraction order selection
%   - non-trivial polarization geometry
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
disp('gdc_demo3.m')

alt_order = false; % Toggle for alternate order truncation.

% Define parameters for grating structure and incident field:
grating_pmt = 2.25; % grating permittivity
wavelength = 1;
h = wavelength; % grating height
w = 1.25*wavelength; % width of squares (half-period along x2, x3 axes)
m_max = 10; % maximum diffraction order index
if 0
    % This L2 setting matches the JOSA reference, but cannot be plotted.
    L2 = 1000; % number of stripes (must be even)
else
    disp('The number of stripes (L2) is reduced to allow plotting.')
    L2 = 26;
end

% Construct grating.
clear grating
% Define grating material permittivities. The third entry
% (sqrt(grating_pmt)) is used for structural blocks straddling the
% material boundaries. This value represents the geometric-mean
% permittivity within such blocks.
grating.pmt = {1,grating_pmt,sqrt(grating_pmt)};
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
% Define the x2 and x3 projections of the first grating period
% (d21,d31) and second grating period (d22,d32). The x2 and x3 axes are
% aligned to unit cell A, and the period vectors define the edges of
% unit cell B.
grating.d21 = w;
grating.d31 = w;
grating.d22 = -w;
grating.d32 = w;
% Construct the stratum. (View the grating plot with a small L2 value
% to follow the construction logic.)
clear stratum
stratum.type = 2;
stratum.thick = h;
stratum.h11 = 1;
stratum.h12 = 0;
stratum.h21 = 0;
stratum.h22 = 1;
stratum.stripe = cell(1,L2);
for l2 = 1:L2
    clear stripe block
    if l2<=L2/2+1
        stripe.type = 1;
        stripe.block = {};
        block.c2 = (1.5-l2)/L2;
        block.pmt_index = 3;
        stripe.block{end+1} = block;
        if l2>1
            block.c2 = -block.c2;
            block.pmt_index = 2;
            stripe.block{end+1} = block;
            if l2<=L2/2
                block.c2 = block.c2+1/L2;
                block.pmt_index = 3;
                stripe.block{end+1} = block;
            end
        end
        if l2<=L2/2
            block.c2 = 1-block.c2;
            block.pmt_index = 1;
            stripe.block{end+1} = block;
        end
    else
        stripe = stratum.stripe{L2-l2+2};
    end
    stripe.c1 = -0.5+(l2-0.5)/L2;
    stratum.stripe{l2} = stripe;
end
clear stripe block
grating.stratum = {stratum};
clear stratum

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength = wavelength;
inc_field.f2 = 0;
inc_field.f3 = 0;

% Specify which diffracted orders are to be retained in the calculations.
order = [];
for m2 = -m_max:m_max
    order(end+1).m2 = m2; %#ok<SAGROW> 
    m1 = -m_max:m_max;
    if alt_order
        % Alternate order truncation: |m1|+|m2|<=m_max
        m1(abs(m1)+abs(m2)>m_max) = [];
    end
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
% Extract the diffraction order indices for the reflected and transmitted
% waves.
R_m1 = [R.m1].';
R_m2 = [R.m2].';
T_m1 = [T.m1].';
T_m2 = [T.m2].';
% Extract the reflection efficiencies (R1...R4) and transmission
% efficiencies (T1...T4) corresponding to the four incident polarization
% states defined in gdc_eff.m.
R1 = [R.eff1].';
R2 = [R.eff2].';
R3 = [R.eff3].';
R4 = [R.eff4].';
T1 = [T.eff1].';
T2 = [T.eff2].';
T3 = [T.eff3].';
T4 = [T.eff4].';

% Compute diffraction efficiencies (R_e2, T_e2) for incident E field
% polarized parallel to x2 axis (unit basis vector e2).
% Define A, B as described in gdc_eff.m such that e2 = A*s+B*p:
%   A*[s2,s3]+B*[p2,p3] = [1,0], [A,B] = [1,0]/[s2,s3;p2,p3]
%   [p2,p3] = [s3,-s2] (GD-Calc.pdf, Eq 4.20 with fInc in -e1 direction)
%   [A,B] = [1,0]/[s2,s3;s3,-s2] = [s2,s3]
A = inc_field.s2;
B = inc_field.s3;
% GD-Calc.pdf, Eqs A.28, A.30:
R_e2 = R1*abs(A)^2 ...
    +R2*abs(B)^2 ...
    +(2*R3-R1-R2)*real(A*conj(B)) ...
    +(2*R4-R1-R2)*imag(A*conj(B));
T_e2 = T1*abs(A)^2 ...
    +T2*abs(B)^2 ...
    +(2*T3-T1-T2)*real(A*conj(B)) ...
    +(2*T4-T1-T2)*imag(A*conj(B));
% Similarly compute efficiencies (R_e3, T_e3) with incident E field
% parallel to e3:
%   [A,B] = [0,1]/[s2,s3;s3,-s2] = [s3,-s2]
A = inc_field.s3;
B = -inc_field.s2;
R_e3 = R1*abs(A)^2 ...
    +R2*abs(B)^2 ...
    +(2*R3-R1-R2)*real(A*conj(B)) ...
    +(2*R4-R1-R2)*imag(A*conj(B));
T_e3 = T1*abs(A)^2 ...
    +T2*abs(B)^2 ...
    +(2*T3-T1-T2)*real(A*conj(B)) ...
    +(2*T4-T1-T2)*imag(A*conj(B));

% Tabulate the diffraction order indices and diffraction efficiencies. Also
% tabulate the fractional energy loss in the grating.
disp(' ');
disp('Diffraction efficiencies (with incident E parallel to e2 or e3)');
disp('R:');
disp(num2str([R_m1 R_m2 R_e2 R_e3]));
disp('T:');
disp(num2str([T_m1 T_m2 T_e2 T_e3]));
disp('Energy loss:');
disp(num2str(1-sum([[R_e2 R_e3]; [T_e2 T_e3]])));

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75];
pmt_display(2).alpha = 1;
pmt_display(3).name = '';
pmt_display(3).color = [.875,.875,.875];
pmt_display(3).alpha = 1;
x_limit = [-0.5*h,-2.5*w,-2.5*w;1.5*h,2.5*w,2.5*w];
h_plot = gdc_plot(grating,1,pmt_display,x_limit);
view(-30,45)
