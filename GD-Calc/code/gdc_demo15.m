% gdc_demo15: Uniperiodic, sinusoidal grating with internal field
% calculation
%
% Test case from section VI.5 of Light Propagation in Periodic Media, by
% Michel Nevière, Evgeny Popov (MarcelDekker, Inc. New York, 2003).
% https://doi.org/10.1201/9781482275919
%
% Sinusoidal aluminum grating, approximaed by 10 lamellar strata: Plot
% cross-section of |H3|, |E1|, and |E2| through grating interior for TM
% polarization.
%
% This demo illsutrates the problem of field singularities at grating edges
% and convergence limitations of the staircase approximation. The
% convergence problem can also exist, in some cases, with lamellar
% gratings; see, for example:
%   Li, L., & Granet, G. (2011). Field singularities at lossless
%   metal-dielectric right-angle edges and their ramifications to the
%   numerical modeling of gratings. JOSA A, 28(5), 738-746.
%   https://doi.org/10.1364/JOSAA.28.000738
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo15.m')

% Define parameters for grating structure and incident field:
grating_pmt = (1.3+7.6i)^2; % grating permittivity (Al)
d = 0.5; % grating period
wavelength = 0.6328;
h = 0.2; % grating height
steps = 10; % number of steps for staircase approximation
% step_partitions = number of stratum partitions per step (for full-field
% sampling)
step_partitions = 20;
m_max = 80; % max diffraction order index
% theta1 = angle between the grating normal (x1 axis) and the incident wave
% vector
theta1 = 40; % deg
% phi1 = angle between the grating-tangential direction normal to the
% grating lines (x2 axis) and the incident wave vector's grating-tangential
% projection (x2, x3 plane)
phi1 = 0; % deg

% Define stratum half-widths (in units of d) by center-sectioning.
c1 = acos(-1+2*((1:steps)-0.5)/steps)/(2*pi);

% Construct grating. (Note: Two grating periods must be specified because
% the grating can generally be biperiodic, although for this example the
% second period (d22,d32) is irrelevant.)
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
stratum.thick = h/steps/step_partitions; % stratum thickness
% The following h11, h12 spec indicates that the stratum's period vector
% matches the first grating period (GD-Calc.pdf, equations 3.22 and 3.23).
stratum.h11 = 1;
stratum.h12 = 0;
for l1 = 1:steps
    stripe.pmt_index = 1; % first stripe's permittivity index
    stripe.c1 = -c1(l1); % first stripe's boundary on positive side
    stratum.stripe{1} = stripe;
    stripe.pmt_index = 2; % second stripe's permittivity index
    stripe.c1 = c1(l1); % second stripe's boundary on positive side
    stratum.stripe{2} = stripe;
    stratum.full_field = [];
    % Put in a zero-thickness stratum to get the boundary field at the
    % bottom of the step.
    grating.stratum{end+1} = stratum;
    grating.stratum{end}.thick = 0;
    % Split the step into step_partitions strata.
    for j = 1:step_partitions
        grating.stratum{end+1} = stratum;
    end
end
clear c1 stratum stripe
% Add strata for sampling the E and H fields in the superstrate.
stratum.type = 0;
stratum.pmt_index = 1;
stratum.thick = h/steps/step_partitions;
stratum.full_field = [];
for l1 = 1:100
    grating.stratum{end+1} = stratum;
end
clear stratum
grating.sub_full_field = [];
grating.sup_full_field = [];

% Define the indicent field.
clear inc_field
inc_field.wavelength = wavelength;
% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector. The grating-normal (x1)
% projection is implicitly f1 = -cosd(theta1)/wavelength.
inc_field.f2 = sind(theta1)*cosd(phi1)/wavelength;
inc_field.f3 = sind(theta1)*sind(phi1)/wavelength;

% Specify which diffracted orders are to be retained in the calculations.
% (m2 is zero because this is a uniperiodic grating - all diffraction
% orders for m2~=0 have zero amplitude.)
clear order
order(1).m2 = 0;
order(1).m1 = -m_max:m_max;

% Run the diffraction calculations (with grating output to capture
% full-field data).
tic
[param_size,scat_field,inc_field,grating] = ...
    gdc(grating,inc_field,order,false);
toc

% Compute the diffraction efficiencies. (Only the reflected waves are
% retained.)
R = gdc_eff(scat_field,inc_field);
% Discard diffracted waves that decay exponentially with distance from the
% grating. These include evanescent waves.
R = R(imag([scat_field.f1r])==0);
% Tabulate the diffraction order indices and diffraction efficiencies for
% TE and TM polarization
disp(' ');
disp('Diffraction efficiencies (m1, eff1, eff2)');
disp('R:');
disp(num2str([[R.m1].' [R.eff1].' [R.eff2].']));

% Calculate spatial field maps.
x2 = (0:.002:d).';
x3 = 0;
x1 = [];
H3p = [];
E1p = [];
E2p = [];
disp(' ')
disp('Doing Fourier transforms ...')
tic
x1_ = 0;
for l1 = 1:length(grating.stratum)
    x1_ = x1_+grating.stratum{l1}.thick;
    x1(end+1) = x1_; %#ok<SAGROW>
    map = EH_map(rmfield(grating.stratum{l1}.full_field,{...
        'ffE1s','ffE2s','ffE3s',...
        'ffH1s','ffH2s','ffH3s',...
        'ffE3p','ffH1p','ffH2p'...
        }),1,x2,x3);
    H3p(:,end+1) = [map.H3p].'; %#ok<SAGROW>
    E1p(:,end+1) = [map.E1p].'; %#ok<SAGROW>
    E2p(:,end+1) = [map.E2p].'; %#ok<SAGROW>
end
toc
[x1,x2] = deal(repmat(x1,length(x2),1),repmat(x2,1,length(x1)));

% Plot the grating.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = '';
pmt_display(2).color = [.75,.75,.75;0,0,0];
pmt_display(2).alpha = 1;
x_limit = [-0.5*h,-d,-d;1.5*h,d,d];
gdc_plot(grating,1,pmt_display,x_limit);
view(127.5,20)

% Show 3-D perspective view of the field magnitude. (To see a
% cross-sectional elevation view, use view(2).)
figure
surface(x2,x1,abs(H3p),'EdgeAlpha',.25,'FaceColor','interp'), view(-45,60);
xlabel('x2');
ylabel('x1');
zlabel('|H3|');
title('|H3|');
Position = get(gcf,'Position');
set(gcf,'Position',[Position(1:3),336]);
figure
surface(x2,x1,abs(E1p),'EdgeAlpha',.25,'FaceColor','interp'), view(-45,60);
xlabel('x2');
ylabel('x1');
zlabel('|E1|');
title('|E1|');
Position = get(gcf,'Position');
set(gcf,'Position',[Position(1:3),336]);
figure
surface(x2,x1,abs(E2p),'EdgeAlpha',.25,'FaceColor','interp'), view(-45,60);
xlabel('x2');
ylabel('x1');
zlabel('|E2|');
title('|E2|');
Position = get(gcf,'Position');
set(gcf,'Position',[Position(1:3),336]);
