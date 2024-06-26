% gdc_demo11: crossed-line grating
%
% Tungsten photonic crystal, based on S. Y. Lin et al, "Highly efficient
% light emission ...", Optics Letters Vol. 28, No. 18, pp. 1683-1685
% (2003).
%
% This script uses the refractive index file W.nk, which must be located
% in a directory specified by the data_path variable. (If data_path is
% empty the refractive index is set to a predetermined constant.)
%
% gdc_demo11 covers the following topics:
%	- biperiodic metallic grating
%   - coordinate break
%   - replication module
%   - harmonic indices
%   - refractive index tables
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

data_path = []; %'../data/'; % can be []

disp(' ')
disp('gdc_demo11.m')

alt_order = true; % Toggle for alternate order truncation.

% Define parameters for grating structure and incident field:
d = 1.5; % grating period
width = 0.5; % grating line width
thick = 0.5; % stratum thickness
rep_count = 4; % replication count (i.e., number of crossed bilayers)
wavelength = 1.825;
m_max = 10; % maximum diffraction order index

if isempty(data_path)
    n = 1.5178+6.4605i;
else
    % Get tabulation wavelengths and refractive indices for metal;
    % interpolate to obtain film permittivity.
    [w,n] = read_nk('W.nk',0.0001,data_path); % tungsten
    n = interp1(w,n,wavelength);
    clear w
end
grating_pmt = n.^2; % grating permittivity
clear n

% Construct grating.
clear grating
grating.pmt = {1,grating_pmt}; % material permittivities
grating.pmt_sub_index = 1; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
grating.d21 = d; % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = d; % second grating period: x3 projection
% Construct the stratum for the grating layer. First construct a cell
% array of strata comprising the replication module.
clear stratum
stratum{1}.type = 1; % module's first stratum is uniperiodic
stratum{1}.thick = thick; % stratum thickness
% The following h11, h12 spec indicates that the stratum's period
% vector matches the first grating period (GD-Calc.pdf, Eqs 3.24
% and 3.25).
stratum{1}.h11 = 1;
stratum{1}.h12 = 0;
clear stripe
stripe.c1 = -0.5*width/d; % first stripe's boundary on positive side
stripe.pmt_index = 1; % first stripe's permittivity index
stratum{1}.stripe{1} = stripe;
stripe.c1 = 0.5*width/d; % second stripe's boundary on positive side
stripe.pmt_index = 2; % second stripe's permittivity index
stratum{1}.stripe{2} = stripe;
clear stripe
stratum{2} = stratum{1}; % Stripes are identical except for orientation.
% The following h11, h12 spec indicates that the stratum's period vector
% matches the second grating period (GD-Calc.pdf, Eq's 3.24 and 3.25).
stratum{2}.h11 = 0;
stratum{2}.h12 = 1;
stratum{3}.type = 3; % module's third stratum is a coordinate break.
stratum{3}.dx2 = d/2; % x2 shift of upper strata
stratum{3}.dx3 = d/2; % x3 shift of upper strata
% Construct the grating's single stratum as a replication module.
grating.stratum{1}.type = 4; % replication module
grating.stratum{1}.stratum = stratum;
grating.stratum{1}.rep_count = rep_count; % replication count
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
% grating. These include evanescent waves.
R = R(imag([scat_field.f1r])==0);
T = T(imag([scat_field.f1t])==0);
% Tabulate the diffraction order indices and diffraction efficiencies for
% the four incident polarization states defined in gdc_eff.m.
disp(' ');
disp('Diffraction efficiencies (m1, m2, eff1, eff2, eff3, eff4)');
disp('R:');
disp(num2str([[R.m1].' [R.m2].' ...
    [R.eff1].' [R.eff2].' [R.eff3].' [R.eff4].']));
disp('T:');
disp(num2str([[T.m1].' [T.m2].' ...
    [T.eff1].' [T.eff2].' [T.eff3].' [T.eff4].']));

% Plot the grating.
clear pmt_display
pmt_display(1).name = ''; % Vacuum
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = ''; % Tungsten
pmt_display(2).color = [1,1,1]*0.75;
pmt_display(2).alpha = 1;
x_limit = [-thick,-(rep_count+1)*thick,-(rep_count+1)*thick;...
    (2*rep_count+1)*thick,(rep_count+1)*thick,(rep_count+1)*thick];
gdc_plot(grating,1,pmt_display,x_limit);
