% gdc_demo10: holographic volume grating
%
% This demo models a holographic transmission volume grating (1st-order
% Bragg diffraction, external surface reflections are neglected).
%
% gdc_demo10 covers the following topics:
%   - parameterization
%   - coordinate break
%   - replication module
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo10.m')

transmission = true; % toggle for transmission or reflection grating
% CAUTION: gdc_plot is about 4X slower with transmission = false.

% Define parameters for grating structure and incident field. (Length
% quantities are in micron units.) The wavelength parameter is vectorized
% in dimension 1.
pmt = 2.25; % grating mean permittivity
% Note: the grating substrate and superstrate also have permittivity pmt;
% this neglects external surface reflections.
dpmt = 0.15; % permittivity modulation amplitude, peak-to-mean
thick = 10.0; % grating layer thickness
theta0 = 0; % incidence angle (deg) inside permittivity-pmt medium
slant = 20; % grating slant angle (deg)
% permittivity discretization steps:
pmt_steps = 32;
% rep_count ("replication count") is the stratification number for the
% grating layer (preferably a power of 2):
rep_count = 64;
Bragg_wavelength = 1.0; % Bragg wavelength (in vacuum)
wavelength = (0.75:0.01:1.25).'; % simulation wavelengths (in vacuum)

if ~transmission
    slant = 90-slant;
    rep_count = 32*rep_count;
    wavelength = (wavelength-Bragg_wavelength)/4+Bragg_wavelength;
end

% Calculate grating period (pd) from Bragg_wavelength.
theta1 = 2*slant-theta0; % geometric diffraction angle
% sind(theta1)-sind(theta0) = Bragg_wavelength/sqrt(pmt)/pd
pd = Bragg_wavelength/sqrt(pmt)/(sind(theta1)-sind(theta0));
disp(['Grating''s surface-tangential period: ' num2str(pd)])

% Construct grating. (Note: Two grating periods must be specified
% because the grating can generally be biperiodic, although for this
% example the second period (d22,d32) is irrelevant.)
clear grating
grating.pmt = num2cell( ...
    [pmt,pmt+cos(pi/pmt_steps*(-pmt_steps+1:2:pmt_steps-1))*dpmt]); ...
    % material permittivities
grating.pmt_sub_index = 1; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index
grating.d21 = pd; % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = pd; % second grating period: x3 projection
% Construct the stratum for the grating layer. First construct a cell
% array of strata for one unit of the replication module.
clear stratum
stratum{1}.type = 1; % module's first stratum is uniperiodic
stratum{1}.thick = thick/rep_count; % module thickness
% The following h11, h12 spec indicates that the stratum's period vector
% matches the first grating period (GD-Calc.pdf, Eq's 3.24 and 3.25).
stratum{1}.h11 = 1;
stratum{1}.h12 = 0;
stripe = repmat({struct('c1',[],'pmt_index',[])},1,pmt_steps);
for j = 1:pmt_steps
    stripe{j}.c1 = j/pmt_steps; % stripe boundary on positive side
    stripe{j}.pmt_index = 1+j; % index into grating.pmt
end
stratum{1}.stripe = stripe;
stratum{2}.type = 3; % module's second stratum is a coordinate break.
stratum{2}.dx2 = -stratum{1}.thick*tand(slant); % x2 shift of upper strata
stratum{2}.dx3 = 0; % x3 shift of upper strata
% Construct the grating layer as a replication module.
grating.stratum{1}.type = 4; % replication module
grating.stratum{1}.stratum = stratum;
grating.stratum{1}.rep_count = rep_count; % replication count

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength = wavelength;
% Include refractive index (sqrt(pmt)) factor in f2 to convert the vacuum
% wavelength to the internal wavelength. (theta0 is the internal incidence
% angle.)
inc_field.f2 = sind(theta0)./(wavelength./sqrt(pmt));
inc_field.f3 = 0;

% Specify which diffracted orders are to be retained in the calculations.
% (m2 is zero because this is a uniperiodic grating - all diffraction
% orders for m2~=0 have zero amplitude.)
clear order
order(1).m2 = 0;
order(1).m1 = 0:1; % zero and first orders only

% Run the diffraction calculations.
tic
[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
toc

% Calculate and plot the diffraction efficiencies via Kogelnik.
[effTE,effTM,transmission] = ...
    Kogelnik(thick,-pd.*cosd(slant),pmt,dpmt,wavelength,theta0,slant+90);
figure, plot(wavelength, effTE,'b--');
hold on, plot(wavelength, effTM,'r--');

% Calculate and plot the diffraction efficiencies via gdc.
[R,T] = gdc_eff(scat_field,inc_field);
% Extract the diffraction order indices for the reflected and transmitted
% waves. (Only the m1 indices matter; the m2's are all zero.)
R_m1 = [R.m1];
T_m1 = [T.m1];
% Extract the reflection efficiencies (R1...R4) and transmission
% efficiencies (T1...T4) corresponding to the four incident polarization
% states defined in gdc_eff.m.
R1 = [R.eff1];
R2 = [R.eff2];
T1 = [T.eff1];
T2 = [T.eff2];
% Display and plot first order.
i = find(T_m1==1);
disp(' ');
p1 = 1:5:length(wavelength); % Select wavelength tabulation indices.
if transmission
    disp(['Wavelengths, diffraction efficiencies ' ...
        '(TE & TM, 1st order, transmission)']);
    disp(num2str([wavelength(p1),T1(p1,i),T2(p1,i)]));
    plot(wavelength,[T1(:,i)],'b');
    plot(wavelength,[T2(:,i)],'r');
    title('Efficiency in 1st order (transmission)');
else
    disp(['Wavelengths, diffraction efficiencies ' ...
        '(TE & TM, 1st order, reflection)']);
    disp(num2str([wavelength(p1),R1(p1,i),R2(p1,i)]));
    plot(wavelength,[R1(:,i)],'b');
    plot(wavelength,[R2(:,i)],'r');
    title('Efficiency in 1st order (reflection)');
end
xlabel('wavelength (micron)');
ylabel('efficiency');
legend('TE (Kogelnik)','TM (Kogelnik)','TE (gdc)','TM (gdc)')
Position = get(gcf,'Position');
set(gcf,'Position',[Position(1:2),560,350])
drawnow

disp('Energy loss (TE & TM):');
disp(num2str([max(abs(1-sum([R1 T1],2))) max(abs(1-sum([R2 T2],2)))]))

% Modify grating parameters for plotting.
pmt_steps = 8;
rep_count = 64;
if ~transmission
    rep_count = 4*rep_count;
end
grating.pmt = num2cell( ...
    [pmt,pmt+cos(pi/pmt_steps*(-pmt_steps+1:2:pmt_steps-1))*dpmt]); ...
stripe = repmat({struct('c1',[],'pmt_index',[])},1,pmt_steps);
for j = 1:pmt_steps
    stripe{j}.c1 = j/pmt_steps; % stripe boundary on positive side
    stripe{j}.pmt_index = 1+j; % index into grating.pmt
end
stratum{1}.stripe = stripe;
stratum{1}.thick = thick/rep_count;
stratum{2}.dx2 = -stratum{1}.thick*tand(slant);
grating.stratum{1}.stratum = stratum;
grating.stratum{1}.rep_count = rep_count;

% Plot the grating.
clear pmt_display
for j = length(grating.pmt):-1:1
    pmt_display(j).name = '';
    pmt_display(j).color = (1+(grating.pmt{j}-pmt)/dpmt)/2*[1,1,1];
    pmt_display(j).alpha = 0.5;
end
[~,param_index] = min(abs(wavelength-Bragg_wavelength));
gdc_plot(grating,param_index,pmt_display,[],[],inf);
