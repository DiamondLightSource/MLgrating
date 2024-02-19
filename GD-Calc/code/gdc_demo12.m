% gdc_demo12: blazed, phase-Fresnel transmission grating
%
% gdc_demo12 covers the following topic:
%	- uniperiodic grating (blazed surface-relief grating, achromatization)
%
% The grating structure in this demo is a blazed surface relief grating on
% a PMMA prism's exit face, which operates to achromatize the prism.
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

disp(' ')
disp('gdc_demo12.m')

% Prism design (initially without grating):

% PMMA permittivity (refractive index squared) at wavelenth(s) w (micron):
% From https://refractiveindex.info/
%   Shelf: ORGANIC, Book: Polymers - PMMA, Page: Sultanova et al. 2009
% Valid range: 0.4368 – 1.052 micron
pmtPMMA = @(w) 1+1.1819*w.^2./(w.^2-0.011313);
dpmtPMMA = @(w) 2*1.1819*(-0.011313)*w./(w.^2-0.011313).^2; % derivative
wB = 0.55; % blaze wavelength, micron
deltaB = 30; % ray deviation angle at wB, deg
pmtB = pmtPMMA(wB); % PMMA permittivity at wB
nB = sqrt(pmtB); % PMMA refractive index at wB
% prism apex angle (deg), GD-Calc_Demo.pdf, Eq B.7:
a = 2*atand(sind(deltaB/2)./(nB-cosd(deltaB/2)));
% The ray angles in GD-Calc_Demo.pdf, Figure B1, are theta0, theta0_,
% theta1, theta1_. A "B" suffix applies to angles at wavelength wB.
theta0B_ = a/2; % from GD-Calc_Demo.pdf, Eq's B.6 and B.4
theta0B = asind(nB.*sind(theta0B_)); % GD-Calc_Demo.pdf, Eq B.1
theta1B = theta0B_; % GD-Calc_Demo.pdf, Eq B.6
theta1B_ = theta0B; % GD-Calc_Demo.pdf, Eq's B.2, B.6, B.1

% Plot prism dispersion without grating compensator.
w = 0.45:0.01:0.65; % wavelengths
pmt = pmtPMMA(w);
n = sqrt(pmt);
theta0 = theta0B;
theta0_ = asind(sind(theta0)./n); % GD-Calc_Demo.pdf, Eq B.1
theta1 = a-theta0_; % GD-Calc_Demo.pdf, Eq B.4
theta1_ = asind(sind(theta1).*n); % GD-Calc_Demo.pdf, Eq B.2
delta = theta0-theta0_+theta1_-theta1; % GD-Calc_Demo.pdf, Eq B.5
figure, plot(w,delta)

% Prism design with grating: Replace second surface with blazed grating
% with period pd and blaze facet angle b; prism facet angle reduced from a
% to a-b. Design the grating for narrow-band achromatism at wavelength wB.

% Get derivatives, evaluated at w = wB
dn = dpmtPMMA(wB)./(2.*nB); % n = sqrt(pmt), dn/dw = (dpmt/dw)/(2*n)
dtheta0_ = -dn.*sind(theta0B)./nB.^2./cosd(theta0B_); ...
    % GD-Calc_Demo.pdf, Eq B.8
b = atand( ...
    ((wB.*dn-nB).*sind(a-theta0B_) ...
    -nB.*wB.*dtheta0_.*cosd(a-theta0B_) ...
    +sind(deltaB+a-theta0B)) ...
    ./((wB.*dn-nB).*cosd(a-theta0B_) ...
    +nB.*wB.*dtheta0_.*sind(a-theta0B_) ...
    +cosd(deltaB+a-theta0B))); % GD-Calc_Demo.pdf, Eq B.14
pd = wB./(sind(deltaB+a-theta0-b)-nB.*sind(a-theta0B_-b)); ...
    % GD-Calc_Demo.pdf, Eq B.12

disp('Prism, grating design parameters:')
disp(['  a = ' num2str(a) ' deg, b = ' num2str(b) ' deg, wB = ' ...
    num2str(wB) ' micron, pd = ' num2str(pd) ' micron'])

% Plot prism dispersion with grating compensator.
theta1_ = asind(sind(theta1-b).*n+w./pd)+b; % GD-Calc_Demo.pdf, Eq B.4
delta = theta0-theta0_+theta1_-theta1; % GD-Calc_Demo.pdf, Eq B.5
hold on, plot(w,delta)
xlabel('wavelength (micron)')
ylabel('deviation angle (deg)')
legend({'uncompensated prism','with grating compensator'})
title('Prism dispersion without and with grating compensator')

% Grating efficiency analysis:
for with_grating = [false,true]
    if with_grating
        % GD-Calc_Demo.pdf, Figure B1(b)
        N = 50; % number of grating strata
        m_max = 50; % Fourier truncation order for diffraction calculation
        aoi = theta1-b; % internal incidence angle at grating
        b_ = 90; % grating facet wall angle, deg
        % Note: The diffraction efficiency can be improved by increasing b_
        % to 116 deg.
    else
        % Use this branch to analyze the prism's exit surface transmission
        % without the grating, GD-Calc_Demo.pdf, Figure B1(a).
        N = 0;
        m_max = 0;
        aoi = theta1;
    end

    % Create the grating struct.
    clear grating
    grating.pmt = {1,pmt}; % complex permittivities
    grating.pmt_sub_index = 1; % index into pmt for substrate (Vacuum)
    grating.pmt_sup_index = 2; % index into pmt for superstrate (PMMA)
    % period vectors, d_1 = [d21,d31] and d_2 = [d22,d32]:
    grating.d21 = pd;
    grating.d31 = 0;
    grating.d22 = 0;
    grating.d32 = pd;
    grating.stratum = {};
    if N>0
        % Construct a right-triangle sawtooth profile (vacuum space below
        % PMMA superstrate) with triangle base angles b_ and b, base
        % coordinates [x1,x2] = [0,0] and [0,pd] and apex coordinates
        % [x1,x2] = apex:
        %   apex = l1*[sind(b_),cosd(b_)] = [0,pd]+l2*[sind(b),-cosd(b)]
        % (l1 = wall profile length, l2 = facet profile length)
        l = [0,pd]/[sind(b_),cosd(b_);-sind(b),cosd(b)]; % [l1,l2]
        apex = l(1)*[sind(b_),cosd(b_)];
        clear stratum
        stratum.type = 1; % uniperiodic
        stratum.thick = apex(1)/N;
        % h11, h12 are defined so that stratum period matches grating
        % period:
        stratum.h11 = 1;
        stratum.h12 = 0;
        stratum.stripe{1}.pmt_index = 2;
        stratum.stripe{2}.pmt_index = 1;
        for j = 1:N
            t = (j-0.5)/N;
            stratum.stripe{1}.c1 = apex(2)/pd*t;
            stratum.stripe{2}.c1 = 1-(1-apex(2)/pd)*t;
            grating.stratum{end+1} = stratum;
        end
    end

    % Define the incident field.
    clear inc_field
    inc_field.wavelength = w;
    inc_field.f2 = sind(aoi)./(w./n); ...
        % "./n" because aoi is inside the index-n medium
    inc_field.f3 = 0;

    % Define Fourier order truncation.
    clear order
    order.m2 = 0;
    order.m1 = -m_max:m_max;

    % Run diffraction simulation.
    tic
    [~,scat_field,inc_field] = gdc(grating,inc_field,order,false);
    toc

    % Extract order-1 diffraction efficiency (for unpolarized
    % illumination).
    [R,T] = gdc_eff(scat_field,inc_field);
    if m_max>0
        j = find([scat_field.m1]==1);
    else
        j = find([scat_field.m1]==0);
    end
    T_ = 0.5*(T(j).eff1+T(j).eff2); % average for unpolarized

    if ~with_grating
        figure, plot(w,T_)
        continue
    end

    % Plot the efficiency (T) versus wavelength (w).
    hold on, plot(w,T_)
    xlabel('wavelength (micron)')
    ylabel('efficiency')
    title('Transmission efficiency without and with grating compensator')
    legend('uncompensated prism','with grating compensator', ...
        'Location','southeast')
    Position = get(gcf,'Position');
    set(gcf,'Position',[Position(1:3),300])

    % Tabulate [w,T] pairs.
    disp('Wavelength, Transmission in order 1:')
    disp(num2str([w.',T_.']))
    % Check energy balance.
    eff1 = [cat(1,R.eff1);cat(1,T.eff1)];
    eff2 = [cat(1,R.eff2);cat(1,T.eff2)];
    disp(['Energy imbalance: ' ...
        num2str(max(max(abs(1-sum(eff1,1))),max(abs(1-sum(eff2,1)))))]);

    % Plot the grating structure.
    x_limit = [-.25;1.25]*[1,1,1]*pd;
    gdc_plot(grating,[],[],x_limit);
    view(0,0)
    axis image
    set(gca,'ZLim',[0,apex(1)])
    set(gca,'ZTick',[0,apex(1)])
    Position = get(gcf,'Position');
    set(gcf,'Position',[Position(1:3),100])
end
