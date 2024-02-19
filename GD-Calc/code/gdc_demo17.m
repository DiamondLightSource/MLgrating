% gdc_demo17: diffractive beam splitter, Gaussian beam
%
% gdc_demo17 covers the following topic:
%	- multiple incident orders
%   - Gaussian incident beam
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%   https://en.wikipedia.org/wiki/Gaussian_beam
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% Define grating structure parameters.
n = 1.5; % refractive index in superstrate and grating
lambda = 0.633; % wavelength, micron (in vacuum)
t = lambda/2/(n-1); % grating thickness (set to extinguish zero order)
theta = 10.0; % first-order diffraction angle (deg)
period = lambda/sind(theta); % grating period
c = 0.5; % (line width)/period

% Define Gaussian beam parameters.
% The beam's spatial amplitude profile is exp(-(x2^2+x3^2)/w0^2).
% The frequency amplitude profile (Fourier transform) is
% exp(-pi^2*(f2^2+f3^2)*w0^2).
NA0 = 0.1; % beam numerical aperture at 1/e^2 intensity
% NA0 = sine of beam convergence angle in superstrate, times refractive
% index.
% NA0/lambda = sqrt(f2^2+f3^2) at 1/e amplitude; pi*(NA0/lambda)*w0 = 1
w0 = lambda/(pi*NA0); % beam radius at 1/e^2 intensity
% weps = beam radius at eps intensity: exp(-weps^2/w0^2)^2 = eps.
weps = sqrt(-log(eps)/2)*w0;
% feps = frequency radius at eps intensity: exp(-pi^2*feps^2*w0^2)^2 = eps.
feps = sqrt(-log(eps)/2)/pi/w0;
% The beam's spatial profile will be periodically tiled in a "supercell"
% array with period 2*weps, and will be frequency-sampled at increments
% 1/(2*weps). Npd = number of grating periods per supercell period.
Npd = 10*ceil(weps/period);
m_max = 10; % maximum diffraction order index

clear grating
grating.pmt = {1,n.^2};
grating.pmt_sub_index = 1;
grating.pmt_sup_index = 2;
grating.d21 = period;
grating.d31 = 0;
grating.d22 = 0;
grating.d32 = period;
clear stratum
stratum.type = 1;
stratum.thick = t;
stratum.h11 = 1;
stratum.h12 = 0;
clear stripe
stripe.c1 = -c/2;
stripe.pmt_index = 1;
stratum.stripe{1} = stripe;
stripe.c1 = c/2;
stripe.pmt_index = 2;
stratum.stripe{2} = stripe;
grating.stratum{1} = stratum;
clear stratum stripe

clear inc_field
inc_field.wavelength = lambda;
% Divide the grating frequency (1/period) into Npd frequency steps per
% incident order.
df = 1/(Npd*period); % freuency sampling step
% f2 uniformly samples the frequency range -0.5/period to 0.5/period at the
% centers of Npd intervals.
inc_field.f2 = (-Npd/2+0.5:Npd/2)*df; % row vector
% f3 is sampled at the same step df and extends to the eps intensity
% threshold.
N = ceil(feps/df);
inc_field.f3 = (-N:N).'*df; % column vector

clear order
order(1).m2 = 0;
order(1).m1 = -m_max:m_max;

tic
[~,scat_field,inc_field] = gdc(grating,inc_field,order,false,true);
toc

% Exclude incident orders with f2 outside of the eps intensity threshold
% (feps).
f2 = cat(1,inc_field.f2); % f2(k_inc,j) = inc_field(k_inc).f2(j)
k_inc = find(any(abs(f2)<=feps,2));
inc_field = inc_field(k_inc);
scat_field = scat_field(k_inc,:);

% Define wave amplitudes (including integration weighting factors, OF) for
% the incident field.
for k_inc = 1:length(inc_field)
    f1 = inc_field(k_inc).f1;
    f2 = inc_field(k_inc).f2;
    f3 = inc_field(k_inc).f3;
    s2 = inc_field(k_inc).s2;
    p2 = inc_field(k_inc).p2;
    % incident obliquity factor, GD-Calc_Demo.pdf, Eq 35:
    OF = -f1*(lambda/n);
    % GD-Calc_Demo.pdf, Eq 41:
    u = exp(-pi^2*w0^2*(f2.^2+f3.^2))./(OF.*sqrt(p2.^2+s2.^2));
    % GD-Calc_Demo.pdf, Eq 39:
    inc_field(k_inc).Asi = u.*p2;
    inc_field(k_inc).Api = -u.*s2;
end

% Plot the incident intensity map Ii.
f2 = [inc_field.f2]; % merged sampling range
f3 = inc_field(1).f3; % same for all inc_field elements because order.m2==0
Asi = [inc_field.Asi];
Api = [inc_field.Api];
Ii = abs(Asi).^2+abs(Api).^2;
figure, surf(lambda/n*f2,lambda/n*f3,Ii,'EdgeColor','none', ...
    'FaceColor','w','FaceLighting','gouraud')
zticks([]), daspect([1,1,2.5]), light('Position',[0,-1,1])
xlabel('\lambda/n*f_{2}')
ylabel('\lambda/n*f_{3}')
title('Incident field')

% Calculate the diffracted field separately for each inc_field element.
for k_inc = 1:size(scat_field,1)
    Asi = inc_field(k_inc).Asi;
    Api = inc_field(k_inc).Api;
    for k = 1:size(scat_field,2)
        % summand in GD-Calc_Demo.pdf, Eq 45:
        scat_field(k_inc,k).Ast = scat_field(k_inc,k).Tss.*Asi ...
            +scat_field(k_inc,k).Tsp.*Api;
        scat_field(k_inc,k).Apt = scat_field(k_inc,k).Tps.*Asi ...
            +scat_field(k_inc,k).Tpp.*Api;
    end
end
scat_field = rmfield(scat_field, ...
    {'m1_inc','m2_inc','Tss','Tsp','Tps','Tpp','Rss','Rsp','Rps','Rpp'});
% All remaining fields of scat_field(1,k), scat_field(2,k), ... except Ast
% and Apt are identical (no dependence on incident field). Sum Ast and Apt
% over incident fields.
for k = 1:size(scat_field,2)
    for k_inc = 2:size(scat_field,1)
        scat_field(1,k).Ast = scat_field(1,k).Ast+scat_field(k_inc,k).Ast;
        scat_field(1,k).Apt = scat_field(1,k).Apt+scat_field(k_inc,k).Apt;
    end
end
scat_field(2:end,:) = [];

% Plot the transmitted far-field intensity map It (including obliquity
% factor).
f2 = [scat_field.f2]; % merged sampling range
f3 = scat_field(1).f3; % same for all scat_field elements
f1t = [scat_field.f1t];
Ast = [scat_field.Ast];
Apt = [scat_field.Apt];
% transmitted obliquity factor, GD-Calc_Demo.pdf, Eq 46:
OF = -real(f1t)*lambda;
% Convert transmitted field's Fourier amplitudes into far-field amplitudes
% (with OF accounting for projected beam area in the propagation
% direction):
Ast = Ast.*OF;
Apt = Apt.*OF;
% far-field intensity map:
It = abs(Ast).^2+abs(Apt).^2;
% Convert evanescent-wave intensities to nan's so they don't plot:
It(imag(f1t(:))~=0) = nan;
figure, surf(lambda*f2,lambda*f3,It,'EdgeColor','none', ...
    'FaceColor','w','FaceLighting','gouraud')
zticks([]), light('Position',[0,-1,0.4]), daspect([1,1,1])
xlabel('\lambda*f_{2}')
ylabel('\lambda*f_{3}')
title('Diffracted field')

