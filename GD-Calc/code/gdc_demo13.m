% gdc_demo13: alignment sensor with internal field calculation
%
% Generate and display a parameterized electromagnetic field map.
% (Adapted from gdc_demo9.m.)
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% If results_dir is nonempty a video of the lateral grating scan and H3p
% field map will be saved to the file [results_dir '/gdc_demo13.avi'].
results_dir = '../results';
if ~isempty(results_dir) && ~exist(results_dir,'dir')
    warning(['Results directory ''' results_dir ''' does not exist.'])
    results_dir = [];
end

disp(' ')
disp('gdc_demo13.m')

% Define parameters for grating structure and incident field. The phase
% plate's lateral shift parameter (dx2) is vectorized in the first
% dimension.
superstrate_pmt = 2.25; % superstrate (phase plate) permittivity
substrate_pmt = (1.37+7.62i)^2; % substrate permittivity (Aluminum)
wavelength = 0.633; % micron
t1 = wavelength/4; % bottom (reflecting) grating thickness
t2 = wavelength*2.5; % air space thickness
t3 = (wavelength/8)/(sqrt(superstrate_pmt)-1); % phase plate thickness
period = 10*wavelength; % grating period
dx2 = period*(0:63).'/64; ...
    % phase plate's lateral shift (vectorized in dimension 1)
m_max = 20; % maximum diffraction order index

% Construct grating.
clear grating
grating.pmt = ...
    {1,substrate_pmt,superstrate_pmt}; % material permittivities
grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 3; % superstrate permittivity index
% Note: Two grating periods must be specified because the grating can
% generally be biperiodic, although for this example the second period
% (d22,d32) is irrelevant.
grating.d21 = period; % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = period; % second grating period: x3 projection
% Construct the stratum for the reflecting grating.
clear stratum
stratum.type = 1; % uniperiodic stratum
stratum.thick = t1; % stratum thickness
% The following h11, h12 spec indicates that the stratum's period vector
% matches the first grating period (GD-Calc.pdf, Eq's 3.24 and 3.25).
stratum.h11 = 1;
stratum.h12 = 0;
clear stripe
stripe.c1 = -0.25; % first stripe's boundary on positive side
stripe.pmt_index = 1; % first stripe's permittivity index
stratum.stripe{1} = stripe;
stripe.c1 = 0.25; % second stripe's boundary on positive side
stripe.pmt_index = 2; % second stripe's permittivity index
stratum.stripe{2} = stripe;
grating.stratum{1} = stratum;
% Construct the stratum for the air space.
clear stratum
stratum.type = 0; % homogeneous stratum
stratum.pmt_index = 1; % stratum's permittivity index
stratum.thick = t2; % stratum thickness
grating.stratum{2} = stratum;
% Construct a (zero-thickness) stratum representing a coordinate break
% (for the phase plate's lateral shift)
clear stratum
stratum.type = 3; % coordinate break
stratum.dx2 = dx2; % x2 shift of all strata above the coordinate break
stratum.dx3 = 0; % x3 shift of all strata above the coordinate break
grating.stratum{3} = stratum;
% Construct the stratum for the phase-plate grating.
clear stratum
stratum.type = 1; % uniperiodic stratum
stratum.thick = t3; % stratum thickness
% The following h11, h12 spec indicates that the stratum's period vector is
% half the grating period (i.e., its fundamental spatial frequency is twice
% that of the grating; see GD-Calc.pdf, Eq's 3.24 and 3.25).
stratum.h11 = 2;
stratum.h12 = 0;
clear stripe
stripe.c1 = -0.25; % first stripe's boundary on positive side
stripe.pmt_index = 1; % first stripe's permittivity index
stratum.stripe{1} = stripe;
stripe.c1 = 0.25; % second stripe's boundary on positive side
stripe.pmt_index = 3; % second stripe's permittivity index
stratum.stripe{2} = stripe;
grating.stratum{4} = stratum;
clear stratum stripe

% Enable full-field calculation.
grating.stratum{2}.full_field = [];

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength = wavelength;
inc_field.f2 = 0;
inc_field.f3 = 0;

% Specify which diffracted orders are to be retained in the calculations.
% (m2 is zero because this is a uniperiodic grating - all diffraction
% orders for m2~=0 have zero amplitude.)
clear order
order(1).m2 = 0;
order(1).m1 = -m_max:m_max;

% Run the diffraction calculations (with grating output to capture
% full-field data).
tic
[~,scat_field,inc_field,grating] = gdc(grating,inc_field,order,false);
toc

% Compute the diffraction efficiencies. (Only the reflected waves are
% retained.)
R = gdc_eff(scat_field,inc_field);
% Keep only orders -3...3.
R = R(abs([R.m1])<=3);
% Extract the diffraction order indices for the reflected waves.
m = [R.m1].';

% Tabulate diffraction efficiencies for unpolarized incident illumination,
% as a function of scan position (dx2) and order index (m).
disp(' ');
disp(['Scan position (dx2/period in col 1) and diffraction efficiency ' ...
    '(R in col''s 2...' num2str(1+length(m)) ')'])
disp(['for diffraction orders ' mat2str(m) ':'])
disp(num2str([[nan,m.'];[dx2,([R.eff1]+[R.eff2])/2]]))

% Generate an animation showing a dx2 parameter scan. Display the grating
% geometry in the top panel; and in the bottom panel display the H3p field
% amplitude at the air-gap top surface.
clear pmt_display
pmt_display(1).name = '';
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = ''; % reflector
pmt_display(2).color = [1,1,1]*0.75;
pmt_display(2).alpha = 1;
pmt_display(3).name = ''; % phase plate
pmt_display(3).color = [1,1,1]*0.875;
pmt_display(3).alpha = 1;
x_limit = [ ...
    -0.5*(t1+t2(1)+t3),-1.5*period,-1.5*period; ...
    1.5*(t1+t2(1)+t3),1.5*period,1.5*period];

x2 = (-1.5:.01:1.5)*period;
x3 = 0;
H3p = cell(size(dx2));
full_field = rmfield(grating.stratum{2}.full_field,{...
    'ffE1s','ffE2s','ffE3s','ffH1s','ffH2s','ffH3s',...
    'ffE1p','ffE2p','ffE3p','ffH1p','ffH2p'}); % keep only ffH3p
for p = 1:length(dx2)
    map = EH_map(full_field,p,x2,x3);
    H3p{p} = [map.H3p];
end

disp(' ')
disp('Making movie ...')

F = struct('cdata',{},'colormap',{}); % movie frames
fig = figure('Visible','off');
for p = 1:2:length(dx2)
    % Plot the grating with lateral shift dx2(p).
    % (The "false" argument makes the plot go to the current figure.)
    subplot(2,1,1);
    gdc_plot(grating,[p,1],pmt_display,x_limit,false);
    set(gca,...
        'XLabel',[],'YLabel',[],'ZLabel',[],...
        'XTick',[],'YTick',[],'ZTick',[])
    
    subplot(2,1,2);
    plot(x2,real(H3p{p}),'b',x2,imag(H3p{p}),'r');
    axis([floor(x2(1)),ceil(x2(end)),-2.5,2.5]);
    xlabel('x2');
    ylabel('H3p');
    legend('real(H3p)','imag(H3p)');
    title('Internal H3p field at top of air space')
    
    F(end+1) = getframe(fig); %#ok<SAGROW> 
end
if isempty(results_dir)
    set(fig,'Visible','on')
    movie(fig,F,5);
else
    vidObj = VideoWriter([results_dir '/gdc_demo13']);
    vidObj.FrameRate = 10;
    open(vidObj);
    for count = 1:5
        writeVideo(vidObj,F);
    end
    close(vidObj)
    close(fig)
    disp(['Video file ' results_dir '/gdc_demo13.avi has been created.'])
end
