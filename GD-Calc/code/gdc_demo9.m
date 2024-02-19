% gdc_demo9: alignment sensor
%
% The grating structure in this demo comprises a phase-plate transmission
% grating in close proximity to a reflection grating. (Both elements are
% uniperiodic, lamellar gratings.) The lateral dislacement of the phase
% plate and the air gap between the gratings are both vectorized
% parameters. The ratio of the +1 to -1 diffraction order efficiencies
% provides a sensitive measure of the lateral displacement.
%
% This demo illustrates animation of the lateral grating scan using the
% MATLAB movie function. The movie can optionally be saved to a video file
% gdc_demo9.avi in directory results_dir. (Set results_dir = [] to disable
% saving.)
%
% This demo also illustrates data writing to a spreadsheet file
% gdc_demo9.xls.xls in directory results_dir (if results_dir is nonempty).
%
% gdc_demo9 covers the following topics:
%	- multilayer grating (uniperiodic and homogeneous layers)
%   - harmonic indices
%   - parameterization
%   - coordinate break
%   - animation
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% If results_dir is nonempty a video of the lateral grating scan will be
% saved to [results_dir '/gdc_demo9.avi']. Also, the results data will be
% written to [results_dir '/gdc_demo9.xls'].
results_dir = '../results';
if ~isempty(results_dir) && ~exist(results_dir,'dir')
    warning(['Results directory ''' results_dir ''' does not exist.'])
    results_dir = [];
end

disp(' ')
disp('gdc_demo9.m')

% Define parameters for grating structure and incident field. The phase
% plate's lateral shift parameter (dx2) is vectorized in the first
% dimension, and the air space parameter (t2) is vectorized in the second
% dimension.
superstrate_pmt = 2.25; % superstrate (phase plate) permittivity
substrate_pmt = (1.37+7.62i)^2; % substrate permittivity (Aluminum)
lambda = 0.633; % wavelength, micron
t1 = lambda/4; % bottom (reflecting) grating thickness
t2 = lambda*[2,2.5,3]; % air space thickness (vectorized in dimension 2)
t3 = (lambda/8)/(sqrt(superstrate_pmt)-1); % phase plate thickness
period = 10*lambda; % grating period
dx2 = period*(0:63).'/64; ...
    % phase plate's lateral shift (vectorized in dimension 1)
m_max = 20; % maximum diffraction order index

% Construct grating.
clear grating
grating.pmt = {1,substrate_pmt,superstrate_pmt}; % material permittivities
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

% Define the indicent field (normal incidence).
clear inc_field
inc_field.wavelength = lambda;
inc_field.f2 = 0;
inc_field.f3 = 0;

% Specify which diffracted orders are to be retained in the calculations.
% (m2 is zero because this is a uniperiodic grating - all diffraction
% orders for m2~=0 have zero amplitude.)
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
% Keep only orders -1:1.
R = R(abs([R.m1])<=1);
% Extract the diffraction order indices for the reflected waves.
m = [R.m1].';
% Construct strings for plot legend.
m_str = {};
for k = 1:length(m)
    m_str{k} = ['order ' num2str(m(k))]; %#ok<SAGROW> 
end
% Take the average of eff1 and eff2 to get the efficiency for unpolarized
% illumination.
R = 0.5*(cat(3,R.eff1)+cat(3,R.eff2));
% R is 3-D with dim's 1, 2, 3 corresponding to dx2, t2, m.

% Tabulate the efficiency data.
disp(['Scan position (dx2/period in col 1) and diffraction efficiency ' ...
    '(R in col''s 2...' num2str(1+length(t2)) ')'])
disp(['for air gap t2 = ' mat2str(t2/lambda) '*lambda:'])
for k = 1:length(m)
    disp([m_str{k} ' reflection efficiencies:'])
    disp(num2str([dx2/period,R(:,:,k)]))
end

if ~isempty(results_dir)
    % Save the results to an xls file.
    xls_file = [[results_dir '/gdc_demo9.xls'] '.xls'];
    onoff = warning('off','MATLAB:xlswrite:AddSheet');
    for k = 1:length(m)
        writecell(m_str(k), ...
            xls_file,'Sheet',k,'Range','A1:A1')
        writecell({'dx2/period'},xls_file,'Sheet',k,'Range','B2:B2')
        writecell({'t2/lambda'},xls_file,'Sheet',k,'Range','A3:A3')
        writematrix(t2/lambda,xls_file,'Sheet',k, ...
            'Range',['C3:' char('C'+length(t2)-1) '3'])
        writematrix([dx2/period,R(:,:,k)],xls_file,'Sheet',k, ...
            'Range',['B4:' char('B'+length(t2)) num2str(3+length(dx2))])
    end
    warning(onoff.state,'MATLAB:xlswrite:AddSheet');
end

% Plot the efficiency data.
figure, hold on
c = 'bgr'; % plot colors for m
c_ = {':','-','--'}; % plot line style for t2
legend_ = {};
for im = 1:length(m)
    for it2 = 1:length(t2)
        plot(dx2/period,R(:,it2,im),[c(im),c_{it2}])
        legend_{end+1} = ['order ' num2str(m(im)) ...
            ', t2=' num2str(t2(it2)/lambda) '\lambda']; %#ok<SAGROW> 
    end
end
title('Efficiency vs dx2');
legend(legend_);
xlabel('scan position dx2/period');
ylabel('efficiency');
axis([0,1.5,0,0.65]);
set(gcf,'Position',[300,400,560,200])
drawnow

% Generate an animation showing a dx2 parameter scan.
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
waitbar_ = waitbar(0,'Making movie, please wait ...');
F = struct('cdata',{},'colormap',{});
fig = figure('Visible','off');
title('Alignment sensor dx2 scan')
p2 = 2;
for p1 = 1:4:length(dx2)
    % Plot the grating with air space t2(p1) and lateral shift dx2(p2).
    % (The "false" argument makes the plot go to the current figure.)
    gdc_plot(grating,[p1,p2],pmt_display,x_limit,false);
    set(gca,...
        'XLabel',[],'YLabel',[],'ZLabel',[],...
        'XTick',[],'YTick',[],'ZTick',[])
    legend(gca,'off')
    F(end+1) = getframe(fig);    %#ok<SAGROW> 
    waitbar(p1/length(dx2),waitbar_);
end
close(waitbar_)
if isempty(results_dir)
    set(fig,'Visible','on')
    movie(fig,F,5);
else
    % Save the video.
    vidObj = VideoWriter([results_dir '/gdc_demo9.avi']);
    vidObj.FrameRate = 10;
    open(vidObj);
    for count = 1:5
        writeVideo(vidObj,F);
    end
    close(vidObj)
    close(fig)
end

% For static display:
% gdc_plot(grating,[12,2],pmt_display,x_limit);
