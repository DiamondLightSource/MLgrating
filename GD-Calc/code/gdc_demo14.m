% gdc_demo14: biperiodic grating with internal field calculation
%
% Biperiodic gold grating comprising square pillars
% period: 1 micron
% depth: 0.5 micron
% wavelength: 0.6328 micron
% Au permittivity (grating): (0.14330+3.6080i)^2
% SiO2 permittivity (substrate): 1.4585^2
% normal incidence
% Plot horizontal cross-sections of |E1|, |E2|, |E3|, |H1|, |H2|, |H3| for
% incident E field polarized parallel to x3 axis.
%
% Refractive index data: https://refractiveindex.info/
%   Shelf: MAIN, Book: Au (Gold), Page: McPeak et al. 2015
%   Shelf: MAIN, Book: SiO2, Page: Malitson 1965
%
% Documentation references:
%   GD-Calc_Demo.pdf and GD-Calc.pdf
%   gdc.m, gdc_plot.m, gdc_eff.m (comment headers)
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% If results_dir is nonempty the computation results will be saved to
% [results_dir '/gdc_demo14.mat']. Also, a video file will be saved to
% [results_dir '/gdc_demo14.avi'].
results_dir = '../results';
if ~isempty(results_dir) && ~exist(results_dir,'dir')
    warning(['Results directory ''' results_dir ''' does not exist.'])
    results_dir = [];
end

disp(' ')
disp('gdc_demo14.m')

% The "if" branch can be toggled to regnerate the plots from previously
% computed data.
if 1
    % Define parameters for grating structure and incident field:
    wavelength = 0.6328;
    grating_pmt = (0.14330+3.6080i)^2; % grating permittivity
    substrate_pmt = 1.4585^2; % substrate permittivity
    h = .5; % grating height
    d = 1; % grating period
    w = 0.25; % pillar width
    m_max = 12; % maximum diffraction order index
    N_plot = 10; % number of full-field samples inside grating

    % Construct grating.
    clear grating
    grating.pmt = {1,1,grating_pmt,substrate_pmt}; ...
        % grating material permittivities
    % (Two vacuum permittivities, pmt{1} and pmt{2}, are used to facilitate
    % plotting.) 
    grating.pmt_sub_index = 4; % substrate permittivity index
    grating.pmt_sup_index = 1; % superstrate permittivity index
    % Define the x2 and x3 projections of the first grating period
    % (d21,d31) and second grating period (d22,d32). The second period is
    % parallel to the x2 axis, and the first period is parallel to the x3
    % axis.
    grating.d21 = 0;
    grating.d31 = d;
    grating.d22 = d;
    grating.d32 = 0;
    grating.stratum = {};
    clear stratum
    % Add strata for sampling the E and H fields in the substrate.
    stratum.type = 0;
    stratum.pmt_index = 4;
    stratum.thick = h/N_plot;
    stratum.full_field = [];
    for l1 = 1:N_plot
        grating.stratum{end+1} = stratum;
    end
    % Grating layer:
    clear stratum
    stratum.type = 2; % biperiodic stratum
    stratum.thick = h; % stratum thickness
    % The following h11 ... h22 spec defines the stratum's period vectors
    % (GD-Calc.pdf, equation 3.20).
    stratum.h11 = 1;
    stratum.h12 = 0;
    stratum.h21 = 0;
    stratum.h22 = 1;
    clear stripe
    stripe.type = 0; % homogeneous stripe
    stripe.c1 = -w/2;
    stripe.pmt_index = 1;
    stratum.stripe{1} = stripe;
    clear stripe block
    stripe.type = 1; % inhomogeneous stripe
    stripe.c1 = w/2;
    block.c2 = -w/2;
    block.pmt_index = 1;
    stripe.block{1} = block;
    block.c2 = w/2;
    block.pmt_index = 3;
    stripe.block{2} = block;
    stratum.stripe{2} = stripe;
    clear stripe block
    stratum.thick = stratum.thick/N_plot;
    stratum.full_field = [];
    for l1 = 1:N_plot
        grating.stratum{end+1} = stratum;
    end
    % Add strata for sampling the E and H fields in the superstrate.
    clear stratum
    stratum.type = 0;
    stratum.pmt_index = 2; % not 1 to facilitate plotting
    stratum.thick = h/N_plot;
    stratum.full_field = [];
    for l1 = 1:N_plot
        grating.stratum{end+1} = stratum;
    end
    clear stratum

    % Define the indicent field (normal incidence).
    clear inc_field
    inc_field.wavelength = wavelength;
    inc_field.f2 = 0;
    inc_field.f3 = 0;

    % Specify which diffracted orders are to be retained in the
    % calculations.
    order = [];
    for m2 = -m_max:m_max
        order(end+1).m2 = m2; %#ok<SAGROW> 
        order(end).m1 = -m_max:m_max;
    end

    % Run the diffraction calculations (with grating output to capture
    % full-field data).
    tic
    [param_size,scat_field,inc_field,grating] = ...
        gdc(grating,inc_field,order,false);
    toc
    % Note: inc_field.p2 = 0, inc_field.p3 = 1; "p" polarization is with
    % the incident field polarized in the e3 direction.

    % Compute the diffraction efficiencies.
    [R,T] = gdc_eff(scat_field,inc_field);
    % Discard diffracted waves that decay exponentially with distance from
    % the grating. These include evanescent waves and, if the substrate's
    % permittivity is not real-valued, all transmitted waves.
    R = R(imag([scat_field.f1r])==0);
    T = T(imag([scat_field.f1t])==0);
    % Tabulate the diffraction order indices and diffraction efficiencies
    % for an incident field polarized in the e3 direction.
    disp(' ');
    disp('Diffraction efficiencies (m1, m2, eff2)');
    disp('R:');
    disp(num2str([[R.m1].' [R.m2].' [R.eff2].']));
    if ~isempty(T)
        disp('T:');
        disp(num2str([[T.m1].' [T.m2].' [T.eff2].']));
    end
    
    % Calculate spatial field maps.
    x2 = (-d/2:.01:d/2);
    x3 = (-d/2:.01:d/2).';
    size_x = max(size(x2),size(x3));
    x2 = repmat(x2,size_x./size(x2));
    x3 = repmat(x3,size_x./size(x3));
    x1 = [];
    E1p = zeros([size_x 0]);
    E2p = zeros([size_x 0]);
    E3p = zeros([size_x 0]);
    H1p = zeros([size_x 0]);
    H2p = zeros([size_x 0]);
    H3p = zeros([size_x 0]);
    disp(' ')
    disp('Doing Fourier transforms ...')
    tic
    x1_ = 0;
    for l1 = 1:length(grating.stratum)
        if isfield(grating.stratum{l1},'thick')
            x1_ = x1_+grating.stratum{l1}.thick;
        end
        if isfield(grating.stratum{l1},'full_field')
            x1(end+1) = x1_-h; ...
                %#ok<SAGROW> % "-h" to make x1 = 0 at base of grating
            map = EH_map(rmfield(grating.stratum{l1}.full_field,{...
                'ffE1s','ffE2s','ffE3s','ffH1s','ffH2s','ffH3s'...
                }),1,x2,x3);
            E1p(:,:,end+1) = reshape([map.E1p],size_x); %#ok<SAGROW> 
            E2p(:,:,end+1) = reshape([map.E2p],size_x); %#ok<SAGROW>
            E3p(:,:,end+1) = reshape([map.E3p],size_x); %#ok<SAGROW>
            H1p(:,:,end+1) = reshape([map.H1p],size_x); %#ok<SAGROW>
            H2p(:,:,end+1) = reshape([map.H2p],size_x); %#ok<SAGROW>
            H3p(:,:,end+1) = reshape([map.H3p],size_x); %#ok<SAGROW>
        end
    end
    toc

    if ~isempty(results_dir)
        disp(' ')
        disp('Saving ...')
        tic
        save([results_dir '/gdc_demo14'])
        toc
    end

else
    assert(~isempty(results_dir),'results_dir must be specified.')
    disp(' ')
    disp('Loading ...')
    tic
    load([results_dir '/gdc_demo14'])
    toc
end

% Plot the grating.
clear pmt_display
pmt_display(1).name = ''; % Vacuum
pmt_display(1).color = [];
pmt_display(1).alpha = 1;
pmt_display(2).name = 'Vac'; % Vacuum (field-sampling strata)
pmt_display(2).color = [0,0,1;0,0,0];
pmt_display(2).alpha = [0.01;0.5];
pmt_display(3).name = 'Au'; % Au
pmt_display(3).color = [1,.85,0;0,0,0];
pmt_display(3).alpha = 1;
pmt_display(4).name = 'SiO2'; % SiO2
pmt_display(4).color = [0,1,0;0,0,0];
pmt_display(4).alpha = [0.1;0.5];
x_limit = [-0.5*h,-1.5*d,-1.3*d;3.5*h,1.5*d,1.3*d];
gdc_plot(grating,1,pmt_display,x_limit);
% Move x1 origin to base of grating.
set(gca,'ZTickLabel', ...
    cellfun(@num2str, ...
    num2cell( ...
    cellfun(@str2num, ...
    get(gca,'ZTickLabel'))-h),'UniformOutput',false));

disp(' ')
disp('Making movie ...')
tic

% Note: Without the drawnow only the first movie frame goes to fig; the
% other frames go to the preceding gdc_plot figure. (MATLAB bug.)
drawnow
F = struct('cdata',{},'colormap',{}); % movie frames
fig = figure('Visible','off');
fig.Position = [100,100,850,425];
maxE1p = max(abs(E1p(:)));
maxE2p = max(abs(E2p(:)));
maxE3p = max(abs(E3p(:)));
maxH1p = max(abs(H1p(:)));
maxH2p = max(abs(H2p(:)));
maxH3p = max(abs(H3p(:)));
for l1 = 1:length(x1)
    x1_ = x1(l1);
    str = sprintf('%5.2f',x1_);
    if x1_<=0
        str = ['x1=' str ' (in SiO2)      '];
    elseif x1_<=h
        str = ['x1=' str ' (in Au grating)'];
    else
        str = ['x1=' str ' (in Vac)       '];
    end
    subplot(2,3,1)
    cla(gca,'reset'), set(gca,'TitleHorizontalAlignment','left')
    surface(x2,x3,abs(E1p(:,:,l1)),'EdgeAlpha',.25,'FaceColor','interp')
    set(gca,'ZLim',[0,maxE1p])
    view(-45,60);
    xlabel('x2');
    ylabel('x3');
    zlabel('|E1|');
    title(['|E1| ' str]);
    subplot(2,3,2)
    cla(gca,'reset'), set(gca,'TitleHorizontalAlignment','left')
    surface(x2,x3,abs(E2p(:,:,l1)),'EdgeAlpha',.25,'FaceColor','interp')
    set(gca,'ZLim',[0,maxE2p])
    view(-45,60);
    xlabel('x2');
    ylabel('x3');
    zlabel('|E2|');
    title(['|E2| ' str]);
    subplot(2,3,3)
    cla(gca,'reset'), set(gca,'TitleHorizontalAlignment','left')
    surface(x2,x3,abs(E3p(:,:,l1)),'EdgeAlpha',.25,'FaceColor','interp')
    set(gca,'ZLim',[0,maxE3p])
    view(-45,60);
    xlabel('x2');
    ylabel('x3');
    zlabel('|E3|');
    title(['|E3| ' str]);
    subplot(2,3,4)
    cla(gca,'reset'), set(gca,'TitleHorizontalAlignment','left')
    surface(x2,x3,abs(H1p(:,:,l1)),'EdgeAlpha',.25,'FaceColor','interp')
    set(gca,'ZLim',[0,maxH1p])
    view(-45,60);
    xlabel('x2');
    ylabel('x3');
    zlabel('|H1|');
    title(['|H1| ' str]);
    subplot(2,3,5)
    cla(gca,'reset'), set(gca,'TitleHorizontalAlignment','left')
    surface(x2,x3,abs(H2p(:,:,l1)),'EdgeAlpha',.25,'FaceColor','interp')
    set(gca,'ZLim',[0,maxH2p])
    view(-45,60);
    xlabel('x2');
    ylabel('x3');
    zlabel('|H2|');
    title(['|H2| ' str]);
    subplot(2,3,6)
    cla(gca,'reset'), set(gca,'TitleHorizontalAlignment','left')
    surface(x2,x3,abs(H3p(:,:,l1)),'EdgeAlpha',.25,'FaceColor','interp')
    set(gca,'ZLim',[0,maxH3p])
    view(-45,60);
    xlabel('x2');
    ylabel('x3');
    zlabel('|H3|');
    title(['|H3| ' str]);
    
    F(end+1) = getframe(fig); %#ok<SAGROW> 
end
Position = get(fig,'Position');
close(fig)
toc
if isempty(results_dir)
    movie(figure('Position',Position),F,1,2);
else
    vidObj = VideoWriter([results_dir '/gdc_demo14']);
    vidObj.FrameRate = 2;
    open(vidObj);
    writeVideo(vidObj,F);
    close(vidObj)
    disp(['Video file ' results_dir '/gdc_demo14.avi has been created.'])
end
