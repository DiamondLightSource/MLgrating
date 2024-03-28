function [grating,time_elapsed] = MLgrating_def_blaz(N,alpha,apex, ...
    num_slices,d,gamma_d,num_periods,materials,datadir,wavelength)

% all lengths (including wavelengths) are defined in nm within code

% start timer
tStart = tic;

% Grating groove period D is calculated in nm (needs line density to be
% provided in lines/mm)
D = 1e6/N;

% Construct grating. (Note: Two grating periods must be specified
% because the grating can generally be biperiodic, although for this
% example the second period (d22,d32) is irrelevant.)

grating.pmt_sub_index = 2; % substrate permittivity index
grating.pmt_sup_index = 1; % superstrate permittivity index

% defines uniperiodic grating with period of D
grating.d21 = D; % first grating period: x2 projection
grating.d31 = 0; % first grating period: x3 projection
grating.d22 = 0; % second grating period: x2 projection
grating.d32 = D; % second grating period: x3 projection

% initialise grating stratum
grating.stratum = {};

stratum.type = 1; % uniperiodic stratum

% The following h11, h12 spec indicates that the stratum's period
% vector matches the first grating period (GD-Calc.pdf, Eq's 3.23-25).
stratum.h11 = 1;
stratum.h12 = 0;

% calculate antiblaze angle from blaze angle and apex angle
antiblaze=180-alpha-apex;

% calculate tangent of blaze and antiblaze angles once to avoid rounding
% errors
tand_alpha=tand(alpha);
tand_antiblaze=tand(antiblaze);

% calculate grating height
h = D./((1./tand_alpha)+(1./tand_antiblaze)); % grating height

% calculate horizontal width of blaze portion of profile in units of D
x = h./(tand_alpha.*D);

% if alpha is defined too small so grating height is less than h then throw
% exception
if x > D/2
    throw('Value of alpha is too small!')
end

% small correction to the multilayer period to ensure it is defined 
% perpendicular to the blazed surface
d_sinalpha=d*sind(alpha);

% calculate arrays containing heights defining multilayer structure (where
% a height of zero corresponds to the top surface of the substrate)
[ML_thickness,H,h_bottoms,h_tops,heights,thicknesses] = ...
    MLgrating_def_coat(h,num_slices,d_sinalpha,gamma_d,num_periods);
    
% find total number of strata
num_strata=numel(heights);

m=1;
for j=1:num_strata
    stratum.thick = thicknesses(j); % stratum thickness
    % positive gradient section
    stripe_max=sum(heights(j) > h_bottoms);
    stripe_min=sum(heights(j)-h > h_bottoms)+1;
    i=1;
    for k = stripe_max:-1:stripe_min
        stripe.c1 = -0.5+(heights(j)-h_bottoms(k)).*x/h; % stripe's boundary on positive side
        if isequal(k,numel(h_bottoms))
            stripe.pmt_index = 1; % Vacuum's permittivity index
        elseif mod(k,2)
            stripe.pmt_index = 3; % 1st layer's permittivity index
        else
            stripe.pmt_index = 4; % 2nd layer's permittivity index
        end
        stratum.stripe{i} = stripe;
        stripe={};
        i=i+1;
    end
    % negative gradient section
    for k = stripe_min:stripe_max
        stripe.c1 = +0.5-(heights(j)-h_bottoms(k)).*(1-x)/h; % stripe's boundary on positive side
        if isequal(k,1)
            stripe.pmt_index = 2; % substrate's permittivity index
        elseif mod(k,2)
            stripe.pmt_index = 4; % 2nd layer's permittivity index
        else
            stripe.pmt_index = 3; % 1st layer's permittivity index
        end
        stratum.stripe{i} = stripe;
        stripe={};
        i=i+1;
    end
    % grating top (handles cases where there is only one stripe in a
    % stratum)
    if ~isequal(i,1)
        if heights(j) >= h_tops(m+1)
            m=m+1;
        end
    elseif isequal(i,1)
        if heights(j) < h
            stripe.pmt_index = 2; % substrate's permittivity index      
        elseif heights(j) >= H
            stripe.pmt_index = 1; % Vacuum's permittivity index
        elseif m == length(h_tops)
            stripe.pmt_index = 1; % Vacuum's permittivity index
        elseif heights(j) < h_tops(m+1)
            if mod(m,2)
                stripe.pmt_index = 3; % 1st layer's permittivity index
            else
                stripe.pmt_index = 4; % 2nd layer's permittivity index
            end
        elseif heights(j) >= h_tops(m+1)
            m=m+1;
            if mod(m,2)
                stripe.pmt_index = 3; % 1st layer's permittivity index
            else
                stripe.pmt_index = 4; % 2nd layer's permittivity index
            end
        end
        stripe.c1 = 0.5; % second stripe's boundary on positive side    
        stratum.stripe{1} = stripe;
        stripe={};
    end
    grating.stratum{end+1} = stratum;
    stratum.stripe={};
end

% read in and define permittivities for all grating materials
grating.pmt = MLgrating_def_perm(datadir,materials,wavelength);

time_elapsed=toc(tStart);