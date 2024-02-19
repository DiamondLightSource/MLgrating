function [ML_thickness,H,h_bottoms,h_tops,heights,thicknesses] = ...
    MLgrating_def_coat(h,num_slices,d,gamma_d,num_periods)

% Carefully define parameters for creating ML grating once to try to
% minimise rounding errors
ML_thickness=num_periods*d;
H=h+ML_thickness;
gamma_dxd=gamma_d*d;
rounding_error=1e-9*d;

% define h_bottoms: the heights of the interfaces within the grooves of the
% laminar grating (zero corresponds to the surface of the substrate)
h_bottoms_1stlayers=0:d:ML_thickness+rounding_error;
h_bottoms_2ndlayers=gamma_dxd:d:ML_thickness+rounding_error;
h_bottoms=sort([h_bottoms_1stlayers h_bottoms_2ndlayers]);

% define h_tops: the heights of the interfaces outside the grooves of the
% laminar grating (zero corresponds to the surface of the substrate)
h_tops=h_bottoms+h;

% % Define heights of all interfaces required to define laminar grating
% h_total=sort([h_bottoms h_tops]);

% approximate stratum thickness
t = h/num_slices; 

if ML_thickness>=h
    % define interfaces for first ML period
    h_strata_bottoms_1stperiod=[0:t:gamma_dxd gamma_dxd:t:d];

    % use number of interfaces in first ML period
    A=numel(h_strata_bottoms_1stperiod);

    % define h_bottoms for all ML periods
    h_strata_bottoms=zeros(1,num_periods*A);
    for i=1:num_periods
        h_strata_bottoms((i-1)*A+1:i*A)=(i-1)*d+h_strata_bottoms_1stperiod;
    end

    % define h_tops: the heights of the interfaces outside the grooves of the
    % laminar grating (zero corresponds to the surface of the substrate)
    h_strata_tops=h_strata_bottoms+h;

    % Define heights of all interfaces required to define laminar grating
    h_strata_total=sort([h_strata_bottoms h_strata_tops]);

    % define thicknesses of stripes required to define laminar grating
    thicknesses=diff(h_strata_total);

    % define average height of each stratum
    heights=h_strata_total(1:end-1)+0.5*thicknesses;

else
    % define average height of each stratum
    heights=0.5*t:t:H;
    
    % define average height of each stratum
    thicknesses=ones(1,numel(heights))*t;
end