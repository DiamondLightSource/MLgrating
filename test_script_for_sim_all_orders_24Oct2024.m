%start timer
tic

% add path to GD-Calc installation
addpath('./GD-Calc/code/')

% define grating parameters
PGM.N=600; % grating line density in lines/mm
A.alpha=0.49; % blaze angle in degrees measured relative to surface
A.apex=175.62; % apex angle in degrees
A.num_slices=10; % number of slices defined to approximate %
                % sloping groove profile using staircase approximation
A.num_periods=1; % number of multilayer periods 
                % (use 1 for single-layer coating)
A.d=30; % multilayer d-spacing in nm (or thickness of the single-layer)
A.gamma_d=1; % ratio of 1st layer thickness in multilayer coating to the 
            % 2nd layer thickness in multilayer (use 1 for single-layer
            % coating)
            
% define PGM parameters
energy = [250:250:3500]; % in eV
PGM.energy=energy(1); % sets PGM.energy to the first element in 'energy'

cff = 3.5;
PGM.cff=cff(1); % sets PGM.energy to the first element in 'energy'

order=1;
PGM.order = order(1); % this defines the simulated diffraction order efficiency

PGM = pgm(PGM); % redefines PGM object as an instance of a pgm class

% define location of optical constants files
datadir=['C:\Users\aw63\OneDrive - Diamond Light Source Ltd\' ...
'Optics\Diamond II\Multilayer gratings\MLgrating\RefInds\'];

% define materials that make up grating. The first element is the
% substrate material, then the two materials that make up the multilayer
% for a single layer coating only two materials need to be defined (the
% substrate material followed by the coating material)
A.materials={};
A.materials = {'Si','Au'};

% convert optical constants from CXRO website to format required by GD-Calc
% (this does not need to be run every time simulations are performed)
% for i=1:numel(materials)
%     MLgrating_write_nk(datadir,materials{i});
% end

% define maximum diffraction order that is simulated
A.m_max = 17; % equivalent to REFLEC maximum

% find number of energies
numE=numel(energy);

% find number of cffs
numcff=numel(cff);

% find number of orders
numorder=numel(order);

%initialise grating
grating={};

% loop over energies of interest
for i=1:numE
    % set PGM object so that it only contains 1 energy value
    PGM.energy=energy(i);
    % define grating
    [grating,def_time(i)] = ...
        MLgrating_def_blaz(PGM.N,A.alpha,A.apex, ...
        A.num_slices,A.d,A.gamma_d,A.num_periods, ...
        A.materials,datadir,PGM.wavelength);
    % perform grating efficiencies simulation
    [eff_s(:,i),eff_p(:,i),eff_d(:,i),eff_c(:,i),output_orders] ...
        = MLgrating_sim_all_orders(grating,PGM,A.m_max);
end

% close open figures if there are any
% close all

% plot colour representation of grating object provided to GD-Calc
% MLgrating_plot(grating,A.materials)

% Example plotting of different diffraction order efficiencies
figure('color','w')
plot(energy,eff_s(0==output_orders,:))
hold all
plot(energy,eff_s(-1==output_orders,:))
xlabel('Energy (eV)')
ylabel('Grating efficiencies')
legend('order = 0','order = -1')
grid on