%start timer
tic

% add path to GD-Calc installation
addpath('./GD-Calc/code/')

% define grating parameters
PGM.N=600; % grating line density in lines/mm
alpha=0.49; % blaze angle in degrees measured relative to surface
apex = 175.62; % apex angle in degrees
num_slices=5; % number of slices defined to approximate %
                % sloping groove profile using staircase approximation
num_periods=1; % number of multilayer periods 
                % (use 1 for single-layer coating)
d=30; % multilayer d-spacing in nm (or thickness of the single-layer)
gamma_d=1; % ratio of 1st layer thickness in multilayer coating to the 
            % 2nd layer thickness in multilayer (use 1 for single-layer
            % coating)
            
% define PGM parameters
energy = 125:1000:10000; % in eV
PGM.cff = 2;
PGM.order = 1; % this defines the simulated diffraction order efficiency
PGM.energy=energy(1); % sets PGM.energy to the first element in 'energy'
PGM = pgm(PGM); % redefines PGM object as an instance of a pgm class

% define location of optical constants files
datadir='./RefInds/';

% define materials that make up grating. The first element is the
% substrate material, then the two materials that make up the multilayer
% for a single layer coating only two materials need to be defined (the
% substrate material followed by the coating material)
materials={};
% materials = {'Si','C','Cr'};
materials = {'Si','Pt'};

% convert optical constants from CXRO website to format required by GD-Calc
% (this does not need to be run every time simulations are performed)
for i=numel(materials)
    MLgrating_write_nk(datadir,materials{i});
end

% define maximum diffraction order that is simulated
m_max = 17; % equivalent to REFLEC maximum

% find number of energies
numE=numel(energy);

%initialise grating
grating={};

% initialise efficiency
eff = zeros(4,numE);

% initialise record of time elapsed
def_time=zeros(1,numE);
sim_time=zeros(1,numE);

% loop over energies of interest
for i=1:numE
    % set PGM object so that it only contains 1 energy
    PGM.energy=energy(i);
    % define grating
    [grating,def_time(i)] = MLgrating_def_blaz(PGM.N,alpha,apex,num_slices,d,gamma_d,num_periods,materials,datadir,PGM.wavelength);
    % perform grating efficiencies simulation
    [eff(1,i),eff(2,i),eff(3,i),eff(4,i),sim_time(i)] = MLgrating_sim(grating,PGM,m_max);
end

% close open figures if there are any
close all

% plot colour representation of grating object provided to GD-Calc
MLgrating_plot(grating,materials)

% Plot simulated efficiencies
figure('color','w')
plot(energy,eff(1,:))
hold all
plot(energy,eff(2,:))
% plot(energy,eff(3,:))
% plot(energy,eff(4,:))
xlabel('Energy (eV)')
ylabel('-1st order grating efficiency')
grid on
legend('sigma','pi')

% end timer
toc