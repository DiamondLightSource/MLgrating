%start timer
tic

% add path to GD-Calc installation
addpath('GD-Calc\code\')

% define grating parameters
PGM.N=400; % grating line density in lines/mm
A.h=13.2; % groove height in nm
A.T=14; % trapezoidal angle in degrees measured relative to surface
A.gamma_D=0.64; % groove width / groove d-spacing (groove measured on top)
A.num_slices=3; % number of slices defined to approximate %
                % sloping groove profile using staircase approximation
% A.num_slices=10; % number of slices defined to approximate %
%                 % sloping groove profile using staircase approximation
A.num_periods=1; % number of multilayer periods 
                % (use 1 for single-layer coating)
A.d=40; % multilayer d-spacing in nm (or thickness of the single-layer)
A.gamma_d=1; % ratio of 1st layer thickness in multilayer coating to the 
            % 2nd layer thickness in multilayer (use 1 for single-layer
            % coating)
        
            
% define PGM parameters
% energy = linspace(50,1500,5); % in eV
% energy = 50:10:15000; % in eV
% energy = 50:100:15000; % in eV
% energy = logspace(log10(50),log10(15000),10); % in eV
energy = logspace(log10(50),log10(15000),512); % in eV
PGM.energy=energy(1); % sets PGM.energy to the first element in 'energy'

% cff = 1.4:0.6:2;
cff = 1.2:0.2:3;
PGM.cff=cff(1); % sets PGM.energy to the first element in 'energy'

% order=1:3;
order=1:5;
PGM.order = order(1); % this defines the simulated diffraction order efficiency
PGM = pgm(PGM); % redefines PGM object as an instance of a pgm class

% define location of optical constants files
datadir='./RefInds/';

% define materials that make up grating. The first element is the
% substrate material, then the two materials that make up the multilayer
% for a single layer coating only two materials need to be defined (the
% substrate material followed by the coating material)
A.materials={};
% materials = {'Si','C','Cr'};
A.materials = {'Si','Au'};

% convert optical constants from CXRO website to format required by GD-Calc
% (this does not need to be run every time simulations are performed)
% for i=1:numel(A.materials)
%     MLgrating_write_nk(datadir,A.materials{i});
% end

% define maximum diffraction order that is simulated
A.m_max = 17; % equivalent to REFLEC maximum
% A.m_max = 3; 

% find number of energies
numE=numel(energy);

% find number of cffs
numcff=numel(cff);

% find number of orders
numorder=numel(order);

%initialise grating
grating={};

% initialise efficiency
eff = zeros(4,numE);

% loop over orders of interest
for k=1:numorder
    % loop over cffs of interest
    for j=1:numcff
        % loop over energies of interest
        for i=1:numE
            % set PGM object so that it only contains 1 energy
            PGM.energy=energy(i);
            % set PGM object so that it only contains 1 energy
            PGM.order=order(k);
            A.order(k).order=-order(k);
            % set PGM object so that it only contains 1 energy
            PGM.cff=cff(j);
            A.order(k).cff(j).cff=cff(j);            
            % define grating
            [grating,def_time(i,j,k)] = MLgrating_def_lam(PGM.N,A.h,A.T,A.gamma_D,A.num_slices,A.d,A.gamma_d,A.num_periods,A.materials,datadir,PGM.wavelength);
            disp(['Simulating at ' num2str(energy(i)) ...
            ' eV, c_ff = ' num2str(cff(j)) ...
            ' and m = '  num2str(order(k)) ' ...'])
            % perform grating efficiencies simulation
            [eff_s,eff_p,eff_d,eff_c,sim_time] = MLgrating_sim(grating,PGM,A.m_max);
            disp(['Elapsed time was ' num2str(sim_time) ' s'])
            A.order(k).cff(j).timestamp{i}=datestr(now);
            A.order(k).cff(j).energy(i)=energy(i);
            A.order(k).cff(j).s(i)=eff_s;
            A.order(k).cff(j).p(i)=eff_p;
            A.order(k).cff(j).d(i)=eff_d;
            A.order(k).cff(j).c(i)=eff_c;
            A.order(k).cff(j).sim_time(i)=sim_time;
        end
%         eval(['A.o' num2str(order(k)) '.cff{' num2str(j) '}.s=eff_s;']);
    end
end

% close open figures if there are any
% close all

% % plot colour representation of grating object provided to GD-Calc
% MLgrating_plot(grating,materials)

txt=jsonencode(A,PrettyPrint=true);
fid = fopen('file.json','w');
fprintf(fid,'%s',txt);

fclose(fid);
% 
% % Plot simulated efficiencies
% figure('color','w')
% plot(energy,eff(1,:,1,1))
% hold all
% plot(energy,eff(1,:,1,2))
% plot(energy,eff(1,:,1,3))
% xlabel('Energy (eV)')
% ylabel('Efficiency')
% grid on
% legend('m=-1','m=-2','m=-3')

% end timer
toc