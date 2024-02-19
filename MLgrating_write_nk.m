function wnk = MLgrating_write_nk(datadir,materialname)
% This function converts refractive indices
% exported from the CXRO website
% (https://henke.lbl.gov/optical_constants/getdb2.html)
% and writes a new file in a format that can be imported by read_nk, 
% a function intrinsic to GD-calc. The input filename containing
% the CXRO data should be in the format 
%
% materialname 'Henke_delta_beta.txt'
%
% and the output filename has the format
% 
% [materialname 'Henke_n_k.txt']
%
% These filenames can be modified by modifying this function or 
% MLgrating_def_perm if preferred

w1=importdata([datadir materialname '_Henke_delta_beta.txt']);

% The energy (in eV), as well as delta and beta are extracted from the
% input file. 
PHOTON.energy=w1.data(:,1);
delta=w1.data(:,2);
beta=w1.data(:,3);

% we convert energy into wavelength (in nm) using the photon class
PHOTON=photon(PHOTON);

% Note that n = 1 - delta + i*beta = n + i*k
n=1-delta;
k=beta;

% GD-Calc requires the optical constants file to go from shorter to longer
% wavelength, so we flip each array up/down and concatenate the 3 arrays
% into an export array
wnk=[flipud(PHOTON.wavelength) flipud(n) flipud(k)];

% We write the optical constants to file in the format required by GD-Calc
writematrix(wnk,[datadir materialname '_Henke_n_k.txt'],'Delimiter',' ');