function pmt = MLgrating_def_perm(datadir,materials,wavelength)

% find number of materials defined
num_mat=numel(materials);

% initialise variables
w=cell(1,num_mat);
n=cell(1,num_mat);
n_interp=zeros(1,num_mat);
perm=cell(1,num_mat);

for i=1:num_mat
    % read in optical constants (n,k) from datafile using GD-Calc intrinsic 
    % function read_nk.
    [w{i},n{i}] = read_nk([materials{i} '_Henke_n_k.txt'],1,datadir);
    % performs linear interpolation of refractive index for given wavelength
    n_interp(i) = interp1(w{i},n{i},wavelength);
    % Calculate material permittivities from refractive index
    perm{i}=n_interp(i)^2;
end

if isequal(num_mat,3)
    pmt = [{1} perm];  % grating material permittivities 
                        % (including 1 for permittivity in vacuum)
elseif isequal(num_mat,2)
    pmt = [{1} perm {1}];  % grating material permittivities 
                        % (including 1 for permittivity in vacuum)
elseif isequal(num_mat,1)
    pmt = [{1} perm {1} {1}];% grating material permittivities 
                        % (including 1 for permittivity in vacuum)
end
    