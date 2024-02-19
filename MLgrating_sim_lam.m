function [E1,E2,E3,E4,time_elapsed]=MLgrating_sim_lam(grating,PGM,m_max)

% start timer
tStart = tic;

% set incoming photon beam geometry using PGM object
inc_field = MLgrating_def_inc_field(PGM);

% Specify which diffracted orders are to be retained in the calculations,
% Eq's 4.12-13 in GD-Calc.pdf. (m2 is zero because this is a uniperiodic
% grating; all diffraction orders for m2~=0 have zero amplitude.)
order(1).m2 = 0;
order(1).m1 = -m_max:m_max;

% calculate grating fields
[~,scat_field,inc_field] = gdc(grating,inc_field,order,false);

% Compute the diffraction efficiencies.
R = gdc_eff(scat_field,inc_field);

% Discard diffracted waves that decay exponentially with distance from the
% grating. These include evanescent waves.
% (Note: "[scat_field.f1r]" is MATLAB syntax for
% "[scat_field(1).f1r, scat_field(2).f1r, ...]".)
R = R(imag([scat_field.f1r])==0);

% Extract the diffraction order indices for the reflected waves. (Only the
% m1 indices matter; the m2's are all zero.)
R_m1 = [R.m1].';

% Extract the reflection efficiencies (R1...R4) corresponding to the four
% incident polarization states defined in gdc_eff.m.
R1 = [R.eff1].';
R2 = [R.eff2].';
R3 = [R.eff3].';
R4 = [R.eff4].';

% output only efficiency for selected diffraction order
E1 = R1(-PGM.order==R_m1);
E2 = R2(-PGM.order==R_m1);
E3 = R3(-PGM.order==R_m1);
E4 = R4(-PGM.order==R_m1);

% stop timer
time_elapsed=toc(tStart);