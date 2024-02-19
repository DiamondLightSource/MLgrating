function [E1,E2,E3,E4,time_elapsed]=MLgrating_sim(grating,PGM,m_max)

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
% incident polarization states defined below:
%
%     R(k_inc,k).eff1 (parameter, real): efficiency for incident field
%     amplitude = s (i.e., A=1, B=0; linear polarization, TE; sigma)
%
%     R(k_inc,k).eff2 (parameter, real): efficiency for incident field
%     amplitude = p (i.e., A=0, B=1; linear polarization, TM; pi)
%
%     R(k_inc,k).eff3 (parameter, real): efficiency for incident field
%     amplitude = (s+p)/sqrt(2) (i.e., A=B=1/sqrt(2); linear polarization,
%     diagonal)
%
%     R(k_inc,k).eff4 (parameter, real): efficiency for incident field
%     amplitude = (s-i*p)/sqrt(2) (i.e., A=1/sqrt(2), B=-i/sqrt(2);
%     right-hand circular polarization)
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