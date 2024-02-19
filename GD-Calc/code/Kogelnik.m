function [effTE,effTM,transmission] = ...
    Kogelnik(d,vpd,pmt,dpmt,lambda,theta,phi)
% Calculate the diffraction efficiency of a dielectric holographic volume
% grating (transmission or reflection) via Kogelnik theory (1st-order
% Bragg diffraction, external surface reflections are neglected):
%   Kogelnik, H. 1969, The Bell System Technical Journal, 48, 2909
%   https://doi.org/10.1002/j.1538-7305.1969.tb01198.x
%
% Syntax:
%   [effTE,effTM] = Kogelnik(d,vpd,pmt,dpmt,lambda,theta,phi);
%
% Inputs: (All inputs can be arrays, size-compatible with implicit
% singleton expansion.)
%
%   d: grating thickness
%
%   vpd: volumetric grating period, signed for selecting between +1st
%   diffraction order (vpd>0) or -1st order (vpd<0)
%
%   pmt: average permittivity, real (Refractive index is n = sqrt(pmt).)
%
%   dpmt: permittivity modulation, peak-to-mean, real (Refractive index
%   modulation is dn = dpmt/(2*n).)
%
%   lambda: wavelength in vacuum
%
%   theta: incidence angle inside grating medium, deg (In relation to the
%   vaccum angle, thetaVac, sind(theta) = sind(thetaVac)/sqrt(pmt).)
%
%   phi: grating slant angle relative to grating surface, deg (equal to
%   90deg plus slant relative to lamellar grating)
%
% Outputs:
%
%   effTE, effTM: diffraction efficiencies in +/-1st diffraction order for
%   TE and TM efficiency
%
%   transmission: true for transmission grating, false for reflection
%   grating
%
% An error is generated if the grating is not unambiguously a transmission
% or reflection grating.

% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% "Eq" references are from Kogelnik.

if ~(isreal(pmt) && isreal(dpmt))
    error('Kogelnik:validation','pmt and dpmt must be real-valued.')
end

n = sqrt(pmt);
dn = dpmt./(2*n);
% Eq 3
K = 2*pi./vpd;
% Eq 5
beta = 2*pi*n./lambda;
% Eq 13
sin_theta = sind(theta);
cos_theta = cosd(theta);
sigma_x = beta.*sin_theta-K.*sind(phi);
sigma_z = beta.*cos_theta-K.*cosd(phi);
% Eq 17
theta_ = (beta.^2-sigma_x.^2-sigma_z.^2)./(2*beta);
% Eq 23
cR = cos_theta;
cS = sigma_z./beta; % positive for transmission, negative for reflection
% Eqs 42, 56
nu = pi*dn.*d./(lambda.*sqrt(cR.*abs(cS)));
xi = theta_.*d./(2*abs(cS));
% Eq 90
nu_ = -nu.*cosd(2*(theta-phi));
if all(cS(:)>0)
    transmission = true;
    % Eq 43
    effTE = sin(sqrt(nu.^2+xi.^2)).^2./(1+(xi./nu).^2);
    effTM = sin(sqrt(nu_.^2+xi.^2)).^2./(1+(xi./nu_).^2);
elseif all(cS(:)<0)
    transmission = false;
    % Eq 57
    effTE = 1./(1+(1-(xi./nu).^2)./sinh(sqrt(nu.^2-xi.^2)).^2);
    effTM = 1./(1+(1-(xi./nu_).^2)./sinh(sqrt(nu_.^2-xi.^2)).^2);
else
    error('Kogelnik:ambiguousGrating', ...
        ['The grating is not unambiguously a transmission or ' ...
        'reflection grating.'])
end
