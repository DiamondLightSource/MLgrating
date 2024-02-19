function [R,T] = gdc_eff(scat_field,inc_field)
% gdc_eff
%
% Compute grating diffraction efficiencies from gdc output. (See gdc.m
% comment header.)
%
% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc)
%   Coupled-Wave Theory for Biperiodic DiffractionGratings
%   (GD-Calc.pdf, version 04-Jun-2022)
% All equation ("Eq") references below are from GD-Calc.pdf.
%
% In the following output specification, s and p are the incident field's
% polarization basis vectors (Eq's 4.17-20). For an incident E field of
% arbitrary complex vector amplitude A*s+B*p, the diffraction efficiency
% (eff) in a particular diffraction order (reflected or transmitted) is
% given by Eq's A.28 and A.30:
%   eff =
%   (eff1.*abs(A).^2 ...
%   +eff2.*abs(B).^2 ...
%   +(2*eff3-eff1-eff2).*real(A.*conj(B)) ...
%   +(2*eff4-eff1-eff2).*imag(A.*conj(B)) ...
%   ./(abs(A).^2+abs(B).^2)
% (eff1 ... eff4 correspond to the data fields of R or T for a particular
% diffraction order.)
%
% Syntax:
%
%   [R,T] = gdc_eff(scat_field,inc_field);
%
% Inputs:
%
%   scat_field, inc_field: gdc outputs, or selected subsets of the gdc
%   outputs. Each inc_field(k_inc).f1 must be real-valued and all-nonzero.
%   (The superstrate medium is lossless and the incident wave is
%   non-evanescent, Eq 4.9. See Notes below.)
%
% Outputs:
%
%   R (struct matrix): data for reflection efficiencies, with rows
%   corresponding to incident orders and column corresponding to
%   diffraction orders.
%
%     R(k_inc,k).m1_inc, m2_inc (integer): incident order indices
%
%     R(k_inc,k).m1, m2 (integer): diffraction order indices
%
%     R(k_inc,k).eff1 (parameter, real): efficiency for incident field
%     amplitude = s (i.e., A = 1, B = 0; linear polarization, TE)
%
%     R(k_inc,k).eff2 (parameter, real): efficiency for incident field
%     amplitude = p (i.e., A = 0, B = 1; linear polarization, TM)
%
%     R(k_inc,k).eff3 (parameter, real): efficiency for incident field
%     amplitude = (s+p)/sqrt(2) (i.e., A = B = 1/sqrt(2); linear
%     polarization, diagonal)
%
%     R(k_inc,k).eff4 (parameter, real): efficiency for incident field
%     amplitude = (s-i*p)/sqrt(2) (i.e., A = 1/sqrt(2), B = -i/sqrt(2);
%     left-hand circular polarization -- see comment after Eq A.32)
%
%   T (struct matrix): data for transmission efficiencies
%
%     T.m1_inc, m2_inc, m1, m2, eff1, eff2, eff3, eff4: same format as R
%
% Notes: If the inc_order option is used in gdc, inc_field can be trimmed
% before calling gdc_eff to eliminate evanescent or uninteresting incident
% orders, provided that the corresponding rows in scat_field are also
% eliminated. For example, assuming the inc_field(k_inc).f1 data are
% scalar,
%   f1Inc = [inc_field.f1];
%   k_inc = find(imag(f1_inc)==0 & f1_inc~=0);
%   [R,T] = gdc_eff(scat_field(k_inc,:),inc_field(k_inc));
% After calling gdc_eff, R and T can be trimmed to eliminate uninteresting
% diffracted orders. For example, evanescent diffracted orders can be
% eliminated as follows (assuming that the scat_field(k_inc,k).f1r and
% scat_field(k_inc,k).f1t data are scalar):
%   R = R(imag([scat_field.f1r])==0);
%   T = T(imag([scat_field.f1t])==0);

% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

for k_inc = 1:length(inc_field)
    f1 = inc_field(k_inc).f1;
    if ~isreal(f1)
        error('gdc_eff:validation', ...
            'inc_field(%d).f1 is not real-valued.',k_inc);
    end
    if any(f1(:)==0)
        error('gdc_eff:validation', ...
            'inc_field(%d).f1 is not all-nonzero.',k_inc);
    end
end

R = reshape(struct(...
    'm1_inc',{scat_field.m1_inc},...
    'm2_inc',{scat_field.m2_inc},...
    'm1',{scat_field.m1},...
    'm2',{scat_field.m2}...
    ),size(scat_field));
T = R;
% scat_field is a matrix. Use index k for flat-array indexing into
% scat_field, R, and T.
for k = 1:numel(scat_field)
    k_inc = find([inc_field.m1_inc]==scat_field(k).m1_inc & ...
        [inc_field.m2_inc]==scat_field(k).m2_inc);
    if isempty(k_inc)
        error('gdc_eff:validation', ...
            'inc_field is missing order (%d,%d)', ...
            scat_field(k).m1_inc,scat_field(k).m2_inc);
    end
    [f1i,f1r,f1t,Rss,Rsp,Rps,Rpp,Tss,Tsp,Tps,Tpp] = ...
        deal(...
        inc_field(k_inc).f1,...
        scat_field(k).f1r,...
        scat_field(k).f1t,...
        scat_field(k).Rss,...
        scat_field(k).Rsp,...
        scat_field(k).Rps,...
        scat_field(k).Rpp,...
        scat_field(k).Tss,...
        scat_field(k).Tsp,...
        scat_field(k).Tps,...
        scat_field(k).Tpp);
    % Eq's A.27
    OF = -real(f1r)./f1i; % obliquity factor
    R(k).eff1 = OF.*(abs(Rss).^2+abs(Rps).^2);
    R(k).eff2 = OF.*(abs(Rsp).^2+abs(Rpp).^2);
    tmp1 = 0.5*(R(k).eff1+R(k).eff2);
    tmp2 = Rss.*conj(Rsp)+Rps.*conj(Rpp);
    R(k).eff3 = tmp1+OF.*real(tmp2);
    R(k).eff4 = tmp1-OF.*imag(tmp2);
    % Eq's A.29
    OF = real(f1t)./f1i;
    T(k).eff1 = OF.*(abs(Tss).^2+abs(Tps).^2);
    T(k).eff2 = OF.*(abs(Tsp).^2+abs(Tpp).^2);
    tmp1 = 0.5*(T(k).eff1+T(k).eff2);
    tmp2 = Tss.*conj(Tsp)+Tps.*conj(Tpp);
    T(k).eff3 = tmp1+OF.*real(tmp2);
    T(k).eff4 = tmp1-OF.*imag(tmp2);
end % for k = 1:numel(scat_field)
