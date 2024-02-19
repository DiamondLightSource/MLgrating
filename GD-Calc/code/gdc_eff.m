function [R,T]=gdc_eff(scat_field,inc_field)
% gdc_eff
%
% Compute grating diffraction efficiencies from gdc output. (See gdc.m
% comment header.) This version of gdc_eff.m is adapted to work with the
% all_inc_order option (5th input argument for gdc).
%
% Syntax:
%
%   [R,T]=gdc_eff(scat_field,inc_field);
%
% Note: In the following output specification, s and p are the incident
% field's polarization basis vectors (GD-Calc.pdf, equations 4.16-19). For
% an incident E field of arbitrary complex amplitude A*s+B*p, the
% diffraction efficiency in a particular diffraction order (reflected or
% transmitted) is
%   (abs(A)^2*eff1 ...
%   +abs(B)^2*eff2 ...
%   +real(conj(A)*B)*(2*eff3-eff1-eff2) ...
%   -imag(conj(A)*B)*(2*eff4-eff1-eff2)) ...
%   /(abs(A)^2+abs(B)^2)
% (eff1 ... eff4 correspond to the data fields of R or T for a particular
% diffraction order.)
%
% inputs:
%
%   scat_field, inc_field: from gdc
%
% outputs:
%
%   R (size-[?,?] struct): data for reflection efficiencies (row vector if
%   all_inc_order is not used).
%
%     R(k_inc,k).m1, m2 (integer): diffraction order indices (k_inc = 1 if
%     all_inc_order is not used)
%
%     R(k_inc,k).eff1 (parameter, real): efficiency for incident field
%     amplitude = s (i.e., A=1, B=0; linear polarization, TE)
%
%     R(k_inc,k).eff2 (parameter, real): efficiency for incident field
%     amplitude = p (i.e., A=0, B=1; linear polarization, TM)
%
%     R(k_inc,k).eff3 (parameter, real): efficiency for incident field
%     amplitude = (s+p)/sqrt(2) (i.e., A=B=1/sqrt(2); linear polarization,
%     diagonal)
%
%     R(k_inc,k).eff4 (parameter, real): efficiency for incident field
%     amplitude = (s-i*p)/sqrt(2) (i.e., A=1/sqrt(2), B=-i/sqrt(2);
%     right-hand circular polarization)
%
%   T (size-[?,?] struct): data for transmission efficiencies
%
%     T.m1, m2, eff1, eff2, eff3, eff4: same format as R
%
% Notes:
%
% (1) The scat_field output from gdc can be trimmed, before calling
% gdc_eff, to eliminate orders that are not of interest (e.g. evanescent
% orders). With the all_inc_order option, inc_field may also (optionally)
% be trimmed to eliminate uninteresting incident orders, provided that the
% corresponding rows in scat_field are also eliminated.
%
% (2) The diffraction efficiency equation for an evanescent incident order
% may include a divide-by-zero. Under this condition a divide-by-zero
% warning will be generated and the corresponding R and T rows will be nan.
%
% (3) The diffraction efficiency equation is formulated to work with
% non-real-value substrate and/or superstrate permittivity. But if an
% order-m incident wave is non-zero then the order-m reflection efficiency
% will not be physically meaningful unless the superstrate permittivity is
% real-valued and the m-th order is non-evanescent in the superstrate. (If
% these conditions do not hold the order-m power in the superstrate is not
% additively separable into incident and reflected terms; it contains
% multiplicative cross-terms that are neglected by gdc_eff.)
%
% (4) The 02/04/2017 update of gdc_eff makes a correction that applies when
% the substrate and/or superstrate permittivity is not real-valued. Eq's
% 4.35-39 in GD-Calc.pdf contain an obliquity factor, real(f1), that must
% be modified when the permittivity (pmt) is not real-valued. The factor
% real(f1) is correct when applied to s polarization (the |Es|^2 term). But
% for p polarization (the |Ep|^2 term) the correct factor is
% real(f1*abs(pmt)/pmt). The following test case illustrates the error in
% the original gdc_eff and the corrected performance with the update:
%   clear grating inc_field order
%   grating.pmt = {1,1+0.5i};
%   grating.pmt_sub_index = 2; grating.pmt_sup_index = 1;
%   grating.d21 = 1; grating.d31 = 0; grating.d22 = 0; grating.d32 = 1;
%   grating.stratum = {};
%   inc_field.wavelength = 1;
%   inc_field.f2 = 0; inc_field.f3 = 0.5;
%   order.m2 = 0; order.m1 = 0;
%   [~,scat_field,inc_field] = gdc(grating,inc_field,order);
%   [R,T]=gdc_eff(scat_field,inc_field)
%   disp('energy balance (should be 1''s):')
%   [R.eff1+T.eff1;R.eff2+T.eff2;R.eff3+T.eff3;R.eff4+T.eff4]
% (If the superstrate permittivity is not real, the above test will show an
% erroneous energy imbalance because of the neglected power cross-terms;
% see Note (3) above.)

% Version 07/13/2019
% Author: Kenneth C. Johnson, KJ Innovation
% https://doi.org/10.24433/CO.7479617.v3
%
% Acknowledgment of Federal support:
%   This material is based in part upon work supported by Sandia National
%   Laboratories under Purchase Order 1765009.

all_inc_order=isfield(inc_field,'m1_inc');

if ~all_inc_order
    R=struct(...
        'm1',{scat_field.m1},...
        'm2',{scat_field.m2});
else
    R=reshape(struct(...
        'm1_inc',{scat_field.m1_inc},...
        'm2_inc',{scat_field.m2_inc},...
        'm1',{scat_field.m1},...
        'm2',{scat_field.m2}...
        ),size(scat_field));
end
T=R;
% scat_field may be a square matrix if all_inc_order. Use index k for
% flat-array indexing into scat_field, R, and T.
warn_div0 = false;
for k=1:numel(scat_field)
    if ~all_inc_order
        k_inc=1;
    else
        k_inc=find([inc_field.m1_inc]==scat_field(k).m1_inc & ...
            [inc_field.m2_inc]==scat_field(k).m2_inc);
        if isempty(k_inc)
            error(['inc_field is missing order (' ...
                num2str(scat_field(k).m1_inc) ',' ...
                num2str(scat_field(k).m2_inc) ')']);
        end
    end
    [f1i,f2i,f3i,f1r,f1t,f2,f3,Rss,Rsp,Rps,Rpp,Tss,Tsp,Tps,Tpp]=...
        repmat_extend(...
        inc_field(k_inc).f1,...
        inc_field(k_inc).f2,...
        inc_field(k_inc).f3,...
        scat_field(k).f1r,...
        scat_field(k).f1t,...
        scat_field(k).f2,...
        scat_field(k).f3,...
        scat_field(k).Rss,...
        scat_field(k).Rsp,...
        scat_field(k).Rps,...
        scat_field(k).Rpp,...
        scat_field(k).Tss,...
        scat_field(k).Tsp,...
        scat_field(k).Tps,...
        scat_field(k).Tpp);
    % Get obliquity factors OF = real(f1) for s polarization,
    % real(f1*abs(pmt)/pmt) for p polarization.
    OFis = real(-f1i);
    OFip = f1i.^2+f2i.^2+f3i.^2; % pmt/wavelength^2
    OFip = real(-f1i.*abs(OFip)./OFip);
    OFrs = real(f1r);
    OFrp = f1r.^2+f2.^2+f3.^2;
    OFrp = real(f1r.*abs(OFrp)./OFrp);
    OFts = real(-f1t);
    OFtp = f1t.^2+f2.^2+f3.^2;
    OFtp = real(-f1t.*abs(OFtp)./OFtp);
    check_div0 = (OFis(:)==0);
    if any(check_div0)
        if ~warn_div0
            warning(['Divide-by-zero in gdc_eff. ' ...
                '(Evanescent incident order?)']);
            warn_div0 = true;
        end
        OFis(check_div0) = nan;
    end
    check_div0 = (OFip(:)==0);
    if any(check_div0)
        if ~warn_div0
            warning(['Divide-by-zero in gdc_eff. ' ...
                '(Evanescent incident order?)']);
            warn_div0 = true;
        end
        OFip(check_div0) = nan;
    end
    %
    % GD-Calc.pdf, equations 4.38, 4.33
    %
    % Esi = 1, Epi = 0 --> Esr = Rss, Epr = Rps:
    R(k).eff1 = (OFrs.*abs(Rss).^2+OFrp.*abs(Rps).^2)./OFis;
    % Esi = 0, Epi = 1 --> Esr = Rsp, Epr = Rpp:
    R(k).eff2 = (OFrs.*abs(Rsp).^2+OFrp.*abs(Rpp).^2)./OFip;
    % Esi = 1/sqrt(2), Epi = 1/sqrt(2) -->
    %   Esr = (Rss+Rsp)/sqrt(2), Epr = (Rps+Rpp)/sqrt(2):
    R(k).eff3 = (OFrs.*abs(Rss+Rsp).^2+OFrp.*abs(Rps+Rpp).^2)./ ...
        (OFis+OFip);
    % Esi = 1/sqrt(2), Epi = -i/sqrt(2) -->
    %   Esr = (Rss-i*Rsp)/sqrt(2), Epr = (Rps-i*+Rpp)/sqrt(2):
    R(k).eff4 = (OFrs.*abs(Rss-1i*Rsp).^2+OFrp.*abs(Rps-1i*Rpp).^2)./ ...
        (OFis+OFip);
    %
    % GD-Calc.pdf, equations 4.39, 4.34
    %
    % Esi = 1, Epi = 0 --> Est = Tss, Ept = Tps:
    T(k).eff1 = (OFts.*abs(Tss).^2+OFtp.*abs(Tps).^2)./OFis;
    % Esi = 0, Epi = 1 --> Est = Tsp, Ept = Tpp:
    T(k).eff2 = (OFts.*abs(Tsp).^2+OFtp.*abs(Tpp).^2)./OFip;
    % Esi = 1/sqrt(2), Epi = 1/sqrt(2) -->
    %   Est = (Tss+Tsp)/sqrt(2), Ept = (Tps+Tpp)/sqrt(2):
    T(k).eff3 = (OFts.*abs(Tss+Tsp).^2+OFtp.*abs(Tps+Tpp).^2)./ ...
        (OFis+OFip);
    % Esi = 1/sqrt(2), Epi = -i/sqrt(2) -->
    %   Est = (Tss-i*Tsp)/sqrt(2), Ept = (Tps-i*+Tpp)/sqrt(2):
    T(k).eff4 = (OFts.*abs(Tss-1i*Tsp).^2+OFtp.*abs(Tps-1i*Tpp).^2)./ ...
        (OFis+OFip);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [varargout]=repmat_extend(varargin)
% [A1,A2,...]=repmat_extend(A1,A2,...);
% Repmat-extend arguments' singleton dimensions to match array sizes.
% (Scalars will not be repmat-extended.) Size compatibility of input
% arguments is assumed.
n=nargin; % >=nargout
s=cell(1,n);
s_=[1,1];
for j=1:n
    s{j}=size(varargin{j});
    len=length(s{j});
    s_(end+1:len)=1;
    s_(1:len)=max(s_(1:len),s{j});
end
ndims_=length(s_);
n=nargout;
varargout=cell(1,n);
for j=1:n
    if all(s{j}==1)
        varargout{j}=varargin{j}; % (Do not repmat-extend scalars.)
    else
        s{j}(end+1:ndims_)=1;
        ext=s_./s{j};
        if all(ext==1)
            varargout{j}=varargin{j};
        else
            varargout{j}=repmat(varargin{j},ext);
        end
    end
end
