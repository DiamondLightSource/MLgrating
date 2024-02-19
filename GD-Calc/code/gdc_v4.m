function [param_size,scat_field,inc_field,grating] = ...
    gdc(grating,inc_field,order,show_progress,inc_order,PadeOrder,RelTol)
%
% gdc (Grating Diffraction Calculator)
%
% Compute the scattered electromagnetic field from a biperiodic grating.
%
% Syntax:
%
%   Compute the scattered electromagnetic field:
%   [param_size,scat_field,inc_field] = gdc(grating,inc_field,order);
%   ...
%   [param_size,scat_field,inc_field,grating] = ...
%       gdc(grating,inc_field,order,...
%       show_progress,inc_order,PadeOrder,RelTol);
%
%   Just validate input(s); do not run computations:
%   gdc(grating,...);
%   param_size = gdc(grating,...);
%
%   Just display version and copyright information:
%   gdc;
%
% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc)
%   Coupled-Wave Theory for Biperiodic DiffractionGratings
%   (GD-Calc.pdf, version 04-Jun-2022)
% All equation ("Eq") references below are from GD-Calc.pdf.
%
% Note on parameterization:
%   In the following interface specification "parameters" are
%   multi-dimensional numeric arrays, which must be non-empty and must all
%   be size-matched except for singleton dimensions. Singleton dimensions
%   are implicitly repmat-expanded, if necessary, to match parameter
%   sizes. See the following MATLAB documentation pages:
%   https://www.mathworks.com/help/matlab/matlab_prog/array-vs-matrix-operations.html
%   https://www.mathworks.com/help/matlab/matlab_prog/compatible-array-sizes-for-basic-operations.html
%
% Note on coordinate breaks:
%   A coordinate break is represented as a type of "stratum" (of zero
%   thickness), in the sense that it is associated with a lateral plane at
%   a particular x1 height in the grating.
%
% Note on replication modules:
%   A "replication module" is a repeating pattern of grating strata, which
%   are collectively treated as a single composite stratum. The module
%   comprises a list of strata, which are stacked from bottom to top in the
%   specified order to form the module, and a replication count
%   (rep_count), which specifies how many copies of the module are to be
%   stacked. The stacking operation's runtime is proportional to the
%   logarithm of rep_count, so a very large number of thin modules may
%   typically be stacked without incurring much computational overhead.
%   (For optimum efficiency, rep_count should preferably be a power of 2.)
%   A replication module's component strata can be of any type
%   (homogeneous, uniperiodic, biperiodic, coordinate break, or nested
%   replication module).
%
% Note on the progress bar:
%   A progress bar will be displayed and will update once per stratum
%   unless it is disabled by setting show_progress = false. (The progress
%   bar can sometimes incur significant runtime overhead.)
%
% Inputs:
%
%   grating (struct): grating specification
%
%     grating.pmt (cell row vector): complex permittivities for grating
%     materials
%
%       grating.pmt{k} (parameter, complex): permittivity for k-th material
%       (must be non-zero; imaginary part must be non-negative.)
%
%     grating.pmt_sub_index, pmt_sup_index (integer): indices into
%     grating.pmt for substrate and superstrate permittivities (Eq's 3.11
%     and 3.12)
%
%     grating.d21, d31, d22, d32 (parameter, real): grating period vectors
%     (Eq's 3.5 and 3.6) The 1st and 2nd period vectors' [x2,x3]
%     coordinates are [d21,d31] and [d22,d32], respectively. The vectors
%     must be linearly independent (Eq 3.7).
%
%     grating.stratum (size-[1,L1] cell array or {}): stratum
%     specifications (L1 = number of strata, may be zero)
%
%       grating.stratum{l1} (struct): specification for stratum l1
%
%         grating.stratum{l1}.type (0, 1, 2, 3, or 4): type of stratum - 0
%         for homogeneous, 1 for uniperiodic, 2 for biperiodic, 3 for
%         coordinate break, 4 for replication module.
%
%         for a homogeneous stratum (stratum type is 0):
%
%           grating.stratum{l1}.thick (parameter, non-negative real):
%           stratum thickness
%
%           grating.stratum{l1}.pmt_index (integer): index into grating.pmt
%           for stratum permittivity
%
%         for a uniperiodic stratum (stratum type is 1):
%
%           grating.stratum{l1}.thick (parameter, non-negative real):
%           stratum thickness
%
%           grating.stratum{l1}.h11, h12 (integer): harmonic indices (Eq's
%           3.24-25), not both zero (Eq 3.23). Use h11 = 1, h12 = 0 if the
%           stratum period matches the first grating period.
%
%           grating.stratum{l1}.stripe (size-[1,L2] cell array): stripe
%           specifications (L2 = number of stripes in stratum l1, must be
%           nonzero)
%
%             grating.stratum{l1}.stripe{l2} (struct): specification for
%             stripe l2 in stratum l1
%
%               grating.stratum{l1}.stripe{l2}.c1 (parameter, real): stripe
%               boundary position, must satisfy monoticity condition (Eq
%               3.35). The stripe's positive-side boundary is at c1 in
%               period-normalized units.
%
%               grating.stratum{l1}.stripe{l2}.pmt_index (integer): index
%               into grating.pmt for stripe permittivity (on negative side
%               of c1 boundary)
%
%         for a biperiodic stratum (stratum type is 2):
%
%           grating.stratum{l1}.thick (parameter, non-negative real):
%           stratum thickness
%
%           grating.stratum{l1}.h11, h12, h21, h22 (integer): harmonic
%           indices (Eq 3.20); must satisfy Eq 3.19. Use h11 = 1, h12 = 0,
%           h21 = 0, h22 = 1 if the stratum periods match the grating
%           periods.
%
%           grating.stratum{l1}.stripe (length-L2 cell array): stripe
%           specifications (L2 = number of stripes in stratum l1, must be
%           nonzero)
%
%             grating.stratum{l1}.stripe{l2} (struct): specification for
%             stripe l2 in stratum l1
%
%               grating.stratum{l1}.stripe{l2}.c1 (parameter, real): stripe
%               boundary position; must satisfy monoticity condition (Eq
%               3.35) The stripe's positive-side boundary is at c1 in
%               period-normalized units.
%
%               grating.stratum{l1}.stripe{l2}.type (0 or 1): type of
%               stripe - 0 for homogeneous, 1 for inhomogeneous
%
%               for a homogeneous stripe (stripe type is 0):
%
%                 grating.stratum{l1}.stripe{l2}.pmt_index (integer): index
%                 into grating.pmt for stripe permittivity (on negative
%                 side of c1 boundary)
%
%               for an inhomogeneous stripe (stripe type is 1):
%
%                 grating.stratum{l1}.stripe{l2}.block (size-[1,L3] cell
%                 array): block specifications (L3 = number of blocks in
%                 stripe l2 of stratum l1, must be nonzero)
%
%                   grating.stratum{l1}.stripe{l2}.block{l3}: specification
%                   for block l3 in stripe l2 of stratum l1
%
%                     grating.stratum{l1}.stripe{l2}.block{l3}.c2
%                     (parameter, real): block boundary position, must
%                     satisfy monoticity condition (Eq 3.37) The block's
%                     positive-side boundary is at c2 in period-normalized
%                     units.
%
%                     grating.stratum{l1}.stripe{l2}.block{l3}.pmt_index
%                     (integer): index into grating.pmt for block
%                     permittivity (on negative side of c2 boundary)
%
%         for a coordinate break (stratum type is 3):
%
%           grating.stratum{l1}.dx2, dx3 (parameter, real): lateral
%           translational offset (Eq 3.41) The coordinate break effects a
%           translational displacement of all higher strata by [dx2,dx3].
%
%         for a replication module (stratum type is 4):
%
%           grating.stratum{l1}.stratum (cell row vector, may be empty):
%           stratum specifications for module's strata (same data format as
%           grating.stratum)
%
%           grating.stratum{l1}.rep_count (non-negative integer): module
%           replication count
%
%   inc_field (struct): incident field specification
%
%     inc_field.wavelength (parameter, positive real): wavelength (same
%     length units as grating.d21, etc.)
%
%     inc_field.f2, f3 (parameter, real): incident field's zero-order
%     tangential spatial frequency vector (reciprocal-length units, Eq's
%     4.1 and 4.2)
%
%   order (nonempty struct row vector): diffraction orders to be retained
%   in calculations, non-empty
%
%     order(k).m2 (integer): m2 index from M2 (Eq 4.12). The m2 indices
%     must be unique, but need not be sorted. One of the specified m2
%     indices must be zero.
%
%     order(k).m1 (nonempty integer row vector): m1 indices corresponding
%     to m2 (i.e., from M1, Eq 4.13), non-empty. For each m2, the
%     associated m1 indices must be unique, but need not be sorted. For
%     m2==0, one of the associated m1 indices must be zero.
%
%   show_progress: true or false, OPTIONAL, default is true. If true, a
%   progress bar (waitbar) will be displayed while gdc is executing.
%
%   inc_order: true or false, or 2-column integer matrix, OPTIONAL, default
%   is [0,0]. Set to inc_order = true to enable multi-wave incident field,
%   or set to a list of [m1,m2] order index pairs to select specific
%   incident orders. (See "Advanced features ... inc_order option" below.)
%
%   PadeOrder: positive even integer, OPTIONAL, default is 6.
%   Approximation order for Padé approximation (n in Eq's D.6 and D.16)
%
%   RelTol: positive real, OPTIONAL, default is eps. Relative accuracy
%   target in Eq D.27
%
% Outputs:
%
%   param_size (integer row vector): Common parameter size, with full
%   repmat expansion. (If the inc_field argument in missing, param_size is
%   determined only on the basis of parameters in the grating input.)
%
%   scat_field (struct row vector): scattered field. Each element
%   scat_field(k) corresponds to one diffraction order (transmitted and
%   reflected fields)
%
%     scat_field(k).m1_inc, m2_inc (integer): incident order indices, both
%     zero unless inc_order option is used. (See "Advanced features ...
%     inc_order option" below.)
%
%     scat_field(k).m1, m2 (integer): diffraction order indices
%
%     scat_field(k).f2, f3 (parameter, real): field's grating-tangential
%     spatial frequencies (Eq 4.7)
%
%     scat_field(k).f1r, f1t (parameter, complex): reflected and
%     transmitted field's grating-normal spatial frequencies (Eq's 4.10,
%     4.11)
%
%     scat_field(k).s2, s3 (parameter, real): grating-tangential basis
%     vector for S polarization (Eq's 4.24, 4.29)
%
%     scat_field(k).p1t, p2t, p3t (parameter, complex): basis vector for P
%     polarization, for transmitted field (Eq 4.30)
%
%     scat_field(k).p1r, p2r, p3r (parameter, complex): basis vector for P
%     polarization, for reflected field (Eq 4.25)
%
%     scat_field(k).Tss, Tsp, Tps, Tpp (parameter, complex): transmission
%     matrix (Eq's 4.33, 4.35)
%
%     scat_field(k).Rss, Rsp, Rps, Rpp (parameter, complex): reflection
%     matrix (Eq's 4.32, 4.34)
%
%   inc_field (struct): incident field specification
%
%     inc_field.m1_inc, m2_inc: incident order indices, both zero unless
%     inc_order option is used. (See "Advanced features ... inc_order
%     option" below.)
%
%     inc_field.wavelength, f2, f3: same as inc_field input argument (Eq's
%     4.1, 4.2)
%
%     inc_field.f1 (parameter, complex): grating-normal spatial frequency
%     (Eq 4.9)
%
%     inc_field.s2, s3 (parameter, real): grating-tangential basis vector
%     for S polarization (Eq 4.19)
%
%     inc_field.p1, p2, p3 (parameter, complex): basis vector for P
%     polarization (Eq 4.20)
%
%
% Advanced features
%
%
% inc_order option:
%
%   The optional input argument inc_order is either a logical or a nonempty
%   two-column matrix of [m1,m2] incident order indices (default = [0,0]).
%   inc_order can be set to true to get the scattered field for all
%   incident orders (not just the zero-order incident wave), or inc_order
%   can be set to a list of specific incident orders. inc_order = false is
%   equivalent to [0,0] (the default), which selects the zero order.
%   inc_order = true selects all [m1,m2] index pairs defined by the order
%   struct. The [m1,m2] pairs defined by inc_order must be unique and
%   within the set defined by the order struct.
%
%   With inc_order specified, the following changes apply:
%
%   The inc_field output struct will be a row vector with elements
%   inc_field(k_inc) corresponding to different incident orders. The order
%   indices corresponding to inc_field(k_inc) will be specified as
%   inc_field(k_inc).m1_inc, inc_field(k_inc).m2_inc.
%
%   The scat_field output struct will be a matrix with rows
%   scat_field(k_inc,:) representing the scattered field for
%   inc_field(k_inc). The incident-order indices for scat_field(k_inc,k)
%   will be stored as scat_field(k_inc,k).m1_inc,
%   scat_field(k_inc,k).m2_inc (identical to inc_field(k_inc).m1_inc,
%   inc_field(k_inc).m2_inc).
%
%
% Full-field option:
%
%   The full electromagnetic field (E and H vectors) can be evaluated in
%   the substrate, in the superstrate, or in any stratum of type 0, 1, or 2
%   (provided that the stratum is not in a replication module). To do a
%   full-field calculation, include "full_field" data fields in the grating
%   struct (as outlined below) and include grating as an output argument:
%
%     [param_size,scat_field,inc_field,grating] = ...
%         gdc(grating,inc_field,order,...);
%
%   Make the following initializations in the grating struct (in any
%   combination - all are optional) to calculate the full field in the
%   substrate, in the superstrate, or in a stratum:
%
%     grating.sub_full_field = []; % To calculate substrate full field.
%
%     grating.sup_full_field = []; % To calculate superstrate full field.
%
%     grating.stratum{l1}.full_field = []; % To calculate stratum-l1 full
%     field.
%
%   In the last line above, the stratum type must be 0, 1, or 2 and the
%   stratum cannot be in a replication module. If any of the above data
%   fields exists, but is initialized to a value other than [], it will be
%   reinitialized to []. After running gdc, each of the above data fields
%   will be set to a struct row vector containing the E and H fields'
%   Fourier amplitudes ffE and ffH (Eq's 5.13 and 5.14), evaluated at x2 =
%   0, x3 = 0, and x1 as follows (with b1 defined in Eq's 3.1-3):
%
%     substrate: x1 = b1[0]-  (below grating)
%
%     superstrate: x1 = b1[L1]+  (above grating)
%
%     in stratum l1: x1 = b1[l1]-  (at top stratum interface)
%
%   The full-field result in the superstrate (grating.sup_full_field)
%   represents the total (incident + diffracted) field, not just the
%   diffracted field.
%
%   Note that the internal field is only calculated at each stratum's top
%   interface. If a vertical profile of the internal field within a grating
%   layer is required, the layer should be split into a number of strata.
%   For example, the grating construction code typically increments the
%   stratum list as follows:
%
%     grating.stratum{end+1} = stratum;
%
%   This line can be replaced by the following in order to sample the
%   stratum layer's internal field at 10 equidistant height locations:
%
%     stratum.full_field = [];
%     stratum.thick = stratum.thick/10;
%     for j = 1:10
%       grating.stratum{end+1} = stratum;
%     end
%
%   After running gdc the full_field data fields will be struct vectors
%   with elements corresponding to the E and H fields' Fourier components,
%   formatted as follows:
%
%     full_field(k).m1, m2 (integer): diffraction order indices
%
%     full_field(k).wavelength (parameter, real): wavelength
%
%     full_field(k).f2, f3 (parameter, real): field's grating-tangential
%     spatial frequencies (Eq 5.8)
%
%     full_field(k).ffE1s, ffE2s, ffE3s (parameter, complex): E-field
%     amplitudes for Fourier order (m1, m2) (Eq 5.13), with incident E
%     field s-polarized and of unit amplitude
%
%     full_field(k).ffH1s, ffH2s, ffH3s (parameter, complex): same as
%     above, but for H field
%
%     full_field(k).ffE1p, ffE2p, ffE3p, ffH1p, ffH2p, ffH3p (parameter,
%     complex): same as above, but for p-polarized incident E field
%
%   If a stratum of type 0, 1, or 2 inside a replication module contains a
%   full_field data field, no error will be generated but the full_field
%   result will be empty ([]).
%
%   (Note: The full-field option is not currently implemented to work with
%   inc_order. Both can be used together, but the full-field data is only
%   generated for the zero-order incident field.)
%

% Version 04-Jun-2022
% Copyright (c) 2005-2022, Kenneth C. Johnson. All rights reserved.
% KJ Innovation https://kjinnovation.com/

if nargout>4 || (nargin==0 && nargout>0)
    error('gdc:nargout','Too many output arguments.')
end
if nargin==0
    if ~any(exist('gdc_engine','file')==[2,6])
        error('gdc:UndefinedFunction', ...
            'gdc_engine is missing or cannot be found.');
    end
    gdc_engine;
    return
end
if (nargout>1 && ~any(nargin==3:7)) || ~any(nargin==1:7)
    error('gdc:nargin','Wrong number of input argument(s).');
end
if ~isstruct(grating) || length(grating)~=1
    error('gdc:validation','Wrong data type or size (grating).');
end
if ~isfield(grating,'pmt')
    error('gdc:validation','Missing data field (grating.pmt).');
end
if ~iscell(grating.pmt) || isempty(grating.pmt) || ~isrow(grating.pmt)
    error('gdc:validation','Wrong data type or size (grating.pmt).');
end
param_size = [1 1];
for k = 1:length(grating.pmt)
    [param_size,err_msg] = check_param(param_size,grating.pmt{k});
    if ~isempty(err_msg)
        error('gdc:validation',[err_msg ' (grating.pmt{%d}).'],k);
    end
    if any(grating.pmt{k}(:)==0) || any(imag(grating.pmt{k}(:))<0)
        error('gdc:validation',['Permittivity must be non-zero; '...
            'imaginary part must be non-negative (grating.pmt{%d}).'],k);
    end
end
if ~isfield(grating,'pmt_sub_index')
    error('gdc:validation','Missing data field (grating.pmt_sub_index).');
end
if ~(check_integer(grating.pmt_sub_index) && ...
        isscalar(grating.pmt_sub_index))
    error('gdc:validation', ...
        'Wrong data type or size (grating.pmt_sub_index).');
end
if grating.pmt_sub_index<1 || grating.pmt_sub_index>length(grating.pmt)
    error('gdc:validation','grating.pmt_sub_index is out of range.');
end
if ~isfield(grating,'pmt_sup_index')
    error('gdc:validation','Missing data field (grating.pmt_sup_index).');
end
if ~(check_integer(grating.pmt_sup_index) && ...
        isscalar(grating.pmt_sup_index))
    error('gdc:validation', ...
        'Wrong data type or size (grating.pmt_sup_index).');
end
if grating.pmt_sup_index<1 || grating.pmt_sup_index>length(grating.pmt)
    error('gdc:validation','grating.pmt_sup_index is out of range.');
end
if ~isfield(grating,'d21')
    error('gdc:validation','Missing data field (grating.d21).');
end
if ~check_real(grating.d21)
    error('gdc:validation','grating.d21 must be real.');
end
[param_size,err_msg] = check_param(param_size,grating.d21);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (grating.d21).']);
end
if ~isfield(grating,'d31')
    error('gdc:validation','Missing data field (grating.d31).');
end
if ~check_real(grating.d31)
    error('gdc:validation','grating.d31 must be real.');
end
[param_size,err_msg] = check_param(param_size,grating.d31);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (grating.d31).']);
end
if ~isfield(grating,'d22')
    error('gdc:validation','Missing data field (grating.d22).');
end
if ~check_real(grating.d22)
    error('gdc:validation','grating.d22 must be real.');
end
[param_size,err_msg] = check_param(param_size,grating.d22);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (grating.d22).']);
end
if ~isfield(grating,'d32')
    error('gdc:validation','Missing data field (grating.d32).');
end
if ~check_real(grating.d32)
    error('gdc:validation','grating.d32 must be real.');
end
[param_size,err_msg] = check_param(param_size,grating.d32);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (grating.d32).']);
end
if any(grating.d21.*grating.d32-grating.d31.*grating.d22==0,'all')
    % violation of Eq 3.7
    error('gdc:validation', ...
        'Grating period vectors must be linearly independent.');
end
if isfield(grating,'sub_full_field')
    grating.sub_full_field = [];
end
if isfield(grating,'sup_full_field')
    grating.sup_full_field = [];
end
if ~isfield(grating,'stratum')
    error('gdc:validation','Missing data field (grating.stratum).');
end
if ~iscell(grating.stratum) || ...
        (~isempty(grating.stratum) && ~isrow(grating.stratum))
    error('gdc:validation','Wrong data type or size (grating.stratum).');
end
L1 = length(grating.stratum);
for l1 = 1:L1
    [param_size,grating.stratum{l1}] = ...
        check_stratum(param_size,grating,grating.stratum{l1},...
        ['grating.stratum{' num2str(l1) '}']);
end
if nargin==1
    return
end
if ~(isstruct(inc_field) && isscalar(inc_field))
    error('gdc:validation','Wrong data type or size (inc_field).');
end
if ~isfield(inc_field,'wavelength')
    error('gdc:validation','Missing data field (inc_field.wavelength).');
end
[param_size,err_msg] = check_param(param_size,inc_field.wavelength);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (inc_field.wavelength).']);
end
if ~check_real(inc_field.wavelength) || any(inc_field.wavelength(:)<=0)
    error('gdc:validation', ...
        'inc_field.wavelength must be real and positive.');
end
if ~isfield(inc_field,'f2')
    error('gdc:validation','Missing data field (inc_field.f2).');
end
[param_size,err_msg] = check_param(param_size,inc_field.f2);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (inc_field.f2).']);
end
if ~check_real(inc_field.f2)
    error('gdc:validation','inc_field.f2 must be real.');
end
if ~isfield(inc_field,'f3')
    error('gdc:validation','Missing data field (inc_field.f3).');
end
[param_size,err_msg] = check_param(param_size,inc_field.f3);
if ~isempty(err_msg)
    error('gdc:validation',[err_msg ' (inc_field.f3).']);
end
if ~check_real(inc_field.f3)
    error('gdc:validation','inc_field.f3 must be real.');
end
if nargin==2
    return
end
if ~(isstruct(order) && isrow(order) && ~isempty(order))
    error('gdc:validation','Wrong data type or size (order).');
end
if ~isfield(order,'m2')
    error('gdc:validation','Missing data field (order(...).m2).');
end
if ~isfield(order,'m1')
    error('gdc:validation','Missing data field (order(...).m1).');
end
for k = 1:length(order)
    m = order(k).m2;
    if ~(check_integer(m) && isscalar(m))
        error('gdc:validation', ...
            'Wrong data type or size (order(%d).m2).',k);
    end
    m = order(k).m1;
    if ~(check_integer(m) && isrow(m) && ~isempty(m))
        error('gdc:validation', ...
            'Wrong data type or size (order(%d).m1).',k);
    end
    if any(diff(sort(m))==0)
        error('gdc:validation','order(%d).m1 elements must be unique.',k);
    end
end
m = [order.m2];
if any(diff(sort(m))==0)
    error('gdc:validation','[order.m2] elements must be unique.');
end
k = find(m==0);
if isempty(k) || ~any(order(k).m1==0)
    error('gdc:validation', ...
        'Missing zero order in specified order indices.');
end
if nargout<=1
    return
end
clear m
if ~any(exist('gdc_engine','file')==[2,6])
    error('gdc:UndefinedFunction', ...
        'gdc_engine is missing or cannot be found.');
end
if nargin<4
    show_progress = true;
elseif ~(islogical(show_progress) && isscalar(show_progress))
    error('gdc:validation','show_progress must be true or false.');
end
if nargin<5 || isequal(inc_order,false)
    inc_order = [0,0];
else
    m = zeros(0,2);
    for k = 1:length(order)
        m_ = order(k).m1.';
        m_(:,2) = order(k).m2;
        m = [m;m_]; %#ok<AGROW>
    end
    if islogical(inc_order)
        if ~isscalar(inc_order)
            error('gdc:validation',['inc_order must be true or false ' ...
                'or a nonempty 2-column integer matrix.']);
        end
        % inc_order = true
        inc_order = m;
    else
        if ~(check_integer(inc_order) && ismatrix(inc_order) && ...
                size(inc_order,2)==2 && ~isempty(inc_order))
            error('gdc:validation',['inc_order must be true or false ' ...
                'or a nonempty 2-column integer matrix.']);
        end
        if any(all(diff(sortrows(inc_order),1,1)==0,2))
            error('gdc:validation','inc_order rows must be unique');
        end
        if ~isempty(setdiff(inc_order,m,'rows'))
            error('gdc:validation',['inc_order rows must be ' ...
                'in the set defined by the order struct.']);
        end
    end
end
% inc_order is a 2-column matrix of [m1,m2] order indices.
if nargin<6
    PadeOrder = 6;
elseif ~(check_integer(PadeOrder) && isscalar(PadeOrder) && ...
        PadeOrder>0 && mod(PadeOrder,2)==0)
    error('gdc:validation','PadeOrder must be a positive even integer.')
end
if nargin<7
    RelTol = eps;
elseif ~(check_real(RelTol) && isscalar(RelTol) && RelTol>0)
    error('gdc:validation','RelTol must be a positive real scalar.')
end
if nargout<4
    [param_size,scat_field,inc_field] = gdc_engine(...
        param_size,grating,inc_field,order, ...
        show_progress,inc_order,PadeOrder,RelTol);
else
    [param_size,scat_field,inc_field,grating] = gdc_engine(...
        param_size,grating,inc_field,order, ...
        show_progress,inc_order,PadeOrder,RelTol);
end

end % gdc


function [param_size,stratum] = ...
    check_stratum(param_size,grating,stratum,stratum_name)
if ~(isstruct(stratum) && isscalar(stratum))
    error('gdc:validation', ...
        ['Wrong data type or size (' stratum_name ').']);
end
if ~isfield(stratum,'type')
    error('gdc:validation', ...
        ['Missing data field (' stratum_name '.type).']);
end
if ~(check_integer(stratum.type) && isscalar(stratum.type))
    error('gdc:validation', ...
        ['Wrong data type or size (' stratum_name '.type).']);
end
if ~any(stratum.type==[0 1 2 3 4])
    error('gdc:validation', ...
        ['Stratum type must be 0, 1, 2, 3, or 4 (' ...
        stratum_name '.type).']);
end
if any(stratum.type==[0,1,2])
    if ~isfield(stratum,'thick')
        error('gdc:validation', ...
            ['Missing data field (' stratum_name '.thick).']);
    end
    [param_size,err_msg] = check_param(param_size,stratum.thick);
    if ~isempty(err_msg)
        error('gdc:validation', ...
            [err_msg ' (' stratum_name '.thick).']);
    end
    if ~check_real(stratum.thick) || any(stratum.thick(:)<0)
        error('gdc:validation', ...
            ['Stratum thickness must be real and non-negative (' ...
            stratum_name '.thick).']);
    end
    if stratum.type==0
        if ~isfield(stratum,'pmt_index')
            error('gdc:validation', ...
                ['Missing data field (' stratum_name '.pmt_index).']);
        end
        if ~(check_integer(stratum.pmt_index) && ...
                isscalar(stratum.pmt_index))
            error('gdc:validation', ...
                ['Wrong data type or size (' stratum_name '.pmt_index).']);
        end
        if stratum.pmt_index<1 || stratum.pmt_index>length(grating.pmt)
            error('gdc:validation', ...
                ['pmt_index is out of range (' ...
                stratum_name '.pmt_index).']);
        end
    else % stratum.type==1 or 2
        if ~isfield(stratum,'h11')
            error('gdc:validation', ...
                ['Missing data field (' stratum_name '.h11).']);
        end
        if ~(check_integer(stratum.h11) && isscalar(stratum.h11))
            error('gdc:validation', ...
                ['Wrong data type or size (' stratum_name '.h11).']);
        end
        if ~isfield(stratum,'h12')
            error('gdc:validation', ...
                ['Missing data field (' stratum_name '.h12).']);
        end
        if ~(check_integer(stratum.h12) && isscalar(stratum.h12))
            error('gdc:validation', ...
                ['Wrong data type or size (' stratum_name '.h12).']);
        end
        if stratum.type==1
            if stratum.h11==0 && stratum.h12==0
                % violation of Eq 3.23
                error('gdc:validation', ...
                    ['h11 and h12 must not both be zero (' ...
                    stratum_name '.h11 and h12).']);
            end
        else
            if ~isfield(stratum,'h21')
                error('gdc:validation', ...
                    ['Missing data field (' stratum_name '.h21).']);
            end
            if ~(check_integer(stratum.h21) && isscalar(stratum.h21))
                error('gdc:validation', ...
                    ['Wrong data type or size (' stratum_name '.h21).']);
            end
            if ~isfield(stratum,'h22')
                error('gdc:validation', ...
                    ['Missing data field (' stratum_name '.h22).']);
            end
            if ~(check_integer(stratum.h22) && isscalar(stratum.h22))
                error('gdc:validation', ...
                    ['Wrong data type or size (' stratum_name '.h22).']);
            end
            if stratum.h11*stratum.h22==stratum.h12*stratum.h21
                % violation of Eq 3.19
                error('gdc:validation', ...
                    ['h11, h12, h21, and h22 must satisfy ' ...
                    'h11*h22~=h12*h21 (' stratum_name ...
                    '.h11, h12, h21, and h22).']);
            end
        end
        if ~isfield(stratum,'stripe')
            error('gdc:validation', ...
                ['Missing data field (' stratum_name '.stripe).']);
        end
        if ~(iscell(stratum.stripe) && isrow(stratum.stripe) ...
                && ~isempty(stratum.stripe))
            error('gdc:validation', ...
                ['Wrong data type or size (' stratum_name '.stripe).']);
        end
        L2 = length(stratum.stripe);
        for l2 = 1:L2
            stripe = stratum.stripe{l2};
            if ~isfield(stripe,'c1')
                error('gdc:validation', ...
                    ['Missing data field (' stratum_name ...
                    '.stripe{%d}.c1).'],l2);
            end
            [param_size,err_msg] = check_param(param_size,stripe.c1);
            if ~isempty(err_msg)
                error('gdc:validation', ...
                    [err_msg ' (' stratum_name '.stripe{%d}.c1).'],l2);
            end
            if ~check_real(stripe.c1)
                error('gdc:validation', ...
                    ['c1 must be real (' stratum_name ...
                    '.stripe{%d}.c1).'],l2);
            end
            if stratum.type==2
                if ~isfield(stripe,'type')
                    error('gdc:validation', ...
                        ['Missing data field (' stratum_name ...
                        '.stripe{%d}.type).'],l2);
                end
                if ~(check_integer(stripe.type) && isscalar(stripe.type))
                    error('gdc:validation', ...
                        ['Wrong data type or size (' stratum_name ...
                        '.stripe{%d}.type).'],l2);
                end
                if ~any(stripe.type==[0 1])
                    error('gdc:validation', ...
                        ['Stripe type must be 0 or 1 (' stratum_name ...
                        '.stripe{%d}.type).'],l2);
                end
            else
                if isfield(stripe,'type')
                    error('gdc:validation', ...
                        ['Extraneous data field (' stratum_name ...
                        '.stripe{%d}.type).'],l2);
                end
            end
            if stratum.type==1 || stripe.type==0
                if ~isfield(stripe,'pmt_index')
                    error('gdc:validation', ...
                        ['Missing data field (' stratum_name ...
                        '.stripe{%d}.pmt_index).'],l2);
                end
                if ~(check_integer(stripe.pmt_index) && ...
                        isscalar(stripe.pmt_index))
                    error('gdc:validation', ...
                        ['Wrong data type or size (' stratum_name ...
                        '.stripe{%d}.pmt_index).'],l2);
                end
                if stripe.pmt_index<1 || ...
                        stripe.pmt_index>length(grating.pmt)
                    error('gdc:validation', ...
                        ['pmt_index is out of range (' stratum_name ...
                        '.stripe{%d}.pmt_index).'],l2);
                end
                if isfield(stripe,'block')
                    error('gdc:validation', ...
                        ['Extraneous data field (' stratum_name ...
                        '.stripe{%d}.block).'],l2);
                end
            else
                if isfield(stripe,'pmt_index')
                    error('gdc:validation', ...
                        ['Extraneous data field (' stratum_name ...
                        '.stripe{%d}.pmt_index).'],l2);
                end
                if ~isfield(stripe,'block')
                    error('gdc:validation', ...
                        ['Missing data field (' stratum_name ...
                        '.stripe{%d}.block).'],l2);
                end
                if ~(iscell(stripe.block) && isrow(stripe.block) ...
                        && ~isempty(stripe.block))
                    error('gdc:validation', ...
                        ['Wrong data type or size (' stratum_name ...
                        '.stripe{%d}.block).'],l2);
                end
                L3 = length(stripe.block);
                for l3 = 1:L3
                    block = stripe.block{l3};
                    if ~isfield(block,'c2')
                        error('gdc:validation', ...
                            ['Missing data field (' stratum_name ...
                            '.stripe{%d}.block{%d}.c2).'],l2,l3);
                    end
                    [param_size,err_msg] = ...
                        check_param(param_size,block.c2);
                    if ~isempty(err_msg)
                        error('gdc:validation', ...
                            [err_msg ' (' stratum_name ...
                            '.stripe{%d}.block{%d}.c2).'],l2,l3);
                    end
                    if ~check_real(block.c2)
                        error('gdc:validation', ...
                            ['c2 must be real (' stratum_name ...
                            '.stripe{%d}.block{%d}.c2).'],l2,l3);
                    end
                    if ~isfield(block,'pmt_index')
                        error('gdc:validation', ...
                            ['Missing data field (' stratum_name ...
                            '.stripe{%d}.block{%d}.pmt_index).'],l2,l3);
                    end
                    if ~(check_integer(block.pmt_index) && ...
                            isscalar(block.pmt_index))
                        error('gdc:validation', ...
                            ['Wrong data type or size (' stratum_name ...
                            '.stripe{%d}.block{%d}.pmt_index).'],l2,l3);
                    end
                    if block.pmt_index<1 || ...
                            block.pmt_index>length(grating.pmt)
                        error('gdc:validation', ...
                            ['pmt_index is out of range (' ...
                            stratum_name ...
                            '.stripe{%d}.block{%d}.pmt_index).'],l2,l3);
                    end
                end
                c2 = stripe.block{end}.c2-1;
                for l3 = 1:L3
                    c2_ = c2;
                    c2 = stripe.block{l3}.c2;
                    if any(c2_>c2,'all')
                        % violation of Eq 3.37
                        error('gdc:validation', ...
                            ['c2 violates monoticity condition (' ...
                            stratum_name ...
                            '.stripe{%d}.block{%d}.c2).'],l2,l3);
                    end
                end
                clear c2 c2_
            end
        end
        c1 = stratum.stripe{end}.c1-1;
        for l2 = 1:L2
            c1_ = c1;
            c1 = stratum.stripe{l2}.c1;
            if any(c1_>c1,'all')
                % violation of Eq 3.35
                error('gdc:validation', ...
                    ['c1 violates monoticity condition (' ...
                    stratum_name '.stripe{%d}.c1).'],l2);
            end
        end
        clear c1 c1_
    end
    if isfield(stratum,'full_field')
        stratum.full_field = [];
    end
end
if stratum.type~=0
    if isfield(stratum,'pmt_index')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.pmt_index).']);
    end
end
if stratum.type~=2
    if stratum.type~=1
        if stratum.type~=0
            if isfield(stratum,'thick')
                error('gdc:validation', ...
                    ['Extraneous data field (' stratum_name '.thick).']);
            end
            if isfield(stratum,'full_field')
                error('gdc:validation', ...
                    ['Extraneous data field (' stratum_name ...
                    '.full_field).']);
            end
        end
        if isfield(stratum,'h11')
            error('gdc:validation', ...
                ['Extraneous data field (' stratum_name '.h11).']);
        end
        if isfield(stratum,'h12')
            error('gdc:validation', ...
                ['Extraneous data field (' stratum_name '.h12).']);
        end
        if isfield(stratum,'stripe')
            error('gdc:validation', ...
                ['Extraneous data field (' stratum_name '.stripe).']);
        end
    end
    if isfield(stratum,'h21')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.h21).']);
    end
    if isfield(stratum,'h22')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.h22).']);
    end
end
if stratum.type==3
    if ~isfield(stratum,'dx2')
        error('gdc:validation', ...
            ['Missing data field (' stratum_name '.dx2).']);
    end
    [param_size,err_msg] = check_param(param_size,stratum.dx2);
    if ~isempty(err_msg)
        error('gdc:validation', ...
            [err_msg ' (' stratum_name '.dx2).']);
    end
    if ~check_real(stratum.dx2)
        error('gdc:validation', ...
            ['dx2 must be real (' stratum_name '.dx2).']);
    end
    if ~isfield(stratum,'dx3')
        error('gdc:validation', ...
            ['Missing data field (' stratum_name '.dx3).']);
    end
    [param_size,err_msg] = check_param(param_size,stratum.dx3);
    if ~isempty(err_msg)
        error('gdc:validation', ...
            [err_msg ' (' stratum_name '.dx3).']);
    end
    if ~check_real(stratum.dx3)
        error('gdc:validation', ...
            ['dx3 must be real (' stratum_name '.dx3).']);
    end
else
    if isfield(stratum,'dx2')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.dx2).']);
    end
    if isfield(stratum,'dx3')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.dx3).']);
    end
end
if stratum.type==4
    if ~isfield(stratum,'stratum')
        error('gdc:validation', ...
            ['Missing data field (' stratum_name '.stratum).']);
    end
    if ~(iscell(stratum.stratum) && ...
            (isrow(stratum.stratum) || isempty(stratum.stratum)))
        error('gdc:validation', ...
            ['Wrong data type or size (' stratum_name '.stratum).']);
    end
    for l1 = 1:length(stratum.stratum)
        [param_size,stratum.stratum{l1}] = ...
            check_stratum(param_size,grating,stratum.stratum{l1},...
            [stratum_name '.stratum{' num2str(l1) '}']);
    end
    if ~isfield(stratum,'rep_count')
        error('gdc:validation', ...
            ['Missing data field (' stratum_name '.rep_count).']);
    end
    if ~(check_integer(stratum.rep_count) && isscalar(stratum.rep_count))
        error('gdc:validation', ...
            ['Wrong data type or size (' stratum_name '.rep_count).']);
    end
    if (stratum.rep_count<0)
        error('gdc:validation', ...
            ['rep_count must be non-negative (' stratum_name ...
            '.rep_count)']);
    end
else
    if isfield(stratum,'stratum')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.stratum).']);
    end
    if isfield(stratum,'rep_count')
        error('gdc:validation', ...
            ['Extraneous data field (' stratum_name '.rep_count).']);
    end
end

end % check_stratum


function ok = check_real(x)
ok = isnumeric(x) && isreal(x);
end


function ok = check_integer(x)
ok = isnumeric(x) && isreal(x) && all(x(:)==fix(x(:)));
end


function [param_size,err_msg] = check_param(param_size,x)
if ~isnumeric(x)
    err_msg = 'Wrong data type';
    return
end
err_msg = '';
s = size(x);
if isequal(s,param_size)
    return
end
if any(s==0)
    err_msg = 'Empty parameter';
    return
end
s(1,end+1:length(param_size)) = 1;
param_size(1,end+1:length(s)) = 1;
if any(s~=param_size & s~=1 & param_size~=1)
    err_msg = 'Incompatible parameter size';
    return
end
param_size = max(param_size,s);
end

