function [param_size,scat_field,inc_field,grating]=...
    gdc(grating,inc_field,order,show_progress,all_inc_order)
%
% gdc (Grating Diffraction Calculator)
%
% Compute the scattered electromagnetic field from a biperiodic grating.
%
% Syntax:
%
%   Compute the scattered electromagnetic field:
%   [param_size,scat_field,inc_field]=gdc(grating,inc_field,order);
%   ...
%   [param_size,scat_field,inc_field]=gdc(grating,inc_field,order,...
%       show_progress,all_inc_order);
%
%   Just validate input(s); do not run computations:
%   gdc(grating,...);
%   param_size=gdc(grating,...);
%
%   Just display version and copyright information:
%   gdc;
%
% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc ®)
%   Coupled-Wave Theory for Biperiodic DiffractionGratings
%   (GD-Calc.pdf)
%
% Note on parameterization:
%   In the following interface specification "parameters" are
%   multi-dimensional numeric arrays, which must be non-empty and must all
%   be size-matched except for singleton dimensions. (Singleton dimensions
%   are implicitly repmat-extended, if necessary, to match parameter
%   sizes.)
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
%   A replication module's component strata can be of any type,
%   (homogeneous, uniperiodic, biperiodic, coordinate break, or nested
%   replication module).
%
% Note on the progress bar:
%   A progress bar will be displayed and will update once per stratum
%   unless it is disabled by setting show_progress = false. CAUTION: The
%   progress bar can sometimes incur significant runtime overhead. The
%   progress bar may not close if an error occurs in gdc_engine - e.g. if
%   there is an out-of-memory error or if you abort the program via CTRL-C.
%   In this case, do a "clear all" so you can close the progress bar
%   manually.
%
% Note on "Exponential overflow" error:
%   This can sometimes occur with a very thick homogeneous layer. The
%   problem can be circumvented by splitting the thick layer into multiple
%   thin layers, e.g.
%     rep_count = 16;
%     grating.stratum{1}.thick = grating.stratum{1}.thick/rep_count;
%     clear stratum
%     stratum.type = 4;
%     stratum.stratum = {grating.stratum{1}};
%     stratum.rep_count = rep_count;
%     grating.stratum{1} = stratum;
%
% Inputs:
%
%   grating (struct): grating specification
%
%     grating.pmt (size-[1,?] cell array): complex permittivities for
%     grating materials
%
%       grating.pmt{k} (parameter, complex): permittivity for k-th material
%       (must be non-zero; imaginary part must be non-negative.)
%
%     grating.pmt_sub_index, pmt_sup_index (integer): indices into
%     grating.pmt for substrate and superstrate permittivities
%     (GD-Calc.pdf, equations 3.11 and 3.12)
%
%     grating.d21, d31, d22, d32 (parameter, real): grating period vectors
%     (GD-Calc.pdf, equations 3.5 and 3.6) The 1st and 2nd period vectors'
%     [x2,x3] coordinates are [d21,d31] and [d22,d32], respectively. The
%     vectors must be linearly independent (GD-Calc.pdf, condition 3.7).
%
%     grating.stratum (size-[1,L1] cell array): stratum specifications (L1
%     = number of strata, may be zero)
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
%           grating.stratum{l1}.h11, h12 (integer): harmonic indices, not
%           both zero (GD-Calc.pdf, condition 3.21)
%
%           grating.stratum{l1}.stripe (size-[1,L2] cell array): stripe
%           specifications (L2 = number of stripes in stratum l1, must be
%           nonzero)
%
%             grating.stratum{l1}.stripe{l2} (struct): specification for
%             stripe l2 in stratum l1
%
%               grating.stratum{l1}.stripe{l2}.c1 (parameter, real): stripe
%               position parameters; must satisfy monoticity condition
%               (GD-Calc.pdf, condition 3.33) The stripe's positive-side
%               boundary is c1 in period-normalized units.
%
%               grating.stratum{l1}.stripe{l2}.pmt_index (integer): index
%               into grating.pmt for stripe permittivity
%
%         for a biperiodic stratum (stratum type is 2):
%
%           grating.stratum{l1}.thick (parameter, non-negative real):
%           stratum thickness
%
%           grating.stratum{l1}.h11, h12, h21, h22 (integer): harmonic
%           indices; must satisfy GD-Calc.pdf, condition 3.17
%
%           grating.stratum{l1}.stripe (length-L2 cell array): stripe
%           specifications (L2 = number of stripes in stratum l1, must be
%           nonzero)
%
%             grating.stratum{l1}.stripe{l2} (struct): specification for
%             stripe l2 in stratum l1
%
%               grating.stratum{l1}.stripe{l2}.c1 (parameter, real): stripe
%               position parameters; must satisfy monoticity condition
%               (GD-Calc.pdf, condition 3.33) The stripe's positive-side
%               boundary is c1 in period-normalized units.
%
%               grating.stratum{l1}.stripe{l2}.type (0 or 1): type of
%               stripe - 0 for homogeneous, 1 for inhomogeneous
%
%               for a homogeneous stripe (stripe type is 0):
%
%                 grating.stratum{l1}.stripe{l2}.pmt_index (integer): index
%                 into grating.pmt for stripe permittivity
%
%               for an inhomogeneous stripe (stripe type is 1):
%
%                 grating.stratum{l1}.stripe{l2}.block (size-[1,L3] cell
%                 array): structural block specifications (L3 = number of
%                 blocks in stripe l2 of stratum l1, must be nonzero)
%
%                   grating.stratum{l1}.stripe{l2}.block{l3}: specification
%                   for block l3 in stripe l2 of stratum l1
%
%                     grating.stratum{l1}.stripe{l2}.block{l3}.c2
%                     (parameter, real): block position parameter, must
%                     satisfy monoticity condition (GD-Calc.pdf, condition
%                     3.35) The block's positive-side boundary is c2 in
%                     period-normalized units.
%
%                     grating.stratum{l1}.stripe{l2}.block{l3}.pmt_index
%                     (integer): index into grating.pmt for block
%                     permittivity
%
%         for a coordinate break (stratum type is 3):
%
%           grating.stratum{l1}.dx2, dx3 (parameter, real): lateral
%           translational offset (GD-Calc.pdf, equation 3.39) The
%           coordinate break effects a translational displacement of all
%           higher strata by [dx2,dx3].
%
%         for a replication module (stratum type is 4):
%
%           grating.stratum{l1}.stratum (size-[1,?] cell array, may be
%           empty): stratum specifications for module's strata (same data
%           format as grating.stratum)
%
%           grating.stratum{l1}.rep_count (non-negative integer): module
%           replication count
%
%   inc_field (struct): incident field specification
%
%     inc_field.wavelength (parameter, positive real): wavelength (same
%     length units as grating.d21, etc.)
%
%     inc_field.f2, f3 (parameter, real): incident field's tangential
%     spatial frequency vector (reciprocal-length units; GD-Calc.pdf,
%     equations 4.1 and 4.2)
%
%   order (size-[1,?] struct): diffraction orders to be retained in
%   calculations, non-empty
%
%     order(k).m2 (integer): m2 index from M2 (GD-Calc.pdf, condition
%     4.11). The m2 indices must be unique, but need not be sorted. One of
%     the specified m2 indices must be zero.
%
%     order(k).m1 (size-[1,?] array of integer): m1 indices corresponding
%     to m2 (i.e., from M1; GD-Calc.pdf, condition 4.12), non-empty. For
%     each m2, the associated m1 indices must be unique, but need not be
%     sorted. For m2=0, one of the associated m1 indices must be zero.
%
%   show_progress: true or false, OPTIONAL, default is true. If true, a
%   progress bar (waitbar) will be displayed while gdc is executing.
%
%   all_inc_order: true or false, OPTIONAL, default is false. Set to
%   all_inc_order = true to enable multi-wave incident field. (See
%   "Advanced features" below.)
%
% Outputs:
%
%   param_size (size-[1,?] integer): Common parameter size, with full
%   repmat extension. (If the inc_field argument in missing, param_size is
%   determined only on the basis of parameters in the grating input.)
%
%   scat_field (size-[1,?] struct): scattered field. Each element
%   scat_field(k) corresponds to one diffracted order (transmitted and
%   reflected fields)
%
%     scat_field(k).m1, m2 (integer): diffraction order indices
%
%     scat_field(k).f2, f3 (parameter, real): field's grating-tangential
%     spatial frequencies (GD-Calc.pdf, equation 4.7)
%
%     scat_field(k).f1r, f1t (parameter, complex): reflected and
%     transmitted field's grating-normal spatial frequencies (GD-Calc.pdf,
%     equations 4.9, 4.10)
%
%     scat_field(k).s2, s3 (parameter, real): grating-tangential basis
%     vector for S polarization (GD-Calc.pdf, equations 4.23, 4.28)
%
%     scat_field(k).p1t, p2t, p3t (parameter, complex): basis vector for P
%     polarization, for transmitted field (GD-Calc.pdf, equation 4.29)
%
%     scat_field(k).p1r, p2r, p3r (parameter, complex): basis vector for P
%     polarization, for reflected field (GD-Calc.pdf, equation 4.24)
%
%     scat_field(k).Tss, Tsp, Tps, Tpp (parameter, complex): transmission
%     matrix (GD-Calc.pdf, equations 4.32, 4.34)
%
%     scat_field(k).Rss, Rsp, Rps, Rpp (parameter, complex): reflection
%     matrix (GD-Calc.pdf, equations 4.31, 4.33)
%
%   inc_field (struct): incident field specification
%
%     inc_field.wavelength, f2, f3: same as inc_field input argument
%
%     inc_field.f1 (parameter, complex): grating-normal spatial frequency
%     (GD-Calc.pdf, equation 4.8)
%
%     inc_field.s2, s3 (parameter, real): grating-tangential basis vector
%     for S polarization (GD-Calc.pdf, equation 4.18)
%
%     inc_field.p1, p2, p3 (parameter, complex): basis vector for P
%     polarization (GD-Calc.pdf, equation 4.19)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advanced features
%
%
% all_inc_order option:
%
%   The optional input argument all_inc_order is a logical (default =
%   false), which can be set to true to get the scattered field for all
%   incident orders (not just the zero-order incident wave). With
%   all_inc_order = true, the following changes apply:
%
%   The inc_field output struct will be a row vector with elements
%   inc_field(k_inc) corresponding to different incident orders. The order
%   indices corresponding to inc_field(k_inc) will be specified as
%   inc_field(k_inc).m1_inc, inc_field(k_inc).m2_inc.
%
%   The scat_field output struct will be a square matrix with rows
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
%     [param_size,scat_field,inc_field,grating]=...
%         gdc(grating,inc_field,order);
%
%   Make the following initializations in the grating struct (in any
%   combination - all are optional) to calculate the full field in the
%   substrate, in a stratum, or in the superstrate:
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
%   will be set to a size-[1,?] struct containing the E and H fields'
%   Fourier amplitudes ffE and ffH (GD-Calc.pdf, equations 5.13 and 5.14),
%   evaluated at x2 = 0, x3 = 0, and x1 as follows (with b1 defined in
%   GD-Calc.pdf, equations 3.1-3):
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
%   After running gdc the full_field structs will be formatted as follows:
%
%     full_field(k).m1, m2 (integer): diffraction order indices
%
%     full_field(k).wavelength (parameter, real): wavelength
%
%     full_field(k).f2, f3 (parameter, real): field's grating-tangential
%     spatial frequencies (GD-Calc.pdf, equation 5.8)
%
%     full_field(k).ffE1s, ffE2s, ffE3s (parameter, complex): E-field
%     amplitudes for Fourier order (m1, m2) (GD-Calc.pdf, equation 5.13),
%     with incident field s-polarized and of unit amplitude
%
%     full_field(k).ffH1s, ffH2s, ffH3s (parameter, complex): same as
%     above, but for H field
%
%     full_field(k).ffE1p, ffE2p, ffE3p, ffH1p, ffH2p, ffH3p (parameter,
%     complex): same as above, but for p-polarized incident field
%
%   If a stratum of type 0, 1, or 2 inside a replication module contains a
%   full_field data field, no error will be generated but the full_field
%   result will be empty ([]).
%
%   (Note: The full-field option is not currently implemented to work with
%   all_inc_order. Both can be used together, but the full-field data is
%   only generated for the zero-order incident field.)
%

% Version 07/13/2019
% Copyright (c) 2005-2019, Kenneth C. Johnson. All rights reserved.
% KJ Innovation
% https://doi.org/10.24433/CO.7479617.v3

if nargout>4 || (nargin==0 && nargout>0)
    error('Too many output arguments.')
end
if nargin==0
    if ~any(exist('gdc_engine','file')==[2,6])
        error(['gdc_engine is missing or cannot be found.' ...
            ' (Download source: codeocean.com)']);
    end
    gdc_engine;
    return
end
if (nargout>1 && ~any(nargin==3:5)) || ~any(nargin==1:5)
    error('Wrong number of input argument(s).');
end
if ~isstruct(grating) || length(grating)~=1
    error('Wrong data type or size (grating).');
end
if ~isfield(grating,'pmt')
    error('Missing data field (grating.pmt).');
end
if ~iscell(grating.pmt) || isempty(grating.pmt) || ~isrow(grating.pmt)
    error('Wrong data type or size (grating.pmt).');
end
param_size=[1 1];
for k=1:length(grating.pmt)
    [param_size,err_msg]=check_param(param_size,grating.pmt{k});
    if ~isempty(err_msg)
        error([err_msg ' (grating.pmt{%d}).'],k);
    end
    if any(grating.pmt{k}(:)==0) || any(imag(grating.pmt{k}(:))<0)
        error(['Permittivity must be non-zero; '...
            'imaginary part must be non-negative (grating.pmt{%d}).'],k);
    end
end
if ~isfield(grating,'pmt_sub_index')
    error('Missing data field (grating.pmt_sub_index).');
end
if ~check_integer(grating.pmt_sub_index) || ...
        length(grating.pmt_sub_index)~=1
    error('Wrong data type or size (grating.pmt_sub_index).');
end
if grating.pmt_sub_index<1 || grating.pmt_sub_index>length(grating.pmt)
    error('grating.pmt_sub_index is out of range.');
end
if ~isfield(grating,'pmt_sup_index')
    error('Missing data field (grating.pmt_sup_index).');
end
if ~check_integer(grating.pmt_sup_index) || ...
        length(grating.pmt_sup_index)~=1
    error('Wrong data type or size (grating.pmt_sup_index).');
end
if grating.pmt_sup_index<1 || grating.pmt_sup_index>length(grating.pmt)
    error('grating.pmt_sup_index is out of range.');
end
if ~isfield(grating,'d21')
    error('Missing data field (grating.d21).');
end
[param_size,err_msg]=check_param(param_size,grating.d21);
if ~isempty(err_msg)
    error([err_msg ' (grating.d21).']);
end
if ~isfield(grating,'d31')
    error('Missing data field (grating.d31).');
end
[param_size,err_msg]=check_param(param_size,grating.d31);
if ~isempty(err_msg)
    error([err_msg ' (grating.d31).']);
end
if ~isfield(grating,'d22')
    error('Missing data field (grating.d22).');
end
[param_size,err_msg]=check_param(param_size,grating.d22);
if ~isempty(err_msg)
    error([err_msg ' (grating.d22).']);
end
if ~isfield(grating,'d32')
    error('Missing data field (grating.d32).');
end
[param_size,err_msg]=check_param(param_size,grating.d32);
if ~isempty(err_msg)
    error([err_msg ' (grating.d32).']);
end
[dg21,dg31,dg22,dg32]=...
    repmat_extend(grating.d21,grating.d31,grating.d22,grating.d32);
if any(dg21(:).*dg32(:)-dg31(:).*dg22(:)==0)
    % violation of GD-Calc.pdf, equation 3.7
    error('Grating period vectors must be linearly independent.');
end
clear dg21 dg31 dg22 dg32
if isfield(grating,'sub_full_field')
    grating.sub_full_field = [];
end
if isfield(grating,'sup_full_field')
    grating.sup_full_field = [];
end
if ~isfield(grating,'stratum')
    error('Missing data field (grating.stratum).');
end
if ~iscell(grating.stratum) || ...
        (~isempty(grating.stratum) && ~isrow(grating.stratum))
    error('Wrong data type or size (grating.stratum).');
end
L1=length(grating.stratum);
for l1=1:L1
    [param_size,grating.stratum{l1}]=...
        check_stratum(param_size,grating,grating.stratum{l1},...
        ['grating.stratum{' num2str(l1) '}']);
end
if nargin==1
    return
end
if ~isstruct(inc_field) || length(inc_field)~=1
    error('Wrong data type or size (inc_field).');
end
if ~isfield(inc_field,'wavelength')
    error('Missing data field (inc_field.wavelength).');
end
[param_size,err_msg]=check_param(param_size,inc_field.wavelength);
if ~isempty(err_msg)
    error([err_msg ' (inc_field.wavelength).']);
end
if ~check_real(inc_field.wavelength) || any(inc_field.wavelength(:)<=0)
    error('inc_field.wavelength must be real and positive.');
end
if ~isfield(inc_field,'f2')
    error('Missing data field (inc_field.f2).');
end
[param_size,err_msg]=check_param(param_size,inc_field.f2);
if ~isempty(err_msg)
    error([err_msg ' (inc_field.f2).']);
end
if ~check_real(inc_field.f2)
    error('inc_field.f2 must be real.');
end
if ~isfield(inc_field,'f3')
    error('Missing data field (inc_field.f3).');
end
[param_size,err_msg]=check_param(param_size,inc_field.f3);
if ~isempty(err_msg)
    error([err_msg ' (inc_field.f3).']);
end
if ~check_real(inc_field.f3)
    error('inc_field.f3 must be real.');
end
if nargin==2
    return
end
if ~isstruct(order) || ~isrow(order)
    error('Wrong data type or size (order).');
end
if ~isfield(order,'m2')
    error('Missing data field (order(...).m2).');
end
if ~isfield(order,'m1')
    error('Missing data field (order(...).m1).');
end
for k=1:length(order)
    if ~check_integer(order(k).m2) || length(order(k).m2)~=1
        error('Wrong data type or size (order(%d).m2).',k);
    end
    if ~check_integer(order(k).m1) || ~isrow(order(k).m1)
        error('Wrong data type or size (order(%d).m1).',k);
    end
    if length(order(k).m1)~=length(unique(order(k).m1))
        error('order(%d).m1 elements must be unique.',k);
    end
end
if length([order.m2])~=length(unique([order.m2]))
    error('[order.m2] elements must be unique.');
end
k=find([order.m2]==0);
if isempty(k) || ~any(order(k).m1==0)
    error('Missing zero order in specified order indices.');
end
if nargout<=1
    return
end
if ~any(exist('gdc_engine','file')==[2,6])
    error(['gdc_engine is missing or cannot be found.' ...
        ' (Download source: codeocean.com)']);
end
if nargin < 4
    show_progress=true;
else
    if ~islogical(show_progress) || ~isscalar(show_progress)
        error('show_progress must be true or false.');
    end
end
if nargin < 5
    all_inc_order=false;
else
    if ~islogical(all_inc_order) || ~isscalar(all_inc_order)
        error('all_inc_order must be true or false.');
    end
end
if nargout < 4
    [param_size,scat_field,inc_field]=gdc_engine(...
        param_size,grating,inc_field,order,show_progress,all_inc_order);
else
    [param_size,scat_field,inc_field,grating]=gdc_engine(...
        param_size,grating,inc_field,order,show_progress,all_inc_order);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [param_size,stratum]=...
    check_stratum(param_size,grating,stratum,stratum_str)
if ~isstruct(stratum) || length(stratum)~=1
    error(['Wrong data type or size (' stratum_str ').']);
end
if ~isfield(stratum,'type')
    error(['Missing data field (' stratum_str '.type).']);
end
if ~check_integer(stratum.type) || length(stratum.type)~=1
    error(['Wrong data type or size (' stratum_str '.type).']);
end
if ~any(stratum.type==[0 1 2 3 4])
    error(['Stratum type must be 0, 1, 2, 3, or 4 (' ...
        stratum_str '.type).']);
end
if any(stratum.type==[0,1,2])
    if ~isfield(stratum,'thick')
        error(['Missing data field (' stratum_str '.thick).']);
    end
    [param_size,err_msg]=check_param(param_size,stratum.thick);
    if ~isempty(err_msg)
        error([err_msg ' (' stratum_str '.thick).']);
    end
    if ~check_real(stratum.thick) || any(stratum.thick(:)<0)
        error(['Stratum thickness must be real and non-negative (' ...
            stratum_str '.thick).']);
    end
    if stratum.type==0
        if ~isfield(stratum,'pmt_index')
            error(['Missing data field (' stratum_str '.pmt_index).']);
        end
        if ~check_integer(stratum.pmt_index) || ...
                length(stratum.pmt_index)~=1
            error(['Wrong data type or size (' ...
                stratum_str '.pmt_index).']);
        end
        if stratum.pmt_index<1 || stratum.pmt_index>length(grating.pmt)
            error(['pmt_index is out of range (' ...
                stratum_str '.pmt_index).']);
        end
    else % stratum.type==1 or 2
        if ~isfield(stratum,'h11')
            error(['Missing data field (' stratum_str '.h11).']);
        end
        if ~check_integer(stratum.h11) || length(stratum.h11)~=1
            error(['Wrong data type or size (' stratum_str '.h11).']);
        end
        if ~isfield(stratum,'h12')
            error(['Missing data field (' stratum_str '.h12).']);
        end
        if ~check_integer(stratum.h12) || length(stratum.h12)~=1
            error(['Wrong data type or size (' stratum_str '.h12).']);
        end
        if stratum.type==1
            if stratum.h11==0 && stratum.h12==0
                % violation of GD-Calc.pdf, condition 3.21
                error(['h11 and h12 must not both be zero (' ...
                    stratum_str '.h11 and h12).']);
            end
        else
            if ~isfield(stratum,'h21')
                error(['Missing data field (' stratum_str '.h21).']);
            end
            if ~check_integer(stratum.h21) || length(stratum.h21)~=1
                error(['Wrong data type or size (' stratum_str '.h21).']);
            end
            if ~isfield(stratum,'h22')
                error(['Missing data field (' stratum_str '.h22).']);
            end
            if ~check_integer(stratum.h22) || length(stratum.h22)~=1
                error(['Wrong data type or size (' stratum_str '.h22).']);
            end
            if stratum.h11*stratum.h22==stratum.h12*stratum.h21
                % violation of GD-Calc.pdf, condition 3.17
                error(['h11, h12, h21, and h22 must satisfy ' ...
                    'h11*h22~=h12*h21 (' ...
                    stratum_str '.h11, h12, h21, and h22).']);
            end
        end
        if ~isfield(stratum,'stripe')
            error(['Missing data field (' stratum_str '.stripe).']);
        end
        if ~iscell(stratum.stripe) || ~isrow(stratum.stripe)
            error(['Wrong data type or size (' stratum_str '.stripe).']);
        end
        L2=length(stratum.stripe);
        for l2=1:L2
            stripe=stratum.stripe{l2};
            if ~isfield(stripe,'c1')
                error(['Missing data field (' ...
                    stratum_str '.stripe{%d}.c1).'],l2);
            end
            [param_size,err_msg]=check_param(param_size,stripe.c1);
            if ~isempty(err_msg)
                error([err_msg ' (' stratum_str '.stripe{%d}.c1).'],l2);
            end
            if ~check_real(stripe.c1)
                error(['c1 must be real (' ...
                    stratum_str '.stripe{%d}.c1).'],l2);
            end
            if stratum.type==2
                if ~isfield(stripe,'type')
                    error(['Missing data field (' ...
                        stratum_str '.stripe{%d}.type).'],l2);
                end
                if ~check_integer(stripe.type) || length(stripe.type)~=1
                    error(['Wrong data type or size (' ...
                        stratum_str '.stripe{%d}.type).'],l2);
                end
                if ~any(stripe.type==[0 1])
                    error(['Stripe type must be 0 or 1 (' ...
                        stratum_str '.stripe{%d}.type).'],l2);
                end
            else
                if isfield(stripe,'type')
                    error(['Extraneous data field (' ...
                        stratum_str '.stripe{%d}.type).'],l2);
                end
            end
            if stratum.type==1 || stripe.type==0
                if ~isfield(stripe,'pmt_index')
                    error(['Missing data field (' ...
                        stratum_str '.stripe{%d}.pmt_index).'],l2);
                end
                if ~check_integer(stripe.pmt_index) || ...
                        length(stripe.pmt_index)~=1
                    error(['Wrong data type or size (' ...
                        stratum_str '.stripe{%d}.pmt_index).'],l2);
                end
                if stripe.pmt_index<1 || ...
                        stripe.pmt_index>length(grating.pmt)
                    error(['pmt_index is out of range (' ...
                        stratum_str '.stripe{%d}.pmt_index).'],l2);
                end
                if isfield(stripe,'block')
                    error(['Extraneous data field (' ...
                        stratum_str '.stripe{%d}.block).'],l2);
                end
            else
                if isfield(stripe,'pmt_index')
                    error(['Extraneous data field (' ...
                        stratum_str '.stripe{%d}.pmt_index).'],l2);
                end
                if ~isfield(stripe,'block')
                    error(['Missing data field (' ...
                        stratum_str '.stripe{%d}.block).'],l2);
                end
                if ~iscell(stripe.block) || ~isrow(stripe.block)
                    error(['Wrong data type or size (' ...
                        stratum_str '.stripe{%d}.block).'],l2);
                end
                L3=length(stripe.block);
                for l3=1:L3
                    block=stripe.block{l3};
                    if ~isfield(block,'c2')
                        error(['Missing data field (' stratum_str ...
                            '.stripe{%d}.block{%d}.c2).'],l2,l3);
                    end
                    if l3==1 || ~isequal(size(block.c2),size_)
                        [param_size,err_msg]=...
                            check_param(param_size,block.c2);
                        if ~isempty(err_msg)
                            error([err_msg ' (' stratum_str ...
                                '.stripe{%d}.block{%d}.c2).'],l2,l3);
                        end
                        size_=size(block.c2);
                    end
                    if ~check_real(block.c2)
                        error(['c2 must be real (' stratum_str ...
                            '.stripe{%d}.block{%d}.c2).'],l2,l3);
                    end
                    if ~isfield(block,'pmt_index')
                        error(['Missing data field (' stratum_str ...
                            '.stripe{%d}.block{%d}.pmt_index).'],l2,l3);
                    end
                    if ~check_integer(block.pmt_index) || ...
                            length(block.pmt_index)~=1
                        error(['Wrong data type or size (' stratum_str ...
                            '.stripe{%d}.block{%d}.pmt_index).'],l2,l3);
                    end
                    if block.pmt_index<1 || ...
                            block.pmt_index>length(grating.pmt)
                        error(['pmt_index is out of range (' ...
                            stratum_str ...
                            '.stripe{%d}.block{%d}.pmt_index).'],l2,l3);
                    end
                end
                c2=stripe.block{end}.c2-1;
                for l3=1:L3
                    c2_=c2;
                    c2=stripe.block{l3}.c2;
                    if isequal(size(c2),size(c2_))
                        if any(c2_(:)>c2(:))
                            % violation of GD-Calc.pdf, condition 3.35
                            error(['c2 violates monoticity condition (' ...
                                stratum_str ...
                                '.stripe{%d}.block{%d}.c2).'],l2,l3);
                        end
                    else
                        [xc2_,xc2]=repmat_extend(c2_,c2);
                        if any(xc2_(:)>xc2(:))
                            % violation of GD-Calc.pdf, condition 3.35
                            error(['c2 violates monoticity condition (' ...
                                stratum_str ...
                                '.stripe{%d}.block{%d}.c2).'],l2,l3);
                        end
                    end
                end
                clear c2 c2_ xc2 xc2_
            end
        end
        c1=stratum.stripe{end}.c1-1;
        for l2=1:L2
            c1_=c1;
            c1=stratum.stripe{l2}.c1;
            if isequal(size(c1),size(c1_))
                if any(c1_(:)>c1(:))
                    % violation of GD-Calc.pdf, condition 3.33
                    error(['c1 violates monoticity condition (' ...
                        stratum_str '.stripe{%d}.c1).'],l2);
                end
            else
                [xc1_,xc1]=repmat_extend(c1_,c1);
                if any(xc1_(:)>xc1(:))
                    % violation of GD-Calc.pdf, condition 3.33
                    error(['c1 violates monoticity condition (' ...
                        stratum_str '.stripe{%d}.c1).'],l2);
                end
            end
        end
        clear c1 c1_ xc1 xc1_
    end
    if isfield(stratum,'full_field')
        stratum.full_field=[];
    end
end
if stratum.type~=0
    if isfield(stratum,'pmt_index')
        error(['Extraneous data field (' stratum_str '.pmt_index).']);
    end
end
if stratum.type~=2
    if stratum.type~=1
        if stratum.type~=0
            if isfield(stratum,'thick')
                error(['Extraneous data field (' stratum_str '.thick).']);
            end
            if isfield(stratum,'full_field')
                error(['Extraneous data field (' stratum_str '.full_field).']);
            end
        end
        if isfield(stratum,'h11')
            error(['Extraneous data field (' stratum_str '.h11).']);
        end
        if isfield(stratum,'h12')
            error(['Extraneous data field (' stratum_str '.h12).']);
        end
        if isfield(stratum,'stripe')
            error(['Extraneous data field (' stratum_str '.stripe).']);
        end
    end
    if isfield(stratum,'h21')
        error(['Extraneous data field (' stratum_str '.h21).']);
    end
    if isfield(stratum,'h22')
        error(['Extraneous data field (' stratum_str '.h22).']);
    end
end
if stratum.type==3
    if ~isfield(stratum,'dx2')
        error(['Missing data field (' stratum_str '.dx2).']);
    end
    [param_size,err_msg]=check_param(param_size,stratum.dx2);
    if ~isempty(err_msg)
        error([err_msg ' (' stratum_str '.dx2).']);
    end
    if ~check_real(stratum.dx2)
        error(['dx2 must be real (' stratum_str '.dx2).']);
    end
    if ~isfield(stratum,'dx3')
        error(['Missing data field (' stratum_str '.dx3).']);
    end
    [param_size,err_msg]=check_param(param_size,stratum.dx3);
    if ~isempty(err_msg)
        error([err_msg ' (' stratum_str '.dx3).']);
    end
    if ~check_real(stratum.dx3)
        error(['dx3 must be real (' stratum_str '.dx3).']);
    end
else
    if isfield(stratum,'dx2')
        error(['Extraneous data field (' stratum_str '.dx2).']);
    end
    if isfield(stratum,'dx3')
        error(['Extraneous data field (' stratum_str '.dx3).']);
    end
end
if stratum.type==4
    if ~isfield(stratum,'stratum')
        error(['Missing data field (' stratum_str '.stratum).']);
    end
    if ~iscell(stratum.stratum) || ...
            (~isempty(stratum.stratum) && ~isrow(stratum.stratum))
        error(['Wrong data type or size (' stratum_str '.stratum).']);
    end
    for l1=1:length(stratum.stratum)
        [param_size,stratum.stratum{l1}]=...
            check_stratum(param_size,grating,stratum.stratum{l1},...
            [stratum_str '.stratum{' num2str(l1) '}']);
    end
    if ~isfield(stratum,'rep_count')
        error(['Missing data field (' stratum_str '.rep_count).']);
    end
    if ~check_integer(stratum.rep_count) || length(stratum.rep_count)~=1
        error(['Wrong data type or size (' stratum_str '.rep_count).']);
    end
    if (stratum.rep_count<0)
        error(['rep_count must be non-negative (' ...
            stratum_str '.rep_count)']);
    end
else
    if isfield(stratum,'stratum')
        error(['Extraneous data field (' stratum_str '.stratum).']);
    end
    if isfield(stratum,'rep_count')
        error(['Extraneous data field (' stratum_str '.rep_count).']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ok=check_real(x)
ok=(isnumeric(x) && isreal(x));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ok=check_integer(x)
ok=(isnumeric(x) && isreal(x) && all(x(:)==fix(x(:))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [param_size,err_msg]=check_param(param_size,x)
err_msg='';
s=size(x);
if isequal(s,param_size)
    if ~isnumeric(x)
        err_msg='Wrong data type';
    end
    return
end
if any(s==0)
    err_msg='Empty parameter';
    return
end
if ~isnumeric(x)
    err_msg='Wrong data type';
    return
end
s(1,end+1:length(param_size))=1;
param_size(1,end+1:length(s))=1;
if any(s~=param_size & s~=1 & param_size~=1)
    err_msg='Mismatched parameter size';
    return
end
param_size=max(param_size,s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout]=repmat_extend(varargin)
% [A1,A2,...]=repmat_extend(A1,A2,...);
% Repmat-extend arguments' singleton dimensions to match array sizes.
% (Scalars will not be repmat-extended.) Size compatibility of input
% arguments is assumed to have already been checked (via check_param).
n=nargin; % >=nargout
match=true;
s=size(varargin{1});
for j=2:n
    if ~isequal(size(varargin{j}),s)
        match=false;
        break
    end
end
if match
    n=nargout;
    varargout=cell(1,n);
    for j=1:n
        varargout{j}=varargin{j};
    end
    return
end
s=cell(1,n);
s_=[1,1];
for j=1:n
    s{j}=size(varargin{j});
    len=length(s{j});
    s_(1,end+1:len)=1;
    s_(1:len)=max(s_(1:len),s{j});
end
ndims_=length(s_);
n=nargout;
varargout=cell(1,n);
for j=1:n
    if all(s{j}==1)
        varargout{j}=varargin{j}; % (Do not repmat-extend scalars.)
    else
        s{j}(1,end+1:ndims_)=1;
        ext=s_./s{j};
        if all(ext==1)
            varargout{j}=varargin{j};
        else
            varargout{j}=repmat(varargin{j},ext);
        end
    end
end
