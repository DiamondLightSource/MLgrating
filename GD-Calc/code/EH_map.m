function map = EH_map(full_field,param_index,x2,x3)
% EH_map
%
% Convert an electromagnetic field's spatial frequency distribution into a
% spatial amplitude map.
%
% Syntax:
%
%   map = EH_map(full_field,param_index,x2,x3);
%
%
% Inputs:
%
%   full_field (struct row vector, nonempty): the electromagnetic field's
%   frequency distribution at a particular height in a grating; returned by
%   gdc.m as grating.sub_full_field (substrate), grating.sup_full_field
%   (superstrate), or grating.stratum{l1}.full_field (in a grating stratum,
%   at the top interface); see gdc.m comment header under "Advanced
%   features".
%
%   param_index (integer row vector,nonempty): multidimensional parameter
%   indices for selecting a parameter combination. param_index must be
%   compatible with parameter sizes (with implicit repmat extension of
%   singleton dimensions).
%
%   x2, x3 (real array, nonempty): x2 and x3 coordinates at which the field
%   will be calculated. x2 and x3 must be size-compatible (with implicit
%   repmat extension of singleton dimensions).
%
% Output:
%
%   map (struct array): spatial sampling of electromagnetic field at x2,
%   x3. size(map) matches max(size(x2),size(x3)), and the elements map(...)
%   include the following struct fields representing the electromagnetic
%   field at coordinated x2(...), x3(...):
%
%     map(...).x2, x3 (scalar real double): spatial sampling coordinates
%
%     map(...).E1s, E2s, E3s (scalar complex double): E field vector for
%     s-polarized incident beam
%
%     map(...).H1s, H2s, H3s: same as above for H field
%
%     map(...).E1p, E2p, E3p, H1p, H2p, H3p: same as above for p-polarized
%     incident beam
%
% Note: If not all of the E and H field results are needed, then remove
% the corresponding fields from full_field to avoid unnecessary
% computations. For example, the following code excerpt will compute only
% the E field quantities:
%   map = EH_map(...
%       rmfield(full_field,...
%       {'ffH1s','ffH2s','ffH3s','ffH1p','ffH2p','ffH3p'}),...
%       param_index,x2,x3);
% The map result will not contain the data fields map.H1s, map.H2s,
% map.H3s, map.H1p, map.H2p, map.H3p.
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc)
%   Coupled-Wave Theory for Biperiodic DiffractionGratings
%   (GD-Calc.pdf, version 04-Jun-2022)
% All equation ("Eq") references below are from GD-Calc.pdf.

% Validate inputs
narginchk(4,4)
if isempty(full_field) || ...
        ~(isstruct(full_field) && isrow(full_field) && ...
        all(isfield(full_field,{'f2','f3'})))
    error('EH_map:validation', ...
        'Wrong type or size, or missing data field (full_field).');
end
flagE1s = isfield(full_field,'ffE1s');
flagE2s = isfield(full_field,'ffE2s');
flagE3s = isfield(full_field,'ffE3s');
flagE1p = isfield(full_field,'ffE1p');
flagE2p = isfield(full_field,'ffE2p');
flagE3p = isfield(full_field,'ffE3p');
flagH1s = isfield(full_field,'ffH1s');
flagH2s = isfield(full_field,'ffH2s');
flagH3s = isfield(full_field,'ffH3s');
flagH1p = isfield(full_field,'ffH1p');
flagH2p = isfield(full_field,'ffH2p');
flagH3p = isfield(full_field,'ffH3p');
if isempty(param_index) || ...
        ~(isnumeric(param_index) && isreal(param_index) && ...
        all(param_index(:)==fix(param_index(:))) && ...
        all(param_index>=1) && isrow(param_index))
    error('EH_map:validation', ...
        'Wrong data type or size (param_index).');
end
if isempty(x2) || ~(isnumeric(x2) && isreal(x2))
    error('EH_map:validation', ...
        'Wrong data type or size (x2).');
end
if isempty(x3) || ~(isnumeric(x3) && isreal(x3))
    error('EH_map:validation', ...
        'Wrong data type or size (x3).');
end
for k = 1:length(full_field)
    [full_field(k).f2,err_msg] = ...
        select_param(param_index,full_field(k).f2);
    if ~isempty(err_msg)
        error('EH_map:validation', ...
            [err_msg ' (full_field(%d).f2).'],k);
    end
    [full_field(k).f3,err_msg] = ...
        select_param(param_index,full_field(k).f3);
    if ~isempty(err_msg)
        error('EH_map:validation', ...
            [err_msg ' (full_field(%d).f3).'],k);
    end
    if flagE1s
        [full_field(k).ffE1s,err_msg] = ...
            select_param(param_index,full_field(k).ffE1s);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffE1s).'],k);
        end
    end
    if flagE2s
        [full_field(k).ffE2s,err_msg] = ...
            select_param(param_index,full_field(k).ffE2s);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffE2s).'],k);
        end
    end
    if flagE3s
        [full_field(k).ffE3s,err_msg] = ...
            select_param(param_index,full_field(k).ffE3s);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffE3s).'],k);
        end
    end
    if flagH1s
        [full_field(k).ffH1s,err_msg] = ...
            select_param(param_index,full_field(k).ffH1s);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffH1s).'],k);
        end
    end
    if flagH2s
        [full_field(k).ffH2s,err_msg] = ...
            select_param(param_index,full_field(k).ffH2s);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffH2s).'],k);
        end
    end
    if flagH3s
        [full_field(k).ffH3s,err_msg] = ...
            select_param(param_index,full_field(k).ffH3s);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffH3s).'],k);
        end
    end
    if flagE1p
        [full_field(k).ffE1p,err_msg] = ...
            select_param(param_index,full_field(k).ffE1p);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffE1p).'],k);
        end
    end
    if flagE2p
        [full_field(k).ffE2p,err_msg] = ...
            select_param(param_index,full_field(k).ffE2p);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffE2p).'],k);
        end
    end
    if flagE3p
        [full_field(k).ffE3p,err_msg] = ...
            select_param(param_index,full_field(k).ffE3p);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffE3p).'],k);
        end
    end
    if flagH1p
        [full_field(k).ffH1p,err_msg] = ...
            select_param(param_index,full_field(k).ffH1p);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffH1p).'],k);
        end
    end
    if flagH2p
        [full_field(k).ffH2p,err_msg] = ...
            select_param(param_index,full_field(k).ffH2p);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffH2p).'],k);
        end
    end
    if flagH3p
        [full_field(k).ffH3p,err_msg] = ...
            select_param(param_index,full_field(k).ffH3p);
        if ~isempty(err_msg)
            error('EH_map:validation', ...
                [err_msg ' (full_field(%d).ffH3p).'],k);
        end
    end
end
size_x2 = size(x2);
size_x3 = size(x3);
size_x2(1,end+1:length(size_x3)) = 1;
size_x3(1,end+1:length(size_x2)) = 1;
if ~all(size_x2==1 | size_x3==1 | size_x2==size_x3)
    error('EH_map:validation', ...
        'Size mismatch (x2, x3).');
end
size_x = max(size_x2,size_x3);
x2 = repmat(x2,size_x./size_x2);
x3 = repmat(x3,size_x./size_x3);
f2 = [full_field.f2]; % row vector
f3 = [full_field.f3];
map.x2 = 0;
map.x3 = 0;
if flagE1s
    map.E1s = 0;
    ffE1s = [full_field.ffE1s].'; % column vector
end
if flagE2s
    map.E2s = 0;
    ffE2s = [full_field.ffE2s].';
end
if flagE3s
    map.E3s = 0;
    ffE3s = [full_field.ffE3s].';
end
if flagH1s
    map.H1s = 0;
    ffH1s = [full_field.ffH1s].';
end
if flagH2s
    map.H2s = 0;
    ffH2s = [full_field.ffH2s].';
end
if flagH3s
    map.H3s = 0;
    ffH3s = [full_field.ffH3s].';
end
if flagE1p
    map.E1p = 0;
    ffE1p = [full_field.ffE1p].';
end
if flagE2p
    map.E2p = 0;
    ffE2p = [full_field.ffE2p].';
end
if flagE3p
    map.E3p = 0;
    ffE3p = [full_field.ffE3p].';
end
if flagH1p
    map.H1p = 0;
    ffH1p = [full_field.ffH1p].';
end
if flagH2p
    map.H2p = 0;
    ffH2p = [full_field.ffH2p].';
end
if flagH3p
    map.H3p = 0;
    ffH3p = [full_field.ffH3p].';
end
map = repmat(map,size_x);
for j = 1:prod(size_x)
    x2_ = x2(j);
    x3_ = x3(j);
    map(j).x2 = x2_;
    map(j).x3 = x3_;
    % Eqs 5.13-14:
    exp_ = exp(2i*pi*(f2*x2_+f3*x3_)); % row vector
    if flagE1s
        map(j).E1s = exp_*ffE1s; % scalar
    end
    if flagE2s
        map(j).E2s = exp_*ffE2s;
    end
    if flagE3s
        map(j).E3s = exp_*ffE3s;
    end
    if flagH1s
        map(j).H1s = exp_*ffH1s;
    end
    if flagH2s
        map(j).H2s = exp_*ffH2s;
    end
    if flagH3s
        map(j).H3s = exp_*ffH3s;
    end
    if flagE1p
        map(j).E1p = exp_*ffE1p;
    end
    if flagE2p
        map(j).E2p = exp_*ffE2p;
    end
    if flagE3p
        map(j).E3p = exp_*ffE3p;
    end
    if flagH1p
        map(j).H1p = exp_*ffH1p;
    end
    if flagH2p
        map(j).H2p = exp_*ffH2p;
    end
    if flagH3p
        map(j).H3p = exp_*ffH3p;
    end
end
end % EH_map


function [x,err_msg] = select_param(param_index,x)
err_msg = '';
s = size(x);
if any(s==0)
    err_msg = 'Empty array';
    return
end
s(1,end+1:length(param_index)) = 1;
param_index(1,end+1:length(s)) = 1;
if any(s~=1 & param_index>s)
    err_msg = 'param_index exceeds array size';
    return
end
param_index(s==1) = 1;
param_index = num2cell(param_index);
x = x(param_index{:});
end % select_param
