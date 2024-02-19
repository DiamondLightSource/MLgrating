function [varargout] = ...
    gdc_plot(grating,param_index,pmt_display,x_limit,new_fig,max_elements)
% gdc_plot
%
% Plot a biperiodic grating structure (for use with gdc).
%
% Syntax:
%
%   gdc_plot(grating);
%   gdc_plot(grating,param_index,pmt_display,x_limit,new_fig,max_elements);
%   [h_plot,param_index,pmt_display,x_limit] = gdc_plot(grating,...);
%
% Inputs: (Any OPTIONAL input can be [], which converts to the default.)
%
%   grating: same as first gdc input
%
%   param_index (integer row vector, OPTIONAL, default = 1):
%   multidimensional parameter indices for selecting which parameter
%   combination to display. param_index must be compatible with parameter
%   sizes, with implicit repmat extension of singleton dimensions. (Run
%   param_size = gdc(grating) to get parameter sizes; ensure that
%   param_index(j)>=1 and either param_index(j)<=param_size(j) or
%   param_size(j)==1.)
%
%   pmt_display (struct row vector, OPTIONAL): display properties for
%   optical materials (permittivities). pmt_display must be size-matched to
%   grating.pmt, and comprises the following fields (all required if
%   pmt_display is specified):
%
%     pmt_display(m).name (string): legend name for material grating.pmt{m}
%     If the string is empty ('') the material will not be shown in the
%     legend. The default name is '' (no legend) if pmt_display is
%     unspecified.
%
%     pmt_display(m).color (either size-[1,3], size-2,3] real, or []): RGB
%     color for displaying material grating.pmt{m}, or [] to suppress
%     displaying the material. Values must be in the range 0 <= RGB <= 1.
%     The default color is [.5,.5,.5] (gray) unless grating.pmt{m} is 1, in
%     which case the default is []. A size-[1,3] color applies to graphic
%     faces, while edges will be black. A size-[2,3] color applies
%     color(1,:) to faces and color(2,:) to edges.
%
%     pmt_display(m).alpha (scalar or length-2 real): transparency value
%     for material grating.pmt{m} The value must be in the range 0 <= alpha
%     <= 1. (0 means transparent; 1 means opaque.) The default alpha is
%     0.5. A scalar alpha applies to graphic faces and edges; a length-2
%     alpha applies alpha(1) to faces and alpha(2) to edges.
%
%   x_limit (size-[2,3] real, OPTIONAL): limits of display rectangle:
%   x_limit = [x1_min, x2_min, x3_min; x1_max, x2_max, x3_max], with
%   x_limit(1,:) < x_limit(2,:). (The grating substrate is the half-space
%   x1<=0.) The default x_limit covers approximately 2 grating periods.
%
%   new_fig (true or false, OPTIONAL, default = true): Set new_fig to false
%   to plot into the current figure (gcf) and axis (gca); otherwise a new
%   figure will be created.
%
%   max_elements (positive integer, OPTIONAL, default = 1000): maximum
%   number of grating structure elements (strata, stripes, blocks) that can
%   be plotted without user confirmation. If the number is exceeded the
%   user will be prompted to continue or cancel the plot.
%
% Outputs:
%
%   h_plot (handle or []): handle to plot window, or [] if the plot was
%   aborted
%
%   param_index, pmt_display ,x_limit: Echo inputs or default arguments.
%   (These values can be selectively edited to modify the plot.)
%
% Notes:
%
%   The default persective is view(3). Use view(2) for a plan view (x2, x3
%   projection).
%
%   Set "axis image" to fit the plot box tightly around the image.
%
%   If the plot has an extreme aspect ratio, set "axis normal" to compress
%   the scale.
%
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innvoation https://kjinnovation.com/

% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc)
%   Coupled-Wave Theory for Biperiodic DiffractionGratings
%   (GD-Calc.pdf, version 04-Jun-2022)
% All equation ("Eq") references below are from GD-Calc.pdf.

narginchk(1,6)
nargoutchk(0,4)

% Validate inputs
if nargin<2 || isequaln(param_index,[])
    % Default param_index:
    param_index = 1;
else
    param_size = gdc(grating);
    if ~(isnumeric(param_index) && isreal(param_index) && ...
            isequal(param_index,fix(param_index)) && isrow(param_index))
        error('gdc_plot:validation', ...
            'Wrong data type or size (param_index).');
    end
    param_size(end+1:length(param_index)) = 1;
    param_index(end+1:length(param_size)) = 1;
    param_index(param_size==1 & param_index>1) = 1;
    if any(param_index<=0) || any(param_index>param_size)
        error('gdc_plot:validation', ...
            'param_index is out of range.');
    end
end
if nargin<3 || isequaln(pmt_display,[])
    % Default pmt_display:
    pmt_display = repmat( ...
        struct('name','','color',[.5,.5,.5],'alpha',.5), ...
        size(grating.pmt));
    for j = 1:length(grating.pmt)
        if isequal(grating.pmt{j},1)
            pmt_display(j).color = [];
        end
    end
else
    if ~isstruct(pmt_display) || ...
            ~isequal(size(pmt_display),size(grating.pmt))
        error('gdc_plot:validation', ...
            'Wrong data type or size (pmt_display).');
    end
    if ~isfield(pmt_display,'name')
        error('gdc_plot:validation', ...
            'Missing data field (pmt_display(m).name).');
    end
    if ~isfield(pmt_display,'color')
        error('gdc_plot:validation', ...
            'Missing data field (pmt_display(m).color).');
    end
    if ~isfield(pmt_display,'alpha')
        error('gdc_plot:validation', ...
            'Missing data field (pmt_display(m).alpha).');
    end
    for m = 1:length(pmt_display)
        if ~ischar(pmt_display(m).name)
            error('gdc_plot:validation', ...
                'Wrong data type (pmt_display(%d).name).',m);
        end
        if ~isequal(pmt_display(m).color,[])
            if ~(isnumeric(pmt_display(m).color) && ...
                    isreal(pmt_display(m).color) && ...
                    (isequal(size(pmt_display(m).color),[1,3]) || ...
                    isequal(size(pmt_display(m).color),[2,3])))
                error('gdc_plot:validation', ...
                    ['Wrong data type or size ' ...
                    '(pmt_display(%d).color).'],m);
            end
            if min(pmt_display(m).color(:))<0 || ...
                    max(pmt_display(m).color(:))>1
                error('gdc_plot:validation', ...
                    ['pmt_display(%d).color must be ' ...
                    'in the range 0 to 1.'],m);
            end
        end
        if ~(isnumeric(pmt_display(m).alpha) && ...
                isreal(pmt_display(m).alpha) && ...
                (isscalar(pmt_display(m).alpha) || ...
                numel(pmt_display(m).alpha)==2))
            error('gdc_plot:validation', ...
                'Wrong data type or size (pmt_display(%d).alpha).',m);
        end
        if any(pmt_display(m).alpha<0) || any(pmt_display(m).alpha>1)
            error('gdc_plot:validation', ...
                'pmt_display(%d).alpha must be in the range 0 to 1.',m);
        end
    end
end

% Grating periods:
d21 = param(grating.d21,param_index);
d31 = param(grating.d31,param_index);
d22 = param(grating.d22,param_index);
d32 = param(grating.d32,param_index);

if nargin<4 || isequaln(x_limit,[])
    % Default x_limit:
    t = 0;
    for j = 1:length(grating.stratum)
        t = t+stratum_thick(grating.stratum{j},param_index);
    end
    x_limit = [t,max(abs([d21,d22])),max(abs([d31,d32]))];
    x_limit(1) = max(x_limit);
    x_limit = [x_limit.*[-.1,-1,-1];x_limit.*[1.1,1,1]];
else
    if ~(isnumeric(x_limit) && isreal(x_limit) && ...
            isequal(size(x_limit),[2,3]))
        error('gdc_plot:validation', ...
            'Wrong data type or size (x_limit)');
    end
    if ~all(x_limit(1,:)<x_limit(2,:))
        error('gdc_plot:validation', ...
            'x_limit must satisfy x_limit(1,:)<x_limit(2,:).')
    end
end

if nargin<5 || isequaln(new_fig,[])
    new_fig = true;
elseif ~(islogical(new_fig) && isscalar(new_fig))
    error('gdc_plot:validation', ...
        'new_fig must be true or false.');
end

if nargin<6 || isequaln(max_elements,[])
    max_elements = 1000;
elseif ~(isnumeric(max_elements) && isreal(max_elements) && ...
        isscalar(max_elements) && max_elements==fix(max_elements) && ...
        max_elements>0)
    error('gdc_plot:validation', ...
        'Invalid input (max_elements).');
end

n = nargout;
varargout = cell(1,n);
if n>=2
    varargout{2} = param_index;
    if n>=3
        varargout{3} = pmt_display;
        if n>=4
            varargout{4} = x_limit;
        end
    end
end

% Warn user if there are too many structure elements in the figure.
N = 0;
for l1 = 1:length(grating.stratum)
    N = N+num_elements(grating.stratum{l1});
end
if N>max_elements && ~isequal( ...
        questdlg(['WARNING: There are ' num2str(N) ...
        ' structure elements. OK to proceed with plot?'], ...
        '','Yes','No','Yes'), ...
        'Yes')
    return
end

% Initialize figure/axis.
if new_fig
    h_plot = figure;
    vis = 'on';
else
    h_plot = gcf;
    vis = get(h_plot,'Visible');
    cla(gca,'reset')
end
set(h_plot,'Visible','off');

if n>=1
    varargout{1} = h_plot;
end

axis([...
    x_limit(1,2),x_limit(2,2),...
    x_limit(1,3),x_limit(2,3),...
    x_limit(1,1),x_limit(2,1)...
    ]);
axis equal
xlabel('x_2');
ylabel('x_3');
zlabel('x_1');
hold on;

% Allocate patch handles for generating legend.
h = cell(size(grating.pmt));

% x_eps = display clearance between structural blocks.
x_eps = max(x_limit(2,:)-x_limit(1,:))*1.0e-5;

% x_max = lateral radius enclosing x_limit
x_max = sqrt(max(x_limit(:,2).^2)+max(x_limit(:,3).^2))+x_eps;

% Display the substrate.
m = grating.pmt_sub_index;
if ~isempty(pmt_display(m).color)
    x0 = x_limit(1,:)+x_eps;
    edge1 = [-x0(1),0,0];
    if edge1(1)>0
        edge2 = [0,x_limit(2,2)-x_eps-x0(2),0];
        edge3 = [0,0,x_limit(2,3)-x_eps-x0(3)];
        h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
            pmt_display(m).color,pmt_display(m).alpha);
        if ~isempty(h_)
            h{m} = h_;
        end
    end
end

% Loop through strata.
x1_top = 0; % cumulative stratum thicknesses
dx2 = 0; % cumulative coordinate breaks
dx3 = 0;
L1 = length(grating.stratum);
for l1 = 1:L1
    [x1_top,dx2,dx3,h] = plot_stratum(...
        x1_top,dx2,dx3,h,grating.stratum{l1},...
        param_index,pmt_display,x_limit,x_eps,x_max,d21,d31,d22,d32);
end

% Construct legend.
h_ = {};
name = {};
for m = 1:length(h)
    if ~isempty(h{m}) && ~isempty(pmt_display(m).name)
        h_{end+1} = h{m}; %#ok<AGROW>
        name{end+1} = pmt_display(m).name; %#ok<AGROW>
    end
end
if ~isempty(h_)
    legend([h_{:}],name{:});
end
view(3)

drawnow
set(h_plot,'Visible',vis);

end % gdc_plot


function p = param(p,i)
% Index into multi-dimensional parameter p with subscript list i; return
% scalar parameter. (Repmat expansion of singleton dimensions in p is
% implicit.)
size_ = size(p);
n = length(size_);
i(n+1:end) = [];
i(end+1:n) = 1;
% length(i)==n
i(size_==1) = 1;
i = num2cell(i);
p = p(i{:});
end % param


function t = stratum_thick(stratum,param_index)
switch stratum.type
    case {0,1,2}
        t = param(stratum.thick,param_index);
    case 4
        t = 0;
        for j = 1:length(stratum.stratum)
            t = t+stratum_thick(stratum.stratum{j},param_index);
        end
        t = t*stratum.rep_count;
    otherwise
        t = 0;
end
end % stratum_thick


function N = num_elements(stratum)
if stratum.type==0
    N = 1;
elseif stratum.type==1
    N = length(stratum.stripe);
elseif stratum.type==2
    N = 0;
    for l2 = 1:length(stratum.stripe)
        if stratum.stripe{l2}.type==0
            N = N+1;
        else
            N = N+length(stratum.stripe{l2}.block);
        end
    end
elseif stratum.type==3
    N = 0;
elseif stratum.type==4
    N = 0;
    for l1 = 1:length(stratum.stratum)
        N = N+num_elements(stratum.stratum{l1});
    end
    N = stratum.rep_count*N;
end
end % num_elements


function [x1_top,dx2,dx3,h] = plot_stratum(x1_top,dx2,dx3,h,stratum,...
    param_index,pmt_display,x_limit,x_eps,x_max,d21,d31,d22,d32)
if stratum.type==0
    % homogeneous stratum
    % (Only the top and bottom stratum surfaces are displayed; the interior
    % appears hollow.)
    x1_bot = x1_top;
    x1_top = x1_bot+param(stratum.thick,param_index);
    m = stratum.pmt_index;
    if ~isempty(pmt_display(m).color)
        x0 = [x1_bot,x_limit(1,2)-x_eps,x_limit(1,3)-x_eps];
        edge1 = [x1_top-x1_bot,0,0];
        if edge1(1)>2*x_eps
            x0(1) = x0(1)+x_eps;
            edge1(1) = edge1(1)-2*x_eps;
        else
            x0(1) = x0(1)+edge1(1)/4;
            edge1(1) = edge1(1)/2;
        end
        edge2 = [0,x_limit(2,2)-x_limit(1,2)+2*x_eps,0];
        edge3 = [0,0,x_limit(2,3)-x_limit(1,3)+2*x_eps];
        h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
            pmt_display(m).color,pmt_display(m).alpha);
        if ~isempty(h_)
            h{m} = h_;
        end
    end
elseif any(stratum.type==[1,2])
    % periodic stratum
    % (Homogeneous stripes are displayed as hollow tubes; structural blocks
    % are displayed as hollow boxes.)
    x1_bot = x1_top;
    x1_top = x1_bot+param(stratum.thick,param_index);
    if stratum.type==1
        % uniperiodic stratum
        % Compute stratum period (Eqs 3.24, 3.25).
        fs1 = [stratum.h11,stratum.h12]/[d21,d22;d31,d32];
        [fs21,fs31] = deal(fs1(1),fs1(2));
        den = fs21^2+fs31^2;
        [ds21,ds31] = deal(fs21/den,fs31/den);
        % Define unit basis vector [e21,e31] perpendicular to the stratum's
        % stripes, and unit basis vector [e22,e32] parallel to the stripes.
        den = sqrt(ds21^2+ds31^2);
        e21 = ds21/den;
        e31 = ds31/den;
        e22 = e31;
        e32 = -e21;
        clear fs1 fs21 fs31 den
    else
        % biperiodic stratum
        % Compute stratum period (Eq 3.20).
        ds = [d21,d22;d31,d32]/...
            [stratum.h11,stratum.h12;stratum.h21,stratum.h22];
        [ds21,ds22,ds31,ds32] = deal(ds(1,1),ds(1,2),ds(2,1),ds(2,2));
        % Define unit basis vector [e21,e31] perpendicular to the stratum's
        % stripes, and unit basis vector [e22,e32] parallel to the stripes.
        den = sqrt(ds22^2+ds32^2);
        e22 = ds22/den;
        e32 = ds32/den;
        tmp = [ds21,ds31]-[ds21,ds31]*[e22;e32]*[e22,e32];
        tmp = tmp/sqrt(tmp(1)^2+tmp(2)^2);
        e21 = tmp(1);
        e31 = tmp(2);
        clear ds den tmp
    end
    % Loop through stripes.
    c1 = param(stratum.stripe{end}.c1,param_index)-1;
    for l2 = 1:length(stratum.stripe)
        stripe = stratum.stripe{l2};
        c1_ = c1;
        c1 = param(stripe.c1,param_index);
        if stratum.type==1 || stripe.type==0
            % homogeneous stripe
            m = stripe.pmt_index;
            if ~isempty(pmt_display(m).color)
                x0 = [x1_bot,c1*[ds21,ds31]+[dx2,dx3]];
                x0(2:3) = x0(2:3)*[e21;e31]*[e21,e31]-x_max*[e22,e32];
                edge1 = [x1_top-x1_bot,0,0];
                edge2 = [0,e22,e32]*2*x_max;
                edge3 = [0,(c1_-c1)*[ds21,ds31]*[e21;e31]*[e21,e31]];
                if edge1(1)>2*x_eps
                    x0(1) = x0(1)+x_eps;
                    edge1(1) = edge1(1)-2*x_eps;
                else
                    x0(1) = x0(1)+edge1(1)/4;
                    edge1(1) = edge1(1)/2;
                end
                sqrt_ = sqrt(edge3(2)^2+edge3(3)^2);
                if sqrt_>2*x_eps
                    x0 = x0+x_eps*edge3/sqrt_;
                    edge3 = edge3-2*x_eps*edge3/sqrt_;
                else
                    x0 = x0+edge3/4;
                    edge3 = edge3/2;
                end
                x0_ = x0;
                while x0(2:3)*[e21;e31]<x_max || ...
                        (x0(2:3)+edge3(2:3))*[e21;e31]<x_max
                    h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
                        pmt_display(m).color,pmt_display(m).alpha);
                    if ~isempty(h_)
                        h{m} = h_;
                    end
                    x0(2:3) = x0(2:3)+[ds21,ds31]*[e21;e31]*[e21,e31];
                end
                x0(2:3) = x0_(2:3)-[ds21,ds31]*[e21;e31]*[e21,e31];
                while x0(2:3)*[e21;e31]>-x_max || ...
                        (x0(2:3)+edge3(2:3))*[e21;e31]>-x_max
                    h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
                        pmt_display(m).color,pmt_display(m).alpha);
                    if ~isempty(h_)
                        h{m} = h_;
                    end
                    x0(2:3) = x0(2:3)-[ds21,ds31]*[e21;e31]*[e21,e31];
                end
            end
        else
            % inhomogeneous stripe
            % Loop through structural blocks.
            c2 = param(stripe.block{end}.c2,param_index)-1;
            for l3 = 1:length(stripe.block)
                block = stripe.block{l3};
                c2_ = c2;
                c2 = param(block.c2,param_index);
                m = block.pmt_index;
                if ~isempty(pmt_display(m).color)
                    x0 = [x1_bot,c1*[ds21,ds31]+c2*[ds22,ds32]+[dx2,dx3]];
                    edge1 = [x1_top-x1_bot,0,0];
                    edge2 = [0,(c2_-c2)*[ds22,ds32]];
                    edge3 = [0,(c1_-c1)*[ds21,ds31]*[e21;e31]*[e21,e31]];
                    if edge1(1)>2*x_eps
                        x0(1) = x0(1)+x_eps;
                        edge1(1) = edge1(1)-2*x_eps;
                    else
                        x0(1) = x0(1)+edge1(1)/4;
                        edge1(1) = edge1(1)/2;
                    end
                    sqrt_ = sqrt(edge2(2)^2+edge2(3)^2);
                    if sqrt_>2*x_eps
                        x0 = x0+x_eps*edge2/sqrt_;
                        edge2 = edge2-2*x_eps*edge2/sqrt_;
                    else
                        x0 = x0+edge2/4;
                        edge2 = edge2/2;
                    end
                    sqrt_ = sqrt(edge3(2)^2+edge3(3)^2);
                    if sqrt_>2*x_eps
                        x0 = x0+x_eps*edge3/sqrt_;
                        edge3 = edge3-2*x_eps*edge3/sqrt_;
                    else
                        x0 = x0+edge3/4;
                        edge3 = edge3/2;
                    end
                    x0_ = x0;
                    while x0(2:3)*[e21;e31]<x_max || ...
                            (x0(2:3)+edge3(2:3))*[e21;e31]<x_max
                        x0__ = x0;
                        while x0(2:3)*[e22;e32]<x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]<x_max
                            h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if ~isempty(h_)
                                h{m} = h_;
                            end
                            x0 = x0+[0,ds22,ds32];
                        end
                        x0(2:3) = x0__(2:3)-[ds22,ds32];
                        while x0(2:3)*[e22;e32]>-x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]>-x_max
                            h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if ~isempty(h_)
                                h{m} = h_;
                            end
                            x0 = x0-[0,ds22,ds32];
                        end
                        x0 = x0__+[0,ds21,ds31];
                    end
                    x0(2:3) = x0_(2:3)-[ds21,ds31];
                    while x0(2:3)*[e21;e31]>-x_max || ...
                            (x0(2:3)+edge3(2:3))*[e21;e31]>-x_max
                        x0__ = x0;
                        while x0(2:3)*[e22;e32]<x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]<x_max
                            h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if ~isempty(h_)
                                h{m} = h_;
                            end
                            x0 = x0+[0,ds22,ds32];
                        end
                        x0(2:3) = x0__(2:3)-[ds22,ds32];
                        while x0(2:3)*[e22;e32]>-x_max || ...
                                (x0(2:3)+edge2(2:3))*[e22;e32]>-x_max
                            h_ = fill_box(x0,edge1,edge2,edge3,x_limit,...
                                pmt_display(m).color,pmt_display(m).alpha);
                            if ~isempty(h_)
                                h{m} = h_;
                            end
                            x0 = x0-[0,ds22,ds32];
                        end
                        x0 = x0__-[0,ds21,ds31];
                    end
                end
            end
        end
    end
elseif stratum.type==3
    % coordinate break
    dx2 = dx2+param(stratum.dx2,param_index);
    dx3 = dx3+param(stratum.dx3,param_index);
else % stratum.type==4
    for count = 1:stratum.rep_count
        for l1 = 1:length(stratum.stratum)
            [x1_top,dx2,dx3,h] = plot_stratum(...
                x1_top,dx2,dx3,h,stratum.stratum{l1},...
                param_index,pmt_display,x_limit,x_eps,x_max,...
                d21,d31,d22,d32);
        end
    end
end
end % plot_stratum


function h = fill_box(x0,edge1,edge2,edge3,x_limit,c,a)
% Display a 3-D rectangle with corner vertex x0(1,1:3) and edge vectors
% edge1(1,1:3), edge2(1,1:3), and edge3(1,1:3). The rectangle is clipped to
% limits x_limit(1,1:3) (low limit) and x_limit(2,1:3) (high limit). The
% rgb display color for faces is c(1,1:3), and for edges is either c(1,1:3)
% if isrow(c), or c(2,1:3) otherwise. The opacity (alpha) for faces is a,
% and for edges is either a if isscalar(a), or a(2) otherwise. Return a
% handle to one of the faces, or [] if all faces are entirely outside the
% limits.
if isrow(c)
    c(2,:) = 0; % c(1,:): faces, c(2,:): edges
end
if isscalar(a)
    a = [a;a]; % a(1): faces, a(2): edges
end
h = [];
h_ = fill_face([...
    x0;...
    x0+edge1;...
    x0+edge1+edge2;...
    x0+edge2...
    ],x_limit,c,a);
if ~isempty(h_)
    h = h_;
end
h_ = fill_face([...
    x0+edge3;...
    x0+edge3+edge1;...
    x0+edge3+edge1+edge2;...
    x0+edge3+edge2...
    ],x_limit,c,a);
if ~isempty(h_)
    h = h_;
end
h_ = fill_face([...
    x0;...
    x0+edge2;...
    x0+edge2+edge3;...
    x0+edge3...
    ],x_limit,c,a);
if ~isempty(h_)
    h = h_;
end
h_ = fill_face([...
    x0+edge1;...
    x0+edge1+edge2;...
    x0+edge1+edge2+edge3;...
    x0+edge1+edge3...
    ],x_limit,c,a);
if ~isempty(h_)
    h = h_;
end
h_ = fill_face([...
    x0;...
    x0+edge3;...
    x0+edge3+edge1;...
    x0+edge1...
    ],x_limit,c,a);
if ~isempty(h_)
    h = h_;
end
h_ = fill_face([...
    x0+edge2;...
    x0+edge2+edge3;...
    x0+edge2+edge3+edge1;...
    x0+edge2+edge1...
    ],x_limit,c,a);
if ~isempty(h_)
    h = h_;
end
end % fill_box


function h = fill_face(x,x_limit,c,a)
% Display a rectangle face with corner vertices x(:,1:3).
h = [];
x = clip_x(x,x_limit(1,1));
if isempty(x)
    return
end
x = -clip_x(-x,-x_limit(2,1));
if isempty(x)
    return
end
x = x(:,[2,3,1]);
x = clip_x(x,x_limit(1,2));
if isempty(x)
    return
end
x = -clip_x(-x,-x_limit(2,2));
if isempty(x)
    return
end
x = x(:,[2,3,1]);
x = clip_x(x,x_limit(1,3));
if isempty(x)
    return
end
x = -clip_x(-x,-x_limit(2,3));
if isempty(x)
    return
end
x = x(:,[2,3,1]);
h = patch(x(:,2),x(:,3),x(:,1),c(1,:),'EdgeColor',c(2,:), ...
    'FaceAlpha',a(1),'EdgeAlpha',a(2));
end % fill_face


function x = clip_x(x,x1_limit)
% Clip rectance face x(:,1:3) to limit x(:,1)>=x1_limit.
if all(x(:,1)>=x1_limit)
    return
end
if all(x(:,1)<=x1_limit)
    x = [];
    return
end
nx = size(x,1);
x_ = [];
xa = x(end,:);
xb = x(1,:);
for j_ = 1:nx
    xc = x(mod(j_,nx)+1,:);
    if xb(1)>=x1_limit
        if xa(1)<x1_limit
            x_(end+1,:) = ...
                xa+(xb-xa)*((x1_limit-xa(1))/(xb(1)-xa(1))); %#ok<AGROW>
        end
        x_(end+1,:) = xb; %#ok<AGROW>
        if xc(1)<x1_limit
            x_(end+1,:) = ...
                xc+(xb-xc)*((x1_limit-xc(1))/(xb(1)-xc(1))); %#ok<AGROW>
        end
    end
    xa = xb;
    xb = xc;
end
x = x_;
end % clip_x
