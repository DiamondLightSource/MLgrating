function MLgrating_plot(grating,materials)

% find number of materials defined
num_mat=numel(materials);

% set up grating colours based on how many materials there are
if isequal(num_mat,3)
    colours=[1 1 0;1 0 0;0 1 0];
elseif isequal(num_mat,2)
    colours=[1 1 0;1 0 0;0 0 0];
elseif isequal(num_mat,1)
    colours=[1 1 0;0 0 0;0 0 0];
end

% initialise pmt_display parameters (MLgrating_def_lam currently ensures 
% that number of materials is always 4, but for more general cases code 
% will need to be modified)
for i=1:4
    pmt_display(i).name = 'Vacuum';
    pmt_display(i).color = [1 1 1];
    pmt_display(i).alpha = 1;
end

for i=1:num_mat
    pmt_display(i+1).name=materials{i};
    pmt_display(i+1).color=colours(i,:);
    pmt_display(i+1).alpha = 1;
end

% find number of strata in grating structure
num_strata=numel(grating.stratum);

% calculate total thickness of grating (sum all thickness of strata)
H=0;
for i=1:num_strata
    H=H+grating.stratum{i}.thick;
end

% extract grating D (d-spacing) from grating structure
D=grating.d21;

% define limits on 3D plot of grating
x_limit = [-0.5*H,-D/2,-D/2;1*H,D+D/2,D/2];

% use GD-Calc intrinsic plotting
gdc_plot(grating,1,pmt_display,x_limit);

% tweak view to hide strange plotting anomaly
view(0,2)

% automatically adjust scales of axes to improve visability
axis normal

% make background of figure white
set(gcf,'color','w')

% draw 3D box around image of grating
box on
set(gca,'BoxStyle','full');