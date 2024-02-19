function inc_field = MLgrating_def_inc_field(PGM)

% set inc_field wavelength to that defined in instance of pgm class
inc_field.wavelength = PGM.wavelength;

% theta3 = angle between the grating lines (x3 axis) and the incident wave
% vector
theta3 = 90;

% f2 and f3 are the grating-tangential (x2 and x3) coordinate projections
% of the incident field's spatial-frequency vector, Eq's 4.1-2 in
% GD-Calc.pdf and Eq 1 in GD-Calc_Demo.pdf. The grating-normal (x1)
% projection is implicitly
%   f1 = -sind(theta3)*cosd(PGM.alpha)/wavelength
inc_field.f2 = sind(theta3)*sind(PGM.alpha)/PGM.wavelength;
inc_field.f3 = cosd(theta3)/PGM.wavelength;