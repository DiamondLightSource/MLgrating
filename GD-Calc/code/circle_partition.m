function x = circle_partition(N)
% x = circle_partition(N)
%
% Partition the first quadrant of the unit circle into N rectangular blocks
% with corner vertices at [x(1),x(end)], [x(2),x(end-1)], ...
% [x(end),x(1)]. x(j) is monotonic decreasing with j.
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

% Documentation reference:
%   Grating Diffraction Calculator (GD-Calc)
%   Demo and Tutorial Guide
%   (GD-Calc_Demo.pdf, version 04-Jun-2022)
%   Appendix A
% All equation ("Eq") references below are from GD-Calc_Demo.pdf.

theta1 = [ ...
    (pi/4)*((N-1/2)*sqrt(3)+1)^(-2/3), ... % estimate from Eq A.14
    (pi/4)*((N+1/2)*sqrt(3)+1)^(-2/3)]; % second point for linear fit
err = [ ...
    target_err(theta1(1),N), ...
    target_err(theta1(2),N)]; % Eq A.4 error
% Find zero crossing of error via linear interpolation:
dtheta1 = diff(theta1);
theta1 = theta1(1)-err(1)/diff(err)*dtheta1;
% Iteratively refine zero crossing:
prev_maxerr = max(abs(err));
prev_derr = diff(err);
while true
    if sign(err(1))~=sign(err(2))
        dtheta1 = dtheta1/10;
    end
    if abs(dtheta1)<=eps*abs(theta1)
        break
    end
    theta1 = theta1+dtheta1*[-0.5,0.5];
    err = [target_err(theta1(1),N),target_err(theta1(2),N)];
    maxerr = max(abs(err));
    if maxerr>=prev_maxerr
        break
    end
    prev_maxerr = maxerr;
    derr = diff(err);
    if abs(derr)>=abs(prev_derr) || abs(derr)<=eps*maxerr
        break
    end
    prev_derr = derr;
    theta1 = theta1(1)-err(1)/derr*dtheta1;
end
theta1 = theta1(1);
area = theta1-0.5*sin(2*theta1); % Eq A.2
theta = [theta1,zeros(1,N-1)];
outside_pt = true;
for j = 2:N
    % Eq A.3
    theta(j) = next(theta(j-1),area,outside_pt);
    outside_pt=~outside_pt;
end
theta = [theta,pi/2-theta(end:-1:1)]; % Eq A.1
x = cos(theta(1:2:end));
end % circle_partition

function theta = next(theta,area,outside_pt)
% Calculate next theta from Eq A.3.
theta_ = theta+2*sqrt(area/sin(2*theta)); % Eq A.5 (estimate)
prev_dtheta_ = [];
c = cos(theta);
s = sin(theta);
while true
    c_ = cos(theta_);
    s_ = sin(theta_);
    % Correct the theta_ estimate via differential approximation to Eq A.3.
    if outside_pt
        % "j odd" branch in Eq A.3
        theta_ = theta_-(c*s_-0.5*(c*s+c_*s_)-0.5*(theta_-theta)-area)/...
            (c*c_-0.5*(c_*c_-s_*s_)-0.5);
    else
        % "j even" branch in Eq A.3
        theta_ = theta_-(c_*s-0.5*(c*s+c_*s_)+0.5*(theta_-theta)-area)/...
            (-s_*s-0.5*(c_*c_-s_*s_)+0.5);
    end
    if isempty(prev_dtheta_)
        prev_dtheta_ = theta_-theta;
    else
        dtheta_ = theta_-prev_theta_;
        if abs(dtheta_)<eps || abs(dtheta_)>=abs(prev_dtheta_)
            break
        end
        prev_dtheta_ = dtheta_;
    end
    prev_theta_ = theta_;
end
theta = theta_;
end % next

function err = target_err(theta1,N)
% Calculate error in Eq A.4
area = theta1-0.5*sin(2*theta1); % Eq A.2
theta = theta1;
outside_pt = true;
for j = 2:N
    % Calculate next theta from Eq A.3.
    theta = next(theta,area,outside_pt);
    outside_pt=~outside_pt;
end
theta_ = next(theta,area,outside_pt);
err = theta+theta_-pi/2;
end % target_err
