classdef pgm < photon
    properties
        order  	% Diffraction order (normally 1)
        N       % Line density of grating (lines/mm)
        cff   	% cos(beta)/cos(alpha) (see definitions below)
    end
    properties (Dependent)
        theta       % incident/reflection angle of plane mirror in deg
                    % to keep pgm outgoing beam horizontal (relative to normal)
        alpha       % Grating incidence angle in deg (relative to normal)
        beta        % Grating diffraction angle in deg (relative to normal)
    end
    properties (Dependent, Hidden)
        B           % parameter (in metres) which enables easier 
                    % calculation of theta, alpha and beta
    end
   methods
       function value=get.B(pgm)
           lambda_q=pgm.order.*pgm.wavelength.*1e-9.*pgm.N.*1e3;
           value=lambda_q.*pgm.cff./(1-pgm.cff.^2);
       end
       function value=get.alpha(pgm)
           value=asind(pgm.B./pgm.cff+sqrt(pgm.B.^2+1));
       end
       function value=get.beta(pgm)
           value=asind(pgm.B.*pgm.cff+sqrt(pgm.B.^2+1));
       end
       function value=get.theta(pgm)
           value=(pgm.alpha+pgm.beta)./2;
       end
       function obj = pgm(val)
           obj@photon(val);
           if exist('val','var')
               for i=7:9
                   eval(['if isfield(val,''' obj.inputpropnames{i} ''') || isprop(val,''' obj.inputpropnames{i} ''')' newline ...
                       'obj.string=repmat(' '''1,''' ',1,i-1);' newline ...
                       'end'])
                   eval(['if isfield(val,''' obj.inputpropnames{i} ''') || isprop(val,''' obj.inputpropnames{i} ''')' newline ...
                       'obj.' obj.inputpropnames{i} '=zeros(' obj.string 'numel(val.' obj.inputpropnames{i} '));' newline ...
                       'obj.' obj.inputpropnames{i} '(' obj.string ':)=val.' obj.inputpropnames{i} ';' newline ...
                       'end'])
               end
           end
       end
   end
end