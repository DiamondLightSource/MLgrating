function [w,n,comment] = read_nk(filename,scale,filepath)
% [w,n,comment] = read_nk(filename,scale,filepath)
%
% Read IMD-format nk file. Each text line in the file can be one of the
% following:
%   (1) an empty string
%   (2) a comment line beginning with ';'
%   (3) a line of the form
%         w  n  k
%       where
%         w = wavelength (Angstroms)
%         n = real part of refractive index
%         k = imaginary part of refractive index (k >= 0)
%
% Syntax:
%   [w,n,comment] = read_nk(filename,scale,filepath)
%
% Inputs:
%
%   filename: file name; string
%
%   scale: wavelength scaling factor; real (optional; default = 1)
%
%   filepath: file path prefix; string (optional; default = './'; filepath
%   must be '/' terminated.)
%
% Outputs:
%
%   w: wavelengths in units of A/scale (e.g., set scale = 0.1 for nm,
%   1e-4 for micron), column vector
%
%   n: complex refractive index (imag(n) >= 0), complex, size-matched to w
%
%   comment: file comment lines, cell row vector of string
%
% Notes:
%
%   A directory of nk files for a variety of materials can be freely
%   downloaded from http://www.rxollc.com/idl/. Download IMD, unpack with
%   the MATLAB function untar, and extract the nk directory.
%
%   The MATLAB function interp1 (linear interpolation) can be used to
%   change the refractive index's wavelength sampling, e.g., from [w,n] to
%   [w_,n_]:
%     n_ = interp1(w,n,w_);
%
% Version 04-Jun-2022
% Author: Kenneth C. Johnson, KJ Innovation https://kjinnovation.com/

if nargin<3
    filepath = './';
    if nargin<2
        scale = 1;
    end
end
fid = fopen([filepath filename]);
if fid==-1
    error(['Could not open file ' filepath filename]);
end
comment = {};
wnk = zeros(3,0);
next_line = fgetl(fid);
while ~isequal(next_line,-1)
    if ~isempty(next_line)
        if isequal(next_line(1),';')
            comment{end+1} = next_line; %#ok<AGROW>
        else
            wnk(:,end+1) = sscanf(next_line,'%g'); %#ok<AGROW>
        end
    end
    next_line = fgetl(fid);
end
w = scale*wnk(1,:).';
n = (wnk(2,:)+1i*wnk(3,:)).';
f = find(diff(w)==0);
if ~isempty(f)
    w(f) = [];
    n(f+1) = (n(f)+n(f+1))/2;
    n(f) = [];
end
fclose(fid);
