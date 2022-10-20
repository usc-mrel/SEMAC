function [csm, cal_im] = estimate_SoS(cal_im)

%   Estimates relative coil sensitivity maps from a set of coil images
%   using Sum-of-Square method
%   INPUT:
%
%     - cal_im    [x,y,coil]  calibration image 
%
%   OUTPUT:
%
%     - csm     [x,y,coil    : Relative coil sensitivity maps]


%----------------------------------------------------------------------
% Calculate the calibration region of k-space
%----------------------------------------------------------------------

%----------------------------------------------------------------------
% Estimate coil sensitivity maps (N1 x N2 x Nc)
%----------------------------------------------------------------------
csm = bsxfun(@rdivide, cal_im, sqrt(sum(abs(cal_im).^2, 3)));
end

