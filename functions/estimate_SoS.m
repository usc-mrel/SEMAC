function csm = estimate_SoS(cal_im,option)

%   Estimates relative coil sensitivity maps from a set of coil images
%   using Sum-of-Square method
%   INPUT:
%
%     - kspace     Nx x Ny x Ncoil
%
%   OUTPUT:
%
%     - csm    Relative coil sensitivity maps (Nx x Ny x Ncoil)

%----------------------------------------------------------------------
% Estimate coil sensitivity maps (N1 x N2 x Nc)
%----------------------------------------------------------------------
csm = bsxfun(@rdivide, cal_im, sqrt(sum(abs(cal_im).^2, 3)));
end

