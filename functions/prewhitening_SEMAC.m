function [kdata_prewhitening] = prewhitening_SEMAC(kdata,inv_L)
%   Calculae noise covariance marix from noise samples
%   Input:
%       - kdata: k-space data (Nx x Ny x Nslice x Ncoil)

%   Output:
%       - Psi: scaled noise coavariance matrix (Ncoil x Ncoil)
%       - inv_L: inverse sqaure roof of the noise covariance (Ncoil x Ncoil)

%   Author: Bochao Li
%   Email: bochaoli@usc.edu

[Nx, Ny, Nz, Nslice, Ncoil] = size(kdata);
kdata_prewhitening = ipermute(reshape(inv_L * reshape(permute(kdata,[5 1 2 3 4]), [Ncoil Nx*Ny*Nz*Nslice]), [Ncoil Nx Ny Nz Nslice]), [5 1 2 3 4]);

end

