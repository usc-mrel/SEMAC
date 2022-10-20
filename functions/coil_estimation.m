function [csm, cal_im] = coil_estimation(option)
%   The function estimates coil sensitivity maps

%   Input:
%       - option:
%           (1) method: coil estimatino method (walsh, sos...)
%           (2) cal_shape: Calibration shape in kspac (Nx x Ny)
%           (3) kdata: k-space (Nx x Ny x Ncoil x Nslice)

%   Output:
%       - csm: coil-senstivity map (Nx x Ny x Ncoil x Nslice)
%       - cal_img: calibration image (Nx x Ny x Ncoil x Nslice)

%   Author: Bochao Li
%   Email: bochaoli@usc.edu

method = option.method; % coil estimation method: walsh, sos
cal_shape = option.cal_shape; % calibration region
kdata = option.kdata;
[Nx, Ny,Nslice, Ncoil] = size(kdata);

cal_data = crop(kdata,[cal_shape Nslice Ncoil]);
cal_data = bsxfun(@times, cal_data, hamming(cal_shape(1)) * hamming(cal_shape(2)).');
cal_im = fft2c(zpad(cal_data, [Nx Ny Nslice Ncoil]));

for  nslice = 1 : size(cal_im,3)
    if strcmp(method, 'walsh')
        csm(:,:,:,nslice) = ismrm_estimate_csm_walsh(squeeze(cal_im(:,:,nslice,:))); % Walsh  (Nx x Ny x Ncoil x Nslice)
    elseif strcmp(method, 'sos')
        csm(:,:,:,nslice) = estimate_SoS(squeeze(cal_im(:,:,nslice,:))); % SoS (Nx x Ny x Ncoil x Nslice)
    end
end

csm = permute(csm,[1 2 4 3]); % Nx x Ny x Nslice x Ncoil


