function [P,LR,LR_k,HR_k] = PSRecon(img, f_k)
%   Phase-Sensitivity Reconstruction
%   Input:
%       kdata: k-space data
%       c_img: image for each coil
%       f_k: truncation factor in k-space(kx,ky,kz) (0~1)
%   Output:
%       P: phase map
%       LR: low-resolution image\

%   Author:
%       Bochao Li
%       E-mail: bochaoli@usc.edu
%--------------------------------------------------------------------------
dim = size(img);
LR_k = zeros(dim);
HR_k = fft2c(img);


k_rang = dim .* f_k;
x_range = [-floor(k_rang(1)/2),ceil(k_rang(1)/2)-1]+floor(dim(1)/2)+1;
y_range = [-floor(k_rang(2)/2),ceil(k_rang(2)/2)-1]+floor(dim(2)/2)+1;
LR_k(x_range(1):x_range(2),y_range(1):y_range(2)) = HR_k(x_range(1):x_range(2),y_range(1):y_range(2));

LR = ifft2c(LR_k);

% Extract phase
P = angle(LR);

end

