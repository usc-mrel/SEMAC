clear

%% Add path
addpath(genpath('./functions'));
addpath(genpath('/server/home/bli/MetalImaging/Code/mrdCode/'))
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    ismrmrd_directory = 'F:\USC\MREL\Tool\ISMRMD';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
    data_parent_directory = '';
elseif strcmp(computer_type,'MACI64')
    ismrmrd_directory = '/Users/bli/ismrmrd';
end
addpath(genpath(ismrmrd_directory));

%% Loading data
data_folder = '/server/home/bli/Data/SEMAC/SEMACforSiemens/SEMACTest_20220526/raw/SEMAC/';
h5_SEMAC = [data_folder,'h5/meas_MID01144_FID33261_SEMAC12_Default.h5'];
noise_SEMAC =  [data_folder,'noise/noise_meas_MID01144_FID33261_SEMAC12_Default.h5'];

if ~exist(h5_SEMAC) || ~exist(noise_SEMAC)
        demo_convert_twix_to_ismrmrd([data_folder]);
end

Read_flags.h5_fileList       = h5_SEMAC;
Read_flags.noise_fileList    = noise_SEMAC;
Read_flags.RemoveOS          = true; % remove oversampling
Read_flags.IgnoreSeg         = true; % concatanate segmemtns
Read_flags.DoAverage         = true; % Do averages (if 'average' was set during data acquistion)
Read_flags.CropPhaseEncoding = true;
Read_flags.Squeeze           = true;
Read_flags.os                = 2; % oversampling rate (Siemens Default value, don't change it)
Read_flags.noisePreWhitening = true;

[c_img, kdata, noise, info] = readSEMAC_ismrmd(Read_flags); % kdata:[Nkx,Nky,Nkz,NSlice,NCoil]

[Nx, Ny,Nz, Nslice,Ncoil] = size(kdata{1});

%% Apply the noise prewhitening matrix on k-space before Recon
% ---- Calculate noise covariance ----
for k = 1:length(noise)
    [Psi, inv_L] = calculate_noise_covariance(noise{k});
    kdata_prew{k} = prewhitening_SEMAC(kdata{k},inv_L);
end

%% Re-order slice 
kdata_prew_coils_slice= zeros(size(kdata_prew{1})); 
kdata_prew_coils_slice(:,:,:,2:2:end,:) = kdata_prew{1}(:,:,:,ceil(Nslice/2+1):end,:);
kdata_prew_coils_slice(:,:,:,1:2:end,:) = kdata_prew{1}(:,:,:,1:ceil(Nslice/2+1)-1,:);


%% Recon SEMAC
SEMAC_coils = fft3c(kdata_prew_coils_slice);

%% CSM estimation
csm_option.method = 'walsh'; % coil estimation method: 'walsh', 'sos'
csm_option.cal_shape = [16 16]; %calibration region
kdata_prew_coils_slice_z = 1/sqrt(Nz)*fftshift(fft(ifftshift(kdata_prew_coils_slice,3),[],3),3);
csm_option.kdata = squeeze(kdata_prew_coils_slice(:,:,7,:,:)); % Please change this number to the kz of 0. (7 for SEMAC12; 4 for SEMAC6)
[csm, cal_im] = coil_estimation(csm_option);

csm = repmat(csm,[1 1 1 1 Nz]);
csm = permute(csm,[1 2 5 3 4]);

SEMAC_CC_Walsh_prewhit =  sum( conj(csm) .* fft2c(kdata_prew_coils_slice_z), 5) ; % Coil combined images ./ sqrt(sum(abs(csm).^2,4))

%% SEMAC combination (SEMAC=12 ; 21 slices)
SEMAC_slice_SOS = zeros(size(SEMAC_CC_Walsh_prewhit,1,2,4));
for ss = 1:21
    if (ss-5 >0) && (ss+ 6<= 21) 
    SEMAC_slice_SOS(:,:,ss) = sqrt(abs(SEMAC_CC_Walsh_prewhit(:,:,12,ss-5)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,11,ss-4)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,10,ss-3)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,9,ss-2)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,8,ss-1)).^2 +...
        abs(SEMAC_CC_Walsh_prewhit(:,:,7,ss)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,6,ss+1)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,5,ss+2)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,4,ss+3)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,3,ss+4)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,2,ss+5)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,1,ss+6)).^2);
    end
end

%% SEMAC combination (SEMAC=6 ; 21 slices)
SEMAC_slice_SOS = zeros(size(SEMAC_CC_Walsh_prewhit,1,2,4));
for ss = 1:21
    if (ss-2 >0) & (ss+3 <= 21) 
    SEMAC_slice_SOS(:,:,ss) = sqrt(abs(SEMAC_CC_Walsh_prewhit(:,:,6,ss-2)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,5,ss-1)).^2 +...
        abs(SEMAC_CC_Walsh_prewhit(:,:,4,ss)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,3,ss+1)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,2,ss+2)).^2+...
        abs(SEMAC_CC_Walsh_prewhit(:,:,1,ss+3)).^2);
    end
end

% %% SEMAC combination for each coil (SEMAC = 6; all slices)
% img_coils = fft2c(kdata_prew_coils_slice_z);
% SEMAC_slice_SOS_coils = zeros(size(img_coils,1,2,4,5));
% for ss = 1:21
%     if (ss-2 >0) & (ss+3 <= 21) 
%     SEMAC_slice_SOS_coils(:,:,ss,:) = sqrt(abs(img_coils(:,:,6,ss-2,:)).^2+...
%         abs(img_coils(:,:,5,ss-1,:)).^2 +...
%         abs(img_coils(:,:,4,ss,:)).^2+...
%         abs(img_coils(:,:,3,ss+1,:)).^2+...
%         abs(img_coils(:,:,2,ss+2,:)).^2+...
%         abs(img_coils(:,:,1,ss+3,:)).^2);
%     end
% end

%% Save resutls
fprintf('Saving SEMAC images in %s\n',[data_folder,'img_semac12'])
save([data_folder,'/img_semac12'],'SEMAC_slice_SOS','-v7.3');
save([data_folder,'/spectral_semac6'],'SEMAC_CC_Walsh_prewhit','-v7.3');

