% demo_cartesian_tse_recon.m
% Written by Nam Gyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 05/21/2021, Last modified: 05/21/2021

%% Notation
%--------------------------------------------------------------------------
% We denote X, Y, Z as the horizontal, vertical, and through-magnet axes,
% respectively in the physical coordinate system, and the corresponding
% coordinates as x, y, and z. We also denote U, V, and W as the readout,
% phase encoding, and slice axes, respectively in the logical coordinate
% system, and the corresponding coordinates as u, v, and w.
%
%                     Anterior
%                        ^ (+y)
%                        |
%           Right        |            Left
%                        |
%                  Head  +--------> (+x)
%                       /
%                      /
%                     /
%               Foot v (+z)
%                      Posterior
%
%--------------------------------------------------------------------------

%% Clean slate
close all; clear; clc;

%% Set directory names
computer_type = computer;
if strcmp(computer_type, 'PCWIN64')
    src_directory = 'F:\USC\MREL\LowField\LungImaging\T2measurement\mrdCode\src';
    ismrmrd_directory = 'F:\USC\MREL\Tool\ISMRMD';
    thirdparty_directory = 'F:\USC\MREL\LowField\LungImaging\T2measurement\mrdCode\thirdparty';
    data_parent_directory = 'F:\USC\MREL\LowField\LungImaging\T2measurement\Data\0524VOL11Lung\RAW';
elseif strcmp(computer_type, 'GLNXA64')
    src_directory = '';
    ismrmrd_directory = '/server/home/nlee/ismrmrd';
    data_parent_directory = '';
end

%% Add paths
addpath(genpath(src_directory));
addpath(genpath(thirdparty_directory));
addpath(genpath(ismrmrd_directory));

%% Define input filenames
%--------------------------------------------------------------------------
% 20200520_apple, axial
%--------------------------------------------------------------------------
data_directory = fullfile(data_parent_directory);
dir_info = dir(fullfile(data_directory, '*BL*.dat')); % Change to whatever pattern you need.
tse_filenames = cell(1,length(dir_info));
for idx = 1:length(tse_filenames)
    dat_filename = dir_info(idx).name;
    dot_loc = strfind(dat_filename, '.');
    tse_filenames{idx} = dat_filename(1:dot_loc-1);
end
%% Define constants
gamma = 2.67522212e+8; % gyromagnetic ratio for 1H [rad/sec/T]

%% Define a fullpath to each filename
nr_echoshifts = length(tse_filenames);
tse_noise_fullpaths = cell(nr_echoshifts,1);
tse_data_fullpaths  = cell(nr_echoshifts,1);
tse_dat_fullpaths   = cell(nr_echoshifts,1);

for idx = 1:nr_echoshifts
    tse_noise_fullpaths{idx} = fullfile(data_directory,'noise', sprintf('noise_%s.h5', tse_filenames{idx}));
    tse_data_fullpaths{idx}  = fullfile(data_directory,'h5', sprintf('%s.h5', tse_filenames{idx}));
    tse_dat_fullpaths{idx}   = fullfile(data_directory, sprintf('%s.dat', tse_filenames{idx}));
end

%% Set flags for Recon
flag.csm_method = 'sos'; % coil estimation methods: sos, walsh
flag.cal_shape = [22 22]; % calibration region in k-space for coil estimation
flag.ignoreSeg = 'false'; % ignore segments (Turbo factos)
flag.TF = 7;

%% Get CSM from the dataset without echo shift
scan_nr = 1; % scan number without echo shift
[im,header,x,y,z,imc,csm] = cartesian_tse_recon(tse_noise_fullpaths{scan_nr}, tse_data_fullpaths{scan_nr}, tse_dat_fullpaths{scan_nr}, flag);
[N1,N2,N3] = size(im);
Nc = size(csm, ndims(csm));
csm = reshape(csm, [N1 N2 N3 Nc]);

%% Reconstruct TSE Cartesian images
im_tse = complex(zeros(N1, N2, N3, nr_echoshifts, 'double'));
for idx = 1:nr_echoshifts
    [im,header,x,y,z,imc,~] = cartesian_tse_recon(tse_noise_fullpaths{idx}, tse_data_fullpaths{idx}, tse_dat_fullpaths{idx}, flag);
    imc = reshape(imc, [N1 N2 N3 Nc]);
    im_tse(:,:,:,idx) = sum(conj(csm) .* imc, 4);
end

%% Display all images
figure, montage(real(im_tse), 'DisplayRange', [], 'Size', [2 7]);
colormap(gray(256));
colorbar;


figure, montage(angle(im_tse)*180/pi, 'DisplayRange', [], 'Size', [2 7]);
colormap(hsv(256));
colorbar;

return

%% ----Generate the ROI mask----
temp = im_tse(:,:,:,7);
temp = abs(im_tse(:,:,:,7))./mean(abs(temp(:)));
temp = temp>1.2;
[L,num]=bwlabel(temp,4);
x=zeros(1,num);
for ii=1:num
    x(ii)=sum(sum(L==ii));
end
[m,ind]=max(x);
bw_img = zeros(size(temp));
% for n = 1:length(x)
%     bw_img = (L==n) | bw_img;
% end
bw_img = (L==ind); % Mask [nx ny nz]
se = strel('line',3,3);
for n_erosion = 1:1
    bw_img = imerode(bw_img,se);
end
bd_bw = bwboundaries(bw_img);

figure,
imagesc(abs(im_tse(:,:,:,7)));colormap(gray)
hold on;
for k = 1:length(bd_bw)
    boundary = bd_bw{k};
    plot(boundary(:,2), boundary(:,1), 'r--', 'LineWidth', 0.5)
end
axis off image

%% Display a signal from an individual voxel
echoshifts = [-15 -10 -7 -5 -3 -1 0 1 3 5 7 10 15 17];
time_range =  header.sequenceParameters.TE + echoshifts;
time_idx = 7:length(echoshifts);
% time_idx = 1:7;
xVal = time_range(time_idx);
figure;
subplot(1,2,1); imagesc(abs(im_tse(:,:,:,7))); axis image;
colormap(gray(256));
colorbar;
for idx = 1:1000
    [x,y,button] = ginput(1);
    row = floor(y)
    col = floor(x)
    subplot(1,2,2);
    plot(xVal.', reshape(real(im_tse(row,col,:,time_idx)), [length(time_idx) 1]));
    if button == 3
        return;
    end
end



%
% figure, montage(abs(im_tse2), 'DisplayRange', [], 'Size', [2 7]);
% colormap(gray(256));
% colorbar;
%
%
% figure, montage(angle(im_tse2)*180/pi, 'DisplayRange', [], 'Size', [2 7]);
% colormap(hsv(256));
% colorbar;
%
%
% %%
% figure('Color', 'w');
% montage(abs(csm), 'DisplayRange', []);
% colormap(parula(256));
% colorbar;
%
% figure('Color', 'w');
% montage(angle(csm)*180/pi, 'DisplayRange', []);
% colormap(hsv(256));
% colorbar;
%
% %%
% figure;
% imagesc(abs(im)); axis image;
% colormap(gray(256));
% colorbar;
%
% figure;
% imagesc(angle(im)*180/pi); axis image;
% colormap(hsv(256));
% colorbar;

%% ----Exp fitting----
% temp  = cat(3,img_List{1:length(img_List)});
% temp = temp(:,75:end,:);
% temp = squeeze(num2cell(temp,[1 2]));
% temp = squeeze(num2cell(c_img_List{nfk},[1 2]));
% flag.fitting = 'c + a * exp(-x/b)'; % T2decay = 'c + PD * exp(-x/T2)'; exp1
flag.fitting = 'exp1'; % T2decay = 'c + PD * exp(-x/T2)'; Exp1
flag.fitting = 'a*exp(b*x) + c*exp(d*x)'; % T2decay = 'c + PD * exp(-x/T2)'; Exp1


temp = squeeze(im_tse(:,:,:,time_idx));
for idx = 1: size(temp,3)
    echoshifts_signal = temp(:,:,idx); %each columne for one TE; each row for one voxel in ROI
    signal_TE(:,idx) = echoshifts_signal(bw_img);
end
tic
parfor k = 1:length(signal_TE)
    k
    f = fit(xVal.'/1000,real(signal_TE(k,1:end)).',flag.fitting);%
    R2(k) = -f.b;
    A(k) = f.a;
%     c(k) = f.c;
end
toc
clear signal_TE

%%  ----Display T2* maps----
% for nfk = 1:length(f_k)
T2Starmap = zeros(size(bw_img));
T2Star = 1./R2;
T2Starmap(bw_img) = T2Star;
T2Starmap(T2Starmap<0) = 0;
T2Starmap(isinf(T2Starmap)|isnan(T2Starmap)) = 0;
% end
% figure,imagesc([T2map.*(bw_img),...
%     reshape(cat(3,T2map_w{1:length(f_k)}),[Nx,Ny*length(f_k)]).*repmat(bw_img,[1,length(f_k)])]*1000,[0 50]);

figure,
imagesc(T2Starmap*1000,[0 50])
h = colorbar;
ylabel(h, 'T2*(ms)','FontWeight','bold')
colormap(hot)
axis image off

%%  ----Display T2up maps----
T2Upmap = zeros(size(bw_img));
T2Up = 1./R2;
T2Upmap(bw_img) = T2Up;
T2Upmap(isinf(T2Upmap)|isnan(T2Upmap)) = 0;
T2Upmap(T2Upmap>0) = 0;

figure,
imagesc(T2Upmap*1000,[-50 0])
h = colorbar;
ylabel(h, '1/(R2-R2'') (ms)','FontWeight','bold')
colormap(hot)
axis image off

%% T2'
T2Prime = 1./ ( (1./T2Star - 1./T2Upmap) /2);
T2Primemap = zeros(size(bw_img));;
T2Prime = 1./((1./T2Starmap(bw_img) - 1./T2Upmap(bw_img))/2);
T2Primemap(bw_img) = T2Prime;
T2Primemap(T2Primemap<0) = 0;
T2Primemap(isinf(T2Primemap)|isnan(T2Primemap)) = 0;
T2primemap_ps = T2Primemap;

figure,
imagesc(T2Primemap*1000,[0 50])
h = colorbar;
ylabel(h, 'T2'' (ms)','FontWeight','bold')
colormap(hot)
axis image off

%% T2
T2 =  1./((1./T2Star_RS15_SoS22_Exp1 + 1./T2Up_RS15_SoS22_Exp1)/2);
T2map = zeros(size(bw_img));
T2map(bw_img) = T2;
T2map(T2map<0) = 0;
T2map(isinf(T2map)|isnan(T2map)) = 0;

figure,
imagesc(T2map*1000,[0 1000])
h = colorbar;
ylabel(h, 'T2 (ms)','FontWeight','bold')
colormap(hot)
axis image off
%%
figure;
surf(x*1e3, y*1e3, z*1e3, abs(im), 'EdgeColor', 'none');
colormap(gray(256));
colorbar;
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
title('Siemens Physical Coordinate System (Device Coordinate System)');

figure;
surf(x*1e3, y*1e3, z*1e3, angle(im)*180/pi, 'EdgeColor', 'none');
colormap(hsv(256));
colorbar;
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
title('Siemens Physical Coordinate System (Device Coordinate System)');

