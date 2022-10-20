function [Noise] = process_noise_ismrmrd_data(noise_fullpath)


%% Read an ismrmrd file
tic; fprintf('Reading an ismrmrd file: %s... ', noise_fullpath);
if exist(noise_fullpath, 'file')
    dset = ismrmrd.Dataset(noise_fullpath, 'dataset');
    fprintf('done! (%6.4f sec)\n', toc);
else
    error('File %s does not exist.  Please generate it.' , noise_fullpath);
end
raw_data = dset.readAcquisition();

%% Sort noise data
is_noise = raw_data.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
meas = raw_data.select(find(is_noise));
nr_repetitions = length(meas.data); % number of repetitions
[nr_samples,nr_channels] = size(meas.data{1});
Noise = complex(zeros(nr_channels, nr_samples, nr_repetitions, 'double'));
for idx = 1:nr_repetitions
    Noise(:,:,idx) = meas.data{idx}.'; % nr_samples x nr_channels => nr_channels x nr_samples
end
end
