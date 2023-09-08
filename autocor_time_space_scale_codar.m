% autocor_time_space_scale_codar

load('/Volumes/home/jmv208/SWARM_data/CODAR_filled_zeroed.mat');

Domain = 'Palmer';
% Palmer ROI
latlim = [-65.5 -64.5];
lonlim = [-66.5 -63.0];

% READ MATFILE file
times = CODAR.dnum; % in days
numFrames = numel(times);

palmer_lon_rho = CODAR.lon; % decimal degree
palmer_lat_rho = CODAR.lat; % decimal degree
palmer_u_east_surface = CODAR.u; % in cm/s
palmer_v_north_surface = CODAR.v; % in cm/s
domain = find((swarm_grid_mask_80p == 0));

for k = 1:length(domain)
    for i = 1:numFrames
        frame_u = palmer_u_east_surface(:,:,i);
        values_u = frame_u(domain);
        frame_v = palmer_v_north_surface(:,:,i);
        values_v = frame_v(domain);
        mag = sqrt(values_u(k)^2 + values_v(k)^2);
        palmer_mag_surface_selected(i)=mag;
    end
    test=(isnan(palmer_mag_surface_selected));
    if sum(test) > 600 % catches bad gridpoints
        noralizedACF_all_mag(:,k) = NaN(size(normalizedACF));
        lags_all_mag(:,k) = NaN(size(lags));
        unnormalizedACF_all_mag(:,k) = NaN(size(unnormalizedACF));
    else
        [normalizedACF, lags] = autocorr(palmer_mag_surface_selected,'NumLags',1000);
        unnormalizedACF = normalizedACF*var(palmer_mag_surface_selected,1);
        normalizedACF_all_mag(:,k) = normalizedACF;
        lags_all_mag(:,k) = lags;
        unnormalizedACF_all_mag(:,k) = unnormalizedACF;
    end
    clear palmer_mag_surface_selected
end
normalizedACF_mag = nanmean(normalizedACF_all_mag,2);

for i = 1:length(normalizedACF_mag)
    SEM = std(normalizedACF_all_mag(i,:), 'omitnan')/sqrt(length(normalizedACF_all_mag(i,:)));
    ts = tinv([0.025  0.975],length(normalizedACF_all_mag(i,:))-1);
    ci = nanmean(normalizedACF_all_mag(i,:),2) + ts.*SEM;
    ci_normalized_acf_all_mag(i,:) = ci;
    ci_dist = abs(ci-nanmean(normalizedACF_all_mag(i,:)));
    ci_dist_normalized_acf_all_mag(i,:) = ci_dist;
end
