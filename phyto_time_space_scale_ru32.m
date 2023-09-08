% PLDP Phytoplankton residence time from RU32 stationary glider
% Jackie Veatch 01Nov2021

% load in version of RU32 data created by 'glider_SWARM_GDA.m' code
load('/Users/jveatch/Documents/MATLAB/SWARM/glider/data/ru32-20200111T1444-profile-sci-delayed_32fc_30cc_bc41.mat');

addpath('/Users/jveatch/Documents/MATLAB/SWARM/glider/code/spt-master/bin/')
ru32 = ru32_20200111T1444_profile_sci_;
time_ru32 = ru32.time;
ru32_start = datestr(epoch2datenum(time_ru32(1)));
ru32_finish =  datestr(epoch2datenum(time_ru32(end)));

%%  get rid of flagged values
ru32_vars=fields(ru32);
for v=1:length(ru32_vars)
    ind_bad = find(ru32.(ru32_vars{v})==999);
    ru32.(ru32_vars{v})(ind_bad) = NaN;
end

%% when data is downloaded from Erdapp, it got re-ordered somehow. This code
% should put them in back into chronological order
ru32.dnum = epoch2datenum(ru32.time);
[~,ru32_order]=sort(ru32.dnum);
for v=1:length(ru32_vars)
    ru32.(ru32_vars{v})=ru32.(ru32_vars{v})(ru32_order,:);
end

%% 

profiles = unique(ru32.profile_id, 'stable');
for i = 1:length(profiles)
    ind = find(ru32.profile_id == profiles(i));
    if length(ind) <= 10
        ru32.prof_lat(i) = NaN;
        ru32.prof_lon(i) = NaN;
        ru32.prof_dnum(i) = NaN;
    else
    lat = ru32.latitude(ind);
    lon = ru32.longitude(ind);
    dnum = ru32.dnum(ind);
    ru32.prof_lat(i) = mean(lat);
    ru32.prof_lon(i) = mean(lon);
    ru32.prof_dnum(i) = mean(dnum);
    end
end 

test = datestr(ru32.prof_dnum);

%% calculate mld from bfrq
addpath('/Users/jveatch/Documents/MATLAB/seawater_ver3_3/');

% take every 10th datapoint to smooth bfrq a bit
ru32.s_red = ru32.salinity(1:10:length(ru32.dnum));
ru32.s_red = double(ru32.s_red);
ru32.t_red = ru32.temperature(1:10:length(ru32.dnum));
ru32.t_red = double(ru32.t_red);
ru32.p_red = ru32.pressure(1:10:length(ru32.dnum));
ru32.p_red = double(ru32.p_red);
ru32.lat_red = ru32.latitude(1:10:length(ru32.dnum));
ru32.lat_red = double(ru32.lat_red);
ru32.profile_red = ru32.profile_id(1:10:length(ru32.dnum));
ru32.depth_red = ru32.depth(1:10:length(ru32.dnum));

ru32.bfrq = sw_bfrq(double(ru32.s_red), double(ru32.t_red), double(ru32.p_red), ru32.lat_red);
ru32.bfrq(length(ru32.s_red))= NaN; % last datapoint is always missing from bfrq

for i = 1:length(ru32.bfrq)
    if ru32.bfrq(i) == Inf
        ru32.bfrq(i) = NaN;
    else
    end
     if ru32.bfrq(i) == -Inf
        ru32.bfrq(i) = NaN;
    else
    end
end

for i = 1:length(profiles)
    prof_num = profiles(i);
    prof = find(ru32.profile_red == prof_num); % find one MLD for each prof
    if length(prof) < 10 % if prof isn't long enough, dont calc MLD
        ru32.mld(i) = NaN;
    else
        x = max(ru32.bfrq(prof)); % depth of max bfrq = MLD
        if isnan(x)
            ru32.mld(i) = NaN;
        else
        ind_max = find(ru32.bfrq(prof) == x);
        ru32.mld(i) = ru32.depth_red(prof(ind_max(1))); % may create a "wavy" phenomenon, indexing the first max
        end
    end
end

%% calculate mixed layer chlorophyll

ru32.ml_chl_avg = NaN(size(profiles));

for i = 1:length(profiles)
    ind = find(ru32.profile_id == profiles(i));
    chl = ru32.chlorophyll_a(ind);
    mixed_layer_ind = find(ru32.depth(ind) < ru32.mld(i));
    chl_ml = chl(mixed_layer_ind);
    if nansum(chl_ml) == 0
        ru32.ml_chl_avg(i) = NaN;
    else
    ru32.ml_chl_avg(i) = nanmean(chl_ml);
    end
end


%% how close are the glider profiles? --> definitely within CODAR resolution
% addpath('/Users/jveatch/Documents/MATLAB/SWARM/ACROBAT/CODE/'); %location of dist_lat_lon
% lat = ru32.prof_lat(100:200);
% lon = ru32.prof_lon(100:200);
% dist = NaN(size(lon));
% 
% for j = 1:length(lat)-1
%     [dist(j), AF, AR] = dist_lat_lon([lat(j); lat(j+1)], [lon(j); lon(j+1)]);
% end
% 
% mean_dist_glider = nanmean(dist);

%% polygon of points on Head of Canyon
ind_first_deplyment = find(ru32.prof_dnum <= datenum('24-Feb-2020 20:00:00'));

first_deployment_lon = ru32.prof_lon(ind_first_deplyment);
first_deployment_lat = ru32.prof_lat(ind_first_deplyment);
first_deployment_time = ru32.prof_dnum(ind_first_deplyment);

xq = [-64.14, -64.01, -64.1, -64.24, -64.14];
yq = [-64.8, -64.83, -64.88, -64.85, -64.8];

ind_hoc = inpolygon(first_deployment_lon, first_deployment_lat, xq,yq);

stationary_points_lon = first_deployment_lon(ind_hoc);
stationary_points_lat = first_deployment_lat(ind_hoc);
stationary_points_time = first_deployment_time(ind_hoc);


%% define a threshold of the bloom

chl = ru32.ml_chl_avg(ind_hoc);
med = nanmedian(chl);
threshold = med + (0.05*med);
bloom_binary = NaN(size(chl));

for i = 1: length(chl)
    if chl(i) >= threshold
        bloom_binary(i) = 1;
    else
        bloom_binary(i) = 0;
    end
    if isnan(chl(i))
        bloom_binary(i) = NaN;
    else
    end
end

ind_non_nan = ~isnan(bloom_binary);
bloom_binary_non_nan = bloom_binary(ind_non_nan);
time_non_nan = stationary_points_time(ind_non_nan);
%% calculate time of bloom patch

% timestamp of size N --> time_non_nan
% logical array 0- no bloom , and 1- bloom of size N -->
% bloom_binary_non_nan
start_time_array = [];
end_time_array = [];
id = []; % id of the bloom patch
duration = []; % the time duration of bloom patch measures by the glider

ID_BLOOM = 0; % the current ID of the bloom patch
T_BLOOM_START = NaN; % the time reference of the BLOOM START
T_BLOOM_END = NaN; % the time reference of the BLOOM_END
BLOOM_DURATION = NaN;
previous = 1;

for i = 2: length(bloom_binary_non_nan) % START OF THE LOOP
current = i;

if bloom_binary_non_nan(previous) == false && bloom_binary_non_nan(current) == false
% NOBLOOM
T_BLOOM_START = NaN;
T_BLOOM_END = NaN;
end

if bloom_binary_non_nan(previous) == false && bloom_binary_non_nan(current) == true
% BLOOM START
T_BLOOM_START = (time_non_nan(previous)+time_non_nan(current))/2;
end

    if bloom_binary_non_nan(previous) == true && bloom_binary_non_nan(current) == false && ~isnan(T_BLOOM_START)
    % BLOOM END
    T_BLOOM_END = (time_non_nan(previous)+time_non_nan(current))/2;

    ID_BLOOM = ID_BLOOM + 1;
    id = [id ID_BLOOM];
    BLOOM_DURATION = T_BLOOM_END - T_BLOOM_START;
    duration = [duration BLOOM_DURATION];
    end_time_array = [end_time_array T_BLOOM_END];
    start_time_array = [start_time_array T_BLOOM_START];
    % reset
    T_BLOOM_START = NaN;
    T_BLOOM_END = NaN;
    BLOOM_DURATION = NaN;

    end

if bloom_binary_non_nan(previous) == true && bloom_binary_non_nan(current) == true
% BLOOM
end
previous = current;

end % END OF THE LOOP

histogram(duration);
mean_patch_duration = nanmean(duration)*24; % convert to hours

%% test plots, first deployment
figure(1)
ru32.ml_chl_avg(find(ru32.ml_chl_avg==0)) = NaN;
t = datetime(stationary_points_time,'ConvertFrom','datenum');
scatter(t, chl, 45, [0, 0.5, 0],'filled');
yline(threshold);
title('stationary glider mixed layer chlorophyll values');
ylabel('average mixed layer chlorophyll');
xtickformat('dd-MMM-yyyy')

figure(2)
time = datetime(time_non_nan,'ConvertFrom','datenum');
plot(time, bloom_binary_non_nan);
ylim([0 1.2]);
xtickformat('dd-MMM-yyyy')
title('stationary glider mixed layer bloom binary');

figure(3)
t_start = datetime(start_time_array,'ConvertFrom','datenum');
plot(t_start, duration);
ylabel('duration of phytoplankton patch (days)');
xlabel('start time of patch');
title('duration of phytoplankton patch over first deployment');

%% small segment
close
figure(1)
ru32.ml_chl_avg(find(ru32.ml_chl_avg==0)) = NaN;
t = datetime(stationary_points_time,'ConvertFrom','datenum');
scatter(t(400:500), chl(400:500), 45, [0, 0.5, 0],'filled');
yline(threshold);
title('stationary glider mixed layer chlorophyll values');
ylabel('average mixed layer chlorophyll');
xtickformat('dd-MMM-yyyy')

%% what are the mag/direction of surface currents in HOC

load('/Volumes/T7_Shield/jmv208/SWARM_data/SWARM_CODAR.mat');
xq = [-64.14, -64.01, -64.1, -64.24, -64.14];
yq = [-64.8, -64.83, -64.88, -64.85, -64.8];
[x,y] = meshgrid(CODAR.lon, CODAR.lat);
ind_hoc_codar = inpolygon(x, y, xq,yq);
ind_time_first_deplyment_CODAR = find(CODAR.dnum >= stationary_points_time(1) & CODAR.dnum <= stationary_points_time(end));
ind_patch_all = [];

% create index of times while patches are present
for i = 1: length(start_time_array)
    ind_patch = find(CODAR.dnum >= start_time_array(i) & CODAR.dnum <= end_time_array(i));
    ind_patch_all = [ind_patch_all, ind_patch];
end

HOC_u = CODAR.u(:,:,ind_patch_all);
for i = 1:length(ind_patch_all)
    u = HOC_u(:,:,i);
    HOC_u_all(:,:,i) = u(ind_hoc_codar);
end
HOC_v = CODAR.v(:,:,ind_patch_all);
for i = 1:length(ind_patch_all)
    v = HOC_v(:,:,i);
    HOC_v_all(:,:,i) = v(ind_hoc_codar);
end

for i = 1:length(ind_patch_all)
    u = HOC_u_all(:,i);
    v = HOC_v_all(:,i);
    mag = (u.^2 + v.^2).^0.5;
    mag_all(i) = mean(mag);
    for j = 1:length(u)
        dir(j) = atan2d(u(j),v(j));
    end
    dir_all(i) = mean(dir);
end

figure(1)
histogram(mag_all);
title('magnitude of CODAR at stationary glider point during patches');
xlabel('magnitude (cm/s)');

figure(2)
histogram(dir_all);
title('direction of CODAR at stationary glider point duration patches');
xlabel('degrees');

figure(3)
compass(HOC_u_all, HOC_v_all);
title('compass of CODAR at station glider during patches');

%% calculate radius of patch using average CODAR velocity during patch * time of patch

for i = 1:length(start_time_array)
    ind_patch = find(CODAR.dnum >= start_time_array(i) & CODAR.dnum <= end_time_array(i));
    u_patch = CODAR.u(:,:,ind_patch); % index in time
    u_patch_hoc = u_patch(ind_hoc_codar); % index in space
    v_patch = CODAR.v(:,:,ind_patch); % index in time
    v_patch_hoc = v_patch(ind_hoc_codar); % index in space
    mag_patch_hoc = (u_patch_hoc.^2 + v_patch_hoc.^2).^0.5;
    patch_radius(i) = duration(i)*24*60*60 * nanmean(mag_patch_hoc);
end


%% plot path of RU32
addpath(genpath('/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code/antarcticaPlotting'))
figure(1)

title('path of glider');
xlabel('longitude');
ylabel('latitude');

scatter(ru32.longitude,ru32.latitude , 25, 'filled');
hold on

addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code'
    
	bathy=load ('antarctic_bathy_2.mat');
	ind2= bathy.depthi==99999;
	bathy.depthi(ind2)=[];
	bathylines1=0:-10:-100;
	bathylines2=0:-200:-1400;
	bathylines=[bathylines2];
	
	[cs, h1] = contour(bathy.loni,bathy.lati, bathy.depthi,bathylines, 'linewidth', .25);
	clabel(cs,h1,'fontsize',6);
	set(h1,'LineColor','black')

tanLand = [240,230,140]./255;
    S1 = shaperead('cst00_polygon_wgs84.shp');
    S2=S1(1:1174);
    ind=[0,find(isnan(S1(1175).X))];
    for x=1:length(ind)-1
        S2(1174+x)=S1(1175);
        S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
        S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
    end
    mapshow(S2,'facecolor', tanLand)

    hold on
 %marking location of CODAR stations
    s(1) = plot(-64.0554167, -64.7741833, 'g^',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(2) = plot(-64.3604167, -64.7871667, 'gs',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(3) = plot(-64.0446333, -64.9183167, 'gd',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');


adelie_lat=[-64.7977; -64.8302; -64.8412; -64.8080; -64.8196; -64.8536; -64.8652; -64.8307; -64.8403; -64.8735]; 
adelie_lon=[-64.1922; -64.0887; -64.1052; -64.2135; -64.2350; -64.1283; -64.1477; -64.2546; -64.2738; -64.1670];

gentoo_lat=[-64.9027; -64.8818; -64.8567; -64.8801; -64.8799; -64.8392; -64.8251; -64.8649; -64.8522; -64.8110]; 
gentoo_lon=[-64.0290; -64.0647; -63.9853; -63.9461; -63.8538; -63.9279; -63.8836; -63.8059; -63.7690; -63.8438];

plot(adelie_lon, adelie_lat, 'g-');
plot(gentoo_lon, gentoo_lat, 'g-'); 
hold on
    ylim([-65.05 -64.75])
    xlim([-64.45 -63.75])
    
title('path of RU32');

%% index for time at head of canyon

%1/11 - 2/19
%2/20-3/11
%2/22 along canyon
%3/4 deep across canyon line

ind_first_deployment = find(ru32.prof_dnum >= datenum('11-Jan-2020 00:00:00') & ru32.prof_dnum <= datenum('20-Feb-2020 00:00:00'));
% plot(ru32.prof_time(ind_first_deployment), ru32.prof_lat(ind_first_deployment));

xlabel('longitude');
ylabel('latitude');

scatter(ru32.prof_lon(ind_first_deployment),ru32.prof_lat(ind_first_deployment), 25, ru32.prof_time(ind_first_deployment), 'filled');
hold on

colorbar()
% caxis([

addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code'
    
	bathy = load('antarctic_bathy_2.mat');
	ind2= bathy.depthi==99999;
	bathy.depthi(ind2)=[];
	bathylines1=0:-10:-100;
	bathylines2=0:-200:-1400;
	bathylines=[bathylines2];
	
	[cs, h1] = contour(bathy.loni,bathy.lati, bathy.depthi,bathylines, 'linewidth', .25);
	clabel(cs,h1,'fontsize',6);
	set(h1,'LineColor','black')

tanLand = [240,230,140]./255;
    S1 = shaperead('cst00_polygon_wgs84.shp');
    S2=S1(1:1174);
    ind=[0,find(isnan(S1(1175).X))];
    for x=1:length(ind)-1
        S2(1174+x)=S1(1175);
        S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
        S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
    end
    mapshow(S2,'facecolor', tanLand)

    hold on
 %marking location of CODAR stations
    s(1) = plot(-64.0554167, -64.7741833, 'g^',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(2) = plot(-64.3604167, -64.7871667, 'gs',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(3) = plot(-64.0446333, -64.9183167, 'gd',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');


adelie_lat=[-64.7977; -64.8302; -64.8412; -64.8080; -64.8196; -64.8536; -64.8652; -64.8307; -64.8403; -64.8735]; 
adelie_lon=[-64.1922; -64.0887; -64.1052; -64.2135; -64.2350; -64.1283; -64.1477; -64.2546; -64.2738; -64.1670];

gentoo_lat=[-64.9027; -64.8818; -64.8567; -64.8801; -64.8799; -64.8392; -64.8251; -64.8649; -64.8522; -64.8110]; 
gentoo_lon=[-64.0290; -64.0647; -63.9853; -63.9461; -63.8538; -63.9279; -63.8836; -63.8059; -63.7690; -63.8438];

plot(adelie_lon, adelie_lat, 'g-');
plot(gentoo_lon, gentoo_lat, 'g-'); 
hold on
    ylim([-65.05 -64.75])
    xlim([-64.45 -63.75])
    
title('path of RU32');


%% plot to check indexing

figure(3)
xlabel('longitude');
ylabel('latitude');

scatter(ru32.prof_lon(ind_hoc),ru32.prof_lat(ind_hoc), 25, 'filled', 'r');
hold on
scatter(ru32.prof_lon(~ind_hoc),ru32.prof_lat(~ind_hoc), 25, 'filled', 'b');

% scatter(ru32.prof_lon, ru32.prof_lat, 25, 'filled', 'b');
% hold on;

addpath '/Users/jveatch/Documents/MATLAB/Particle_Track_Code/Matlab_Code'
    
	bathy = load('antarctic_bathy_2.mat');
	ind2= bathy.depthi==99999;
	bathy.depthi(ind2)=[];
	bathylines1=0:-10:-100;
	bathylines2=0:-200:-1400;
	bathylines=[bathylines2];
	
	[cs, h1] = contour(bathy.loni,bathy.lati, bathy.depthi,bathylines, 'linewidth', .25);
	clabel(cs,h1,'fontsize',6);
	set(h1,'LineColor','black')

tanLand = [240,230,140]./255;
    S1 = shaperead('cst00_polygon_wgs84.shp');
    S2=S1(1:1174);
    ind=[0,find(isnan(S1(1175).X))];
    for x=1:length(ind)-1
        S2(1174+x)=S1(1175);
        S2(1174+x).X=S2(1174+x).X(ind(x)+1:ind(x+1));
        S2(1174+x).Y=S2(1174+x).Y(ind(x)+1:ind(x+1));
    end
    mapshow(S2,'facecolor', tanLand)

    hold on
 %marking location of CODAR stations
    s(1) = plot(-64.0554167, -64.7741833, 'g^',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(2) = plot(-64.3604167, -64.7871667, 'gs',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');
    s(3) = plot(-64.0446333, -64.9183167, 'gd',...
        'markersize', 12,...
        'markerfacecolor', 'green',...
        'markeredgecolor', 'black');


adelie_lat=[-64.7977; -64.8302; -64.8412; -64.8080; -64.8196; -64.8536; -64.8652; -64.8307; -64.8403; -64.8735]; 
adelie_lon=[-64.1922; -64.0887; -64.1052; -64.2135; -64.2350; -64.1283; -64.1477; -64.2546; -64.2738; -64.1670];

gentoo_lat=[-64.9027; -64.8818; -64.8567; -64.8801; -64.8799; -64.8392; -64.8251; -64.8649; -64.8522; -64.8110]; 
gentoo_lon=[-64.0290; -64.0647; -63.9853; -63.9461; -63.8538; -63.9279; -63.8836; -63.8059; -63.7690; -63.8438];

plot(adelie_lon, adelie_lat, 'g-');
plot(gentoo_lon, gentoo_lat, 'g-'); 
hold on
    ylim([-65.05 -64.75])
    xlim([-64.45 -63.75])
    
title('index stationary points RU32');
