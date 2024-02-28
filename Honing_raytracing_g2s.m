clear
month=4;
date=1;
R12=1; 
hour=3;
minute=40;

UT = [2015,month,date,hour,minute];
speed_Nf_light = 2.99792458e8;
R12 = R12;
doppler_flag = 0;

raytrace_params.origin_lat = 52.16;             % latitude of the start point of rays
raytrace_params.origin_lon = -106.53;      % longitude of the start point of rays
raytrace_params.origin_alt =0.494; % altitude of the start point of rays
raytrace_params.tol=[1e-8 0.005 1];
raytrace_params.nhops=2;% number of hops attempted
raytrace_params.init_delta_elev=1;% elevation angle perturb
raytrace_params.init_delta_bear=10;% bearing angle perturb
raytrace_params.fr=17.5; % frequency in MHz

Newton_params.elev_tol=0.13; % lowest elevation angle attempted
Newton_params.iteration_num=10; % number of Newton's iteration
Newton_params.rece_dis=100; % size of receiving zone (m)

ht_start = 90;          % start height for ionospheric grid (km)
ht_inc = 5;             % height increment step length (km)
num_ht = 100;
lat_start = 0;
lat_inc = 1;
num_lat = 90;
lon_start= -140;
lon_inc = 1;
num_lon = 90;

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc =5;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = ...
            gen_iono_grid_3d(UT, R12, iono_grid_parms, ...
            geomag_grid_parms, doppler_flag);

save(['ionosphere','_',num2str(month),'_',num2str(date),'_',num2str(hour),'.mat'],'iono_pf_grid', 'iono_pf_grid_5', 'collision_freq', 'Bx', 'By', 'Bz');

%%
clear
raytrace_params.origin_lat = 52.16;             % latitude of the start point of rays
raytrace_params.origin_lon = -106.53;            % longitude of the start point of rays
raytrace_params.origin_alt =0.494; % altitude of the start point of rays
raytrace_params.tol=[1e-8 0.005 1];
raytrace_params.nhops=2;% number of hops attempted
raytrace_params.init_delta_elev=1;% elevation angle perturb
raytrace_params.init_delta_bear=10;% bearing angle perturb
raytrace_params.fr=17.5; % frequency in MHz

Newton_params.elev_tol=0.13; % lowest elevation angle attempted
Newton_params.iteration_num=10; % number of Newton's iteration
Newton_params.rece_dis=100; % size of receiving zone (m)

ht_start = 90;          % start height for ionospheric grid (km)
ht_inc = 5;             % height increment step length (km)
num_ht = 100;
lat_start = 0;
lat_inc = 1;
num_lat = 90;
lon_start= -140;
lon_inc = 1;
num_lon = 90;

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ...
    ht_start, ht_inc, num_ht ];

B_ht_start = ht_start;          % start height for geomagnetic grid (km)
B_ht_inc =5;                  % height increment (km)
B_num_ht = ceil(num_ht .* ht_inc ./ B_ht_inc);
B_lat_start = lat_start;
B_lat_inc = 1.0;
B_num_lat = ceil(num_lat .* lat_inc ./ B_lat_inc);
B_lon_start = lon_start;
B_lon_inc = 1;
B_num_lon = ceil(num_lon .* lon_inc ./ B_lon_inc);

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, ...
    B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht];

load('E:\projectnew\raytracing_subodh\2015-04-01\G2S\ionosphere_4_1_3.mat')
lat=h5read('E:\projectnew\raytracing_subodh\2015-04-01\G2S\RRI_20150401_033844_034241_lv1_13.1.0.h5','/CASSIOPE Ephemeris/Geographic Latitude (deg)');
lon=h5read('E:\projectnew\raytracing_subodh\2015-04-01\G2S\RRI_20150401_033844_034241_lv1_13.1.0.h5','/CASSIOPE Ephemeris/Geographic Longitude (deg)');
alt=h5read('E:\projectnew\raytracing_subodh\2015-04-01\G2S\RRI_20150401_033844_034241_lv1_13.1.0.h5','/CASSIOPE Ephemeris/Altitude (km)')*1000;% spacecraft height (m)

satellite_location.lon=lon;
satellite_location.lat=lat;
satellite_location.alt=alt;
% load('D:\Doppler_newton_method\Dec_8\Eclips_inosphere_final_with_B_Dec_8_6_30.mat')
% Bx=B_x;
% By=B_y;
% Bz=B_z;
%Ne_5=zeros(size(Ne));
%  Ne(90-76:90-70,35:38,:)=Ne(90-76:90-70,35:38,:)*2;
%      Ne=Ne*factor;
 Ne= iono_pf_grid.^2 / 80.6164e-6;
 Ne_5 = iono_pf_grid_5.^2 / 80.6164e-6;

% Set OX_mode to -1
raytrace_params.OX_mode = 1;
[~, ray_O_receive] = Honing_raytracing(satellite_location, ...
    Newton_params, raytrace_params, iono_grid_parms, geomag_grid_parms, Ne, Ne_5, collision_freq, ...
    Bx, By, Bz);

% Set OX_mode to 1
raytrace_params.OX_mode = -1;
[~, ray_X_receive] = Honing_raytracing(satellite_location, ...
    Newton_params, raytrace_params, iono_grid_parms, geomag_grid_parms, Ne, Ne_5, collision_freq, ...
    Bx, By, Bz);

% Save both results in a single call
save('E:\projectnew\raytracing_subodh\2015-04-01\G2S\ground2satellite1', 'ray_O_receive', 'ray_X_receive');
%%
 load('E:\projectnew\raytracing_subodh\2015-04-01\G2S\ray_received_2015_4_1_40.mat')
 for k=1:111
    for l=1:98
        x1=ray_O_receive(k).lat;
        x2=ray_X_receive(l).lat;
        y1=ray_O_receive(k).lon;
        y2=ray_X_receive(l).lon;

        %[xx,yy]=meshgrid(x,y);
        O=ray_O_receive(k).phase_path();
        X=ray_X_receive(l).phase_path();
    
        %surf(xx,yy,z);
        plot3(x1,y1,O,"Color",'r');
        hold on
        plot3(x2,y2,X,"color", 'b');
        xlabel('Latitude')
        ylabel('Longitude')
        zlabel('phase path')
        title('X and O Rays received at NJIT(2018-03-20)')
        legend('O-ray','X-ray','Location','bestoutside')
    end
end

