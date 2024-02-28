%
%Name:
% Honing_raytracing
%
%Purpose:
%This is a horning algorithm to do ray tracing between ground and satellite.
%This function is based on PHaLARP function raytrace_3d.m.
%This function requires preloaded ionosphere and geomagnetic information.
%The grid size information can be found in "Input".
%In order to use this function, PHaLARP_4.5.1 is needed and could be downloaded
%from:  https://www.dst.defence.gov.au/our-technologies/pharlap-provision-high-frequency-raytracing-laboratory-propagation-studies
%
% 
% Calling sequence:
%  ray_receive=Horning_raytracing(satellite_lat,satellite_lon,satellite_alt,...
%                                 rece_dis,elev_tol,iteration_num,raytrace_params,...
%                                 iono_grid_parms,Ne,Ne_5,collision_freq,Bx,By,Bz)
%Input:
%   satellite_location  - 3*M matrix with latitude, longitude, and altitude.
%      .lat            - 1*M (where M is the number of locations) array of
%                        geodetic (WGS84) latitude (-90 to 90 degrees) of
%                        the satellite
%      .lon            - 1*M (where M is the number of locations) array of
%                        geodetic (WGS84) longitude (-180 to 180 degrees) of
%                        the satellite
%      .alt             - 1*M (where M is the number of locations) array of
%                        geodetic (WGS84) altitude (m) of the satellite
%    Newton_params      - 1*3 structure, with fields: 
%      .rece_dis        - receiving zone distance (m) from the satellite
%      .elev_tol        - The lowest initial elevation angle.
%      .iteration_num   - Newton's method iteration number normally the value
%                        is 10-20
%      raytrace_params - 1*7 structure containing the parameters which
%                       define the raytracing initial conditions
%
%     .origin_lat      - geodetic (WGS84) latitude (-90 to 90 degrees) of
%                       start point of rays (All rays have the same origin)
%     .origin_long     - geodetic (WGS84) longitude (-180 to 180 degrees) of
%                       start point of rays
%     .origin_height   - WGS84 height above sea-level of a start point of rays which
%                       must be below the start of the ionosphere (km). (If the start of
%                       the ray tracing is inside the ionosphere, then the initial
%                       state vector of the ray must be specified.)
%     .init_delta_elev - initial elevation angle perturbation (degree)
%                        normally set as 0-1
%     .init_delta_bear - initial bearing angle perturbation (degree)
%                        normally set as 1-10
%     .fr              - wave frequency of the rays (MHz)
%     .OX_mode         - polarization mode of rays: 1 = O, -1 = X, 0 = no field
%                       (All rays have the same polarization mode)
%     .nhops           - number of hops to attempt (same for all rays)
%     .tol             - a 1 or 3 element vector controlling ODE solver precision.
%        tol(1)  =  ODE solver tolerance, valid values 1e-12 to 1e-2
%        tol(2)  =  ODE solver minimum step size to consider (0.001 to 1 km)
%        tol(3)  =  ODE solver maximum step size to consider (1 to 100 km)
%     Ne               - 3d grid of ionospheric electron density (electrons/cm^3).
%                      Maximum grid size is 701 X 701 X 401 (WGS84 lat, lon,
%                      height). The ionospheric data and their gradients must
%                      be smooth and continuous. See also note 2 below on memory.
%     Ne_5             - 3d grid of ionospheric electron density (electrons/cm^3)
%                      5 minutes later. This array is used to calculate Doppler
%                      shift. If Doppler shift is not required, then this array
%                      can be set to iono_en_grid. The size of grid must be
%                      identical to iono_en_grid.
%     collision_freq - 3d grid of effective electron collision frequency (Hz).
%                      This is used to calculate the total and deviation
%                      absorption (see notes). If absorption is not required
%, then the array elements can be set to zero. Size of
%                      grid must be identical to iono_en_grid.
%    iono_grid_parms - 9x1 vector containing the parameters which define the
%                       ionospheric grid :
%             (1) geodetic latitude (degrees) of the start of grid - must be in the
%                 range -90 to 90 degrees
%             (2) latitude step (degrees)
%             (3) number of latitudes
%             (4) geodetic longitude (degrees) of the start of grid - must be in
%                 the range -180 to 180 degrees
%             (5) longitude step (degrees)
%             (6) number of longitudes
%             (7) geodetic height (km) start
%             (8) height step (km)
%             (9) number of heights
%
%     Bx, By, Bz        - 3D grids of x, y and z components of the geomagnetic
%                         field (Tesla) as a function of WGS84 latitude,
%                         longitude and height. The maximum grid size is
%                         101 X 101 X 201 (lat, lon, height). The geomagnetic field
%                         data and their gradients must be smooth and continuous.
%
%     geomag_grid_parms - 9x1 vector containing the parameters which define the
%                         ionospheric grid.
%                         array format:
%                       iono_grid_parms = [lat_start, lat_inc, num_lat, ...
%                                          lon_start, lon_inc, num_lon, ...
%                                          ht_start, ht_inc, num_ht];
%             (1) geodetic latitude (degrees) of start of grid - must be in the
%                 range -90 to 90 degrees
%             (2) latitude step (degrees)
%             (3) number of latitudes
%             (4) geodetic longitude (degrees)) of start of grid - must be in
%                 the range -180 to 180 degrees
%             (5) longitude step (degrees)
%             (6) number of longitudes
%             (7) geodetic height (km) start
%             (8) height step (km)
%             (9) number of heights
%    Note: Use exactly the same variable name and order for the iono_grid_parms input
%          as shown in the array format.
%
%Output:
% ray_receive(:) - 1 * M structure contained received ray information (total of
%                  M rays). For each ray, the data ends at the receiving point
%                  M is the number of locations of the satellite. Each
%                  field is a  1 * N array containing information of each
%                  ray. The fields of this structure are:
%       .initial_elev            - initial elevation of the ray (degrees)
%       .initial_bearing         - initial bearing of the ray (degrees)
%       .frequency               - carrier frequency of the ray (MHz)
%       .lat                     - geodetic (WGS84) latitude (degrees)
%       .lon                     - geodetic (WGS84) longitude (degrees)
%       .height                  - the height of ray above WGS84 ellipsoid (km)
%       .phase_path              - phase path (km)
%       .group_range             - group range (km)
%       .loca_num                - index number of satellite location along 
%                                  sataltite_location.lat array.
%       .frequency               - carrier frequency of the ray (MHz)
%       .refractive_index        - refractive index
%       .electron_density        - electron density (e- / cm^3)
%       .absorption              - absorption in dB
%
% Reference: James, H. G. (2006). Effects on transionospheric HF propagation
%            observed by ISIS at middle and auroral latitudes. Advances in
%            Space Research, 38(11), 2303-2312.

function [ray_receive_V1,ray_receive]=Honing_raytracing(satellite_location,...
    Newton_params,raytrace_params,iono_grid_parms,geomag_grid_parms,Ne,Ne_5,...
    collision_freq,Bx,By,Bz)

% speed_of_light = 2.99792458e8;
origin_lat=raytrace_params.origin_lat;
origin_lon=raytrace_params.origin_lon;
origin_ht=raytrace_params.origin_alt;% ground transmitter location
satellite_lat=satellite_location.lat;
satellite_lon=satellite_location.lon;
satellite_alt=satellite_location.alt;
rece_dis=Newton_params.rece_dis;
iteration_num=Newton_params.iteration_num;
elev_tol=Newton_params.elev_tol;% set up raytracing parameters.
% find initial elevationa and bearing angle
wgs84 = wgs84Ellipsoid;
[az,elev,~] = geodetic2aer(satellite_lat,satellite_lon,satellite_alt,origin_lat,origin_lon,origin_ht*1000,wgs84);

L=1;
roun_number=nan(1,length(satellite_lat)); % evaluate the iteration efficiency
for loca_num=1:length(satellite_lat)
    initial_elevs=elev(loca_num);
    initial_bear=az(loca_num);
    delta_elevs=raytrace_params.init_delta_elev;
    delta_bearing=raytrace_params.init_delta_bear;% initial elevs and bear variation need to adjust with different conditions.
    if elev(loca_num)<elev_tol %start tratracing procese with initial elevation threshold
        continue
    end
    iono_en_grid = Ne;
    iono_en_grid_5 = Ne_5;
    tic % record time
    for Newton_count=1:20  % maximum 20 rounds of Newton's iteration
        if Newton_count>1 % Use the smallest distance ray from the first-round Newton's iteration as the new initial
            delta_elevs=rand(1);
            delta_bearing=0.5*rand(1);
            [~,resolve_ray_index]=min(sR(:));
            initial_elevs=elevs_array(resolve_ray_index);
            initial_bear=bearing_array(resolve_ray_index);
        end
        elevs_matrix=nan(iteration_num,4);
        bearing_matrix=nan(iteration_num,4);
        sR=nan(iteration_num,4);
        for m=1:iteration_num % Number of Newton's iteration. The distance between the spacecraft and the ray will converge after several iterations.
            elevs=[initial_elevs,initial_elevs+delta_elevs];
            bearing_angle=[initial_bear,initial_bear+delta_bearing];% set up elevation and bearing angle
            elevs_matrix(m,:)=[initial_elevs,initial_elevs+delta_elevs,initial_elevs,initial_elevs+delta_elevs];
            bearing_matrix(m,:)=[initial_bear,initial_bear,initial_bear+delta_bearing,initial_bear+delta_bearing];
            elevs_array=elevs_matrix(:);
            bearing_array=bearing_matrix(:);% record elevation and bearing for updating "inital_elevs" and "initial_bear"
            fr=raytrace_params.fr;
            freqs = ones(size(elevs))*fr;
            i=1;
            ray_inf=struct('lat',[],'lon',[],'elev',[],'bearing',[],'absorption',[],'loca_num'...
                ,[],'electron_density',[],'refractive_index',[],'phase_path',[],'range',[],'group_range',[],'geometric_path_length',[],'Bx',[],'By',[],'Bz',[],'wavenorm_B_angle',[],'polariz_mag',[],'wave_Efield_tilt',[],'volume_polariz_tilt',[]);
            repmat(ray_inf,[1,4]);
            lat_receive=nan([1,4]);
            lon_receive=nan([1,4]);
            for bearnumber=1:length(bearing_angle)% call raytracing_3d at different bearing angle
                ray_bears =ones(size(elevs))*bearing_angle(bearnumber);
                nhops = raytrace_params.nhops;                  % number of hops
                tol =  raytrace_params.tol;      % ODE solver tolerance and min max stepsizes
                
                OX_mode = raytrace_params.OX_mode;
                [ray_data, ray, ~] = ...
                    raytrace_3d(origin_lat, origin_lon, origin_ht, elevs, ray_bears, freqs, ...
                    OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
                    collision_freq, iono_grid_parms, Bx, By, Bz, ...
                    geomag_grid_parms);
                for elevsnumber=1:length(elevs)
                    ray_O_height=ray(elevsnumber).height;
                    height_instability=satellite_alt(loca_num)-ray_O_height*1000;% distance between ray point and satalite altitude plane
                    [~,endpoint]=min(abs(height_instability));% find the minimum distance point(the potential receiving point)
                    %record selected rays information with cutoff at the potential receiving point. 
                    ray_inf(i).lat=ray(elevsnumber).lat(1:endpoint);
                    ray_inf(i).lon=ray(elevsnumber).lon(1:endpoint);
                    ray_inf(i).height=ray(elevsnumber).height(1:endpoint);
                    ray_inf(i).elev=ray(elevsnumber).initial_elev;
                    ray_inf(i).bearing=ray(elevsnumber).initial_bearing;
                    ray_inf(i).absorption=ray(elevsnumber).absorption(1:endpoint);
                    ray_inf(i).group_range=ray(elevsnumber).group_range(1:endpoint);
                    ray_inf(i).loca_num=loca_num;
                    ray_inf(i).electron_density=ray(elevsnumber).electron_density(1:endpoint);
                    ray_inf(i).refractive_index=ray(elevsnumber).refractive_index(1:endpoint);
                    ray_inf(i).phase_path=ray(elevsnumber).phase_path(1:endpoint);
                    ray_inf(i).geometric_path_length=ray_data(elevsnumber).geometric_path_length;
                    ray_inf(i).Bx=ray(elevsnumber).geomag_x(1:endpoint);
                    ray_inf(i).By=ray(elevsnumber).geomag_y(1:endpoint);
                    ray_inf(i).Bz=ray(elevsnumber).geomag_z(1:endpoint);
                    ray_inf(i).wavenorm_B_angle=ray(elevsnumber).wavenorm_B_angle(1:endpoint);
                    ray_inf(i).polariz_mag=ray(elevsnumber).polariz_mag(1:endpoint);
                    ray_inf(i).wave_Efield_tilt=ray(elevsnumber).wave_Efield_tilt(1:endpoint);
                    ray_inf(i).volume_polariz_tilt=ray(elevsnumber).volume_polariz_tilt(1:endpoint);
                    lon_receive(i)=ray_inf(i).lon(end);
                    lat_receive(i)=ray_inf(i).lat(end);% receiving point latitude and longitude
                    i=i+1;
                end
            end
            [~,~,newslantRange] = geodetic2aer(satellite_lat(loca_num),satellite_lon(loca_num),satellite_alt(loca_num),lat_receive,lon_receive,satellite_alt(loca_num),wgs84);
            sR(m,:)=newslantRange;
            [M,receive_ray_index]=min(newslantRange); % find the minimum distance with assumption that ray receiving point at the same altitude as satellite.
            if M>rece_dis %slantRange larger than threshold continue the iteration
                % The iteration equations from James, H. G. (2006) paper
                slope1=(lon_receive(1)-lon_receive(3))/(lat_receive(1)-lat_receive(3));
                slope2=(lon_receive(1)-lon_receive(2))/(lat_receive(1)-lat_receive(2));
                delta_elevs=delta_elevs*(satellite_lon(loca_num)-lon_receive(1)-slope1*(satellite_lat(loca_num)-lat_receive(1)))/(lon_receive(2)-lon_receive(1)-slope1*(lat_receive(2)-lat_receive(1)));
                delta_bearing=delta_bearing*(satellite_lon(loca_num)-lon_receive(1)-slope2*(satellite_lat(loca_num)-lat_receive(1)))/(lon_receive(3)-lon_receive(1)-slope2*(lat_receive(3)-lat_receive(1)));
            else
                %slantRange smaller than threshold the received ray is
                %found
                ray_receive_V1(L).lat=ray_inf(receive_ray_index).lat;% save the received ray information
                ray_receive_V1(L).lon=ray_inf(receive_ray_index).lon;
                ray_receive_V1(L).height=ray_inf(receive_ray_index).height;
                ray_receive_V1(L).elevs=ray_inf(receive_ray_index).elev;
                ray_receive_V1(L).bearing=ray_inf(receive_ray_index).bearing;
                ray_receive_V1(L).loca_num=ray_inf(receive_ray_index).loca_num;
                ray_receive_V1(L).absorption=ray_inf(receive_ray_index).absorption;
                ray_receive_V1(L).group_range=ray_inf(receive_ray_index).group_range;
                ray_receive_V1(L).electron_density=ray_inf(receive_ray_index).electron_density;
                ray_receive_V1(L).phase_path=ray_inf(receive_ray_index).phase_path;
                ray_receive_V1(L).refractive_index=ray_inf(receive_ray_index).refractive_index;
                ray_receive_V1(L).geometric_path_length=ray_inf(receive_ray_index).geometric_path_length;
                ray_receive_V1(L).Bx=ray_inf(receive_ray_index).Bx;
                ray_receive_V1(L).By=ray_inf(receive_ray_index).By;
                ray_receive_V1(L).Bz=ray_inf(receive_ray_index).Bz;
                ray_receive_V1(L).wavenorm_B_angle=ray_inf(receive_ray_index).wavenorm_B_angle;
                ray_receive_V1(L).polariz_mag=ray_inf(receive_ray_index).polariz_mag;
                ray_receive_V1(L).wave_Efield_tilt=ray_inf(receive_ray_index).wave_Efield_tilt;
                ray_receive_V1(L).volume_polariz_tilt=ray_inf(receive_ray_index).volume_polariz_tilt;
                ray_receive_V1(L).range=M;
                L=L+1;
                break % slantRange smaller than the threshold break from iteration loop
            end
        end
        if M<rece_dis % slantRange smaller than the threshold break from 
            break
        end
    end
    roun_number(loca_num)=Newton_count;
    if Newton_count==20
        fprintf('\n initial condition need change \n');
    end
    toc
    loca_num % moniter the process 
end

% select the receiving 
N=1;
for i=1:length(ray_receive_V1) % restrict the receiving point height. 
    loca_num=ray_receive_V1(i).loca_num;
    alt_diff=satellite_alt(loca_num)-ray_receive_V1(i).height(end)*1000;
    if abs(alt_diff)<=rece_dis
        ray_receive(N)=ray_receive_V1(i);
        N=N+1;
    end
end
if isempty(ray_receive(1).lat)
fprintf('\n No ray is received.\n');
end
end