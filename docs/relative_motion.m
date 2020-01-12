%% Import data and create variables
clear all
close all
clc 

location = uigetdir;
Runs = dir([location,'\*.csv']); 
numruns = length(Runs);

filenames = {Runs(:).name}.';
accel_data = cell(1,numruns);

for i = numruns:-1:1
    fullname = Runs(i).name; 
    accel_data{i} = xlsread(fullname); 
    accel_data{i}(:,3:2:end)=[];  % get rid of extra time columns
    accel_data{i}(1:3,:) = []; % get rid of headers
end

%% Filter, get velocities and displacements, plot.
clc

accel_num = 4; % Change this depending on # of accels.
accel_dirs = 3*accel_num; % Change this if unidirectional.
vels = cell(1,numruns+1);
disps = cell(1,numruns+1);
accel_filt = cell(1,numruns);
vel_filt = cell(1,numruns);
t = cell(1,numruns);
g = 9.80665;

fpass = 20; % Change to adjust filtered frequences in Hz.
inc = 0.00015625; % Change based on sampling increment in seconds.
fs = 1/inc; % Sampling rate in Hz.

% Filter and get velocity
for i = 1:numruns
    accel_filt{i} = highpass(accel_data{i}(:,2:end),fpass,fs);
    t{i} = accel_data{i}(:,1);
    disps{i}(:,1) = t{i};
    vels{i}(:,1) = t{i};
    for j = 1:accel_dirs
        accel_filt{i}(:,j) = accel_filt{i}(:,j) - mean(accel_filt{i}(:,j));
        dir_accel = accel_filt{i}(:,j) * g; % accel. in each direction in m/s2
        vels{i}(:,j+1) = cumtrapz(t{i},dir_accel);
    end
end

% Filter and get displacement
for i = 1:numruns
    vel_filt{i} = highpass(vels{i}(:,2:end),fpass,fs);
    t{i} = accel_data{i}(:,1);
    for j = 1:accel_dirs
        vel_filt{i}(:,j) = vel_filt{i}(:,j) - mean(vel_filt{i}(:,j));
        dir_vel = vel_filt{i}(:,j); % vel. in each direction in m/s
        disps{i}(:,j+1) = cumtrapz(t{i}, dir_vel);
    end
end

% Plot accel, vel, displacement
    % Plot accel.
    figure(1)
    hold on
    plot(t{2}, accel_filt{2}(:,2:7)*g)
    title('Acceleration')
    xlabel('time (s)')
    ylabel('m/s^2')

    % Plot velocity.
    figure(2)
    hold on
    plot(t{2}, vel_filt{2}(:,2:7))
    title('Velocity')
    xlabel('time (s)')
    ylabel('m/s')

    % Plot displacement.
    figure(3)
    hold on
    plot(t{2}, disps{2}(:,2:7))
    title('Displacement')
    xlabel('time (s)')
    ylabel('m')



%% Get relative displacements of interest

% Rear camera mount vs. central lidar relative displacements in car ref. frame (x, y, z).
% driver side (x+) , up (y+), front of car (z+)
% all units converted to mm

[lidar_dx, lidar_dy, lidar_dz, cam_dx, cam_dy, cam_dz, dx, dy, dz] = deal(cell(1,numruns));
[max_t_dx, max_t_dy, max_t_dz] = deal(cell(1,numruns));

% Get individual directional displacements and compare. 
for i = 1:numruns
    % Get individual directional displacements. Change from accel frame to car ref. frame.
    lidar_dx{i} = disps{i}(:,2) * 1000;
    lidar_dy{i} = -disps{i}(:,3) * 1000; 
    lidar_dz{i} = -disps{i}(:,4) * 1000;
    cam_dx{i} = cosd(36) * disps{i}(:,5) * 1000;
    cam_dy{i} = disps{i}(:,6) * 1000;
    cam_dz{i} = cosd(36) * disps{i}(:,7) * 1000;
    
    % Relative displacements over time.
    dx{i} = abs(lidar_dx{i} - cam_dx{i});
    dy{i} = abs(lidar_dy{i} - cam_dy{i});
    dz{i} = abs(lidar_dz{i} - cam_dz{i});
    
    % Maximum relative displacements correlated over time for each run.
    max_t_dx{i} = max(dx{i});
    max_t_dy{i} = max(dy{i});
    max_t_dz{i} = max(dz{i});
    
end

% Maximum relative displacement correlated over time for all runs.
max_dx = max(cell2mat(max_t_dx));
max_dy = max(cell2mat(max_t_dy));
max_dz = max(cell2mat(max_t_dz));

max_relative_displacements_mm = [max_dx,max_dy,max_dz]

% Test plot relative displacements.
    figure
    hold on
    plot(t{2}, dx{2})
    plot(t{2}, dy{2})
    plot(t{2}, dz{2})
    title('Relative Displacements for One Run')
    legend('rel_x', 'rel_y', 'rel_z')
    xlabel('time')
    ylabel('mm')


%% Get rotations from relative displacement

% Measured distances from LiDAR to camera in mm.
x = .52 * 1000;
y = .18 * 1000;
z = 1.84 * 1000;

% Angle of rotation about x, y, z in degrees.
[theta_x, theta_y, theta_z] = deal(cell(1,numruns));
[max_t_theta_x, max_t_theta_y, max_t_theta_z] = deal(cell(1,numruns));

% Get rotation angle over time.
for i = 1:numruns
    % Relative rotations correlated over time.
    theta_x{i} = atand((y + dy{i}) ./ (z + dz{i})) - atand(y ./ (z + dz{i}));
    theta_y{i} = atand((x + dx{i}) ./ (z + dz{i})) - atand(x ./ (z + dz{i}));
    theta_z{i} = atand((y + dy{i}) ./ (x + dx{i})) - atand(y ./ (x + dx{i}));
    
    % Maximum relative rotations correlated over time for each run.
    max_t_theta_x{i} = max(theta_x{i});
    max_t_theta_y{i} = max(theta_y{i});
    max_t_theta_z{i} = max(theta_z{i});
    
end

% Maximum relative rotation correlated over time for all runs.
max_theta_x = max(cell2mat(max_t_theta_x));
max_theta_y = max(cell2mat(max_t_theta_y));
max_theta_z = max(cell2mat(max_t_theta_z));

max_relative_rotation_degrees = [max_theta_x,max_theta_y,max_theta_z]

% Test plot relative rotation.
    figure
    hold on
    plot(t{2}, theta_x{2})
    plot(t{2}, theta_y{2})
    plot(t{2}, theta_z{2})
    title('Relative Rotations for One Run')
    legend('rel_x', 'rel_y', 'rel_z')
    xlabel('time')
    ylabel('degrees')


