%% NOTES - TODO
results = [];

%% Define inputs - Agarwal verification studies
folder = 'robottip';  % Cylinder, Simple, PlateAnchor or RobotTip
object = 'shaft1';  % Name of stl
triangle_size_calculation = 'Normal';  % 'Fine', 'Normal', 'Rough', 'VeryRough'
triangle_size_visualization = 'Rough';  % 'Fine', 'Normal', 'Rough', 'VeryRough'
rotation = true;  % true or false
linear_velocity = 0.000;  % linear velocity in m/s
direction_angle_xz = -90 * pi / 180;  % angle between direction and x-z-axis
direction_angle_y = -90 * pi / 180;  % angle between direction and y-axis
angular_velocity = [0, 0, -4/3*pi];  % angular velocity in rad/s
rho_c = 1450;  % bulk density of the sand in kg/m³   
mu_int = 1.07;  % internal friction coefficient of the sand
mu_surf = 0.4;  % intruder-surface interaction coefficient
gravity = 9.81;  % gravity in m/s²
depth = 0.12;  % in m

%direction_vector = [1 1 0];
direction_vector = [round(cos(direction_angle_xz), 15) round(cos(direction_angle_y), 15) round(sin(direction_angle_xz), 15)];

%% Plot options
show_geometry = false;
show_direction = false;

show_f_quiver = false;
show_alpha = false;

show_f_scatter = true;
show_f_scatterxyz = false;

show_linear_f = false;

saveFigures = false;

unit_test = false;

%% Read .stl file
TRG = stlread(strcat('./', folder, '/Models/', object, triangle_size_calculation, '.stl'));  % Mesh size for calculation
TRG_visual = stlread(strcat('./', folder, '/Models/', object, triangle_size_visualization, '.stl'));  % Mesh size for force plots

TRG = rotate_triangulation_x(TRG, -0);  % Rotate TRG object
TRG_visual = rotate_triangulation_x(TRG_visual, -0);

TRG = move_triangulation_z(TRG, depth);  % Align bottom of object with depth input
TRG_visual = move_triangulation_z(TRG_visual, depth);

points = (incenter(TRG)').';
normals = (faceNormal(TRG)').';
area = (generateArea(TRG.Points', TRG.ConnectivityList')).';

%% Compute forces using 3D-RFT function
[c_inc, v_norm_vec, F, f, forces_x, forces_y, forces_z, forces, T, torque_x, torque_y, torque_z, alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = RFT3Dfunc(points, normals, area, rotation, angular_velocity, linear_velocity, direction_vector, rho_c, mu_int, mu_surf, gravity, unit_test);

%% Plots
RFTPlots

col_results = [depth; forces_x; forces_y; forces_z; torque_x; torque_y; torque_z];
results = horzcat(results, col_results);

function areaarray = generateArea(Points,List)
    % Compute the side lengths of the triangles once
    a = vecnorm(Points(:, List(1,:)) - Points(:, List(2,:)));
    b = vecnorm(Points(:, List(2,:)) - Points(:, List(3,:)));
    c = vecnorm(Points(:, List(1,:)) - Points(:, List(3,:)));
    s = (a + b + c) / 2;
    areaarray = sqrt(s .* (s - a) .* (s - b) .* (s - c));
end

function TRG = rotate_triangulation_x(TRG, theta)
    % Create the rotation matrix
    R = [1 0 0; 0 round(cosd(theta), 15) -round(sind(theta), 15); 0 round(sind(theta), 15) round(cosd(theta), 15)];
    % Rotate the points in the triangulation object
    Points = TRG.Points * R;
    
    % Create a new triangulation object with the rotated points and the same connectivity list
    TRG = triangulation(TRG.ConnectivityList, Points);
end

function TRG = move_triangulation_z(TRG, depth)
    minZ = min(TRG.Points(:, 3));
    Points = TRG.Points;
    Points(:, 3) = Points(:, 3) - minZ - depth*1000;
    
    % Create a new triangulation object with the shifted points and the same connectivity list
    TRG = triangulation(TRG.ConnectivityList, Points);
end

