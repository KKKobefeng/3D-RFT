% Iterate through depth 
clear
close
tips = {'tipnr4', 'tipnr5', 'tipnr2', 'tipnr7', 'tipnr8'};

figure
hold on;

start_depth = 0;
end_depth = 0.13;
step_size = 0.01;
num_steps = (end_depth - start_depth)./step_size;

for index = 1:1:5
%% Define inputs - Agarwal verification studies
folder = 'robottip';        % Cylinder, Simple, PlateAnchor or RobotTip
object = tips{index};          % Name of stl
triangle_size_calculation = 'normal';   % 'Fine', 'Normal', 'Rough', 'VeryRough'
triangle_size_visualization = 'veryrough';  % 'Fine', 'Normal', 'Rough', 'VeryRough'
rotation = true;                        % true or false
linear_velocity = 0.01;                % linear velocity in m/s
direction_angle_xz = -90 * pi / 180;    % angle between direction and x-z-axis
direction_angle_y = -90 * pi / 180;     % angle between direction and y-axis
angular_velocity = [0, 0, -8.59106529209622*pi];     % angular velocity in rad/s
rho_c = 1607;               % bulk density of the sand in kg/m³   
mu_int = 1.07;              % internal friction coefficient of the sand
mu_surf = 0.1787;             % intruder-surface interaction coefficient
gravity = 9.81;             % gravity in m/s²

%direction_vector = [1 1 0];
direction_vector = [round(cos(direction_angle_xz), 15) round(cos(direction_angle_y), 15) round(sin(direction_angle_xz), 15)];

unit_test = false;

%% Read .stl file
TRG = stlread(strcat('./', folder, '/Models/', object, triangle_size_calculation, '.stl'));  % Mesh size for calculation
TRG_visual = stlread(strcat('./', folder, '/Models/', object, triangle_size_visualization, '.stl'));  % Mesh size for force plots

TRG = rotate_triangulation_x(TRG, -90);  % Rotate TRG object
TRG_visual = rotate_triangulation_x(TRG_visual, -90);

%% Plot options
show_geometry = false;
show_direction = false;

show_f_quiver = false;
show_alpha = false;

show_f_scatter = 0;
show_f_scatterxyz = false;

show_linear_f = false;

saveFigures = false;


sim_data_dtos = SimDataTransferObject.empty(0,num_steps);

step = 1;
depths = start_depth:step_size:end_depth;
z = zeros(1, num_steps);

for depth = start_depth:step_size:end_depth
    TRG = move_triangulation_z(TRG, depth);  % Align bottom of object with depth input
    TRG_visual = move_triangulation_z(TRG_visual, depth);

    points = (incenter(TRG)').';
    normals = (faceNormal(TRG)').';
    area = (generate_area(TRG.Points', TRG.ConnectivityList')).';

    [c_inc, v_norm_vec, F, f, forces_x, forces_y, forces_z, forces, T, torque_x, torque_y, torque_z, alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = RFT3Dfunc(points, normals, area, rotation, angular_velocity, linear_velocity, direction_vector, rho_c, mu_int, mu_surf, gravity, unit_test);
    sim_data_dtos(step) = SimDataTransferObject(points, c_inc, v_norm_vec, F, f, forces_x, forces_y, forces_z, forces, T, torque_x, torque_y, torque_z, alpha_gen, alpha_gen_n, alpha_gen_t, alpha, TRG, TRG_visual);
    z(step) = sim_data_dtos(step).torque_z;
    step = step + 1;
end

for_excel = [depths; z];

points = sim_data_dtos(end).points;
c_inc = sim_data_dtos(end).c_inc;
v_norm_vec = sim_data_dtos(end).v_norm_vec;
F = sim_data_dtos(end).F;
f = sim_data_dtos(end).f;
forces_x = sim_data_dtos(end).forces_x;
forces_y = sim_data_dtos(end).forces_y;
forces_z = sim_data_dtos(end).forces_z;
forces = sim_data_dtos(end).forces;
T = sim_data_dtos(end).T;
torque_x = sim_data_dtos(end).torque_x;
torque_y = sim_data_dtos(end).torque_y;
torque_z = sim_data_dtos(end).torque_z;
alpha_gen = sim_data_dtos(end).alpha_gen;
alpha_gen_n = sim_data_dtos(end).alpha_gen_n;
alpha_gen_t = sim_data_dtos(end).alpha_gen_t;
alpha = sim_data_dtos(end).alpha;
TRG = sim_data_dtos(end).TRG;
TRG_visual = sim_data_dtos(end).TRG_visual;

plots_publication


plot(depths, abs(z), 'LineWidth', 1.5);
end

legend(tips , "Location", "northwest");