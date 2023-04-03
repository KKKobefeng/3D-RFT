clear
close all

%% Intruder geometry
folder = 'robottip';                                    % cylinder, simple, robottip
object = 'tipnr3';                                    % name of stl
triangle_size_calculation = 'normal';                    % 'Fine', 'Normal', 'Rough', 'VeryRough'
triangle_size_visualization = 'veryrough';              % 'Fine', 'Normal', 'Rough', 'VeryRough'

%% Physical Properties
rho_c = 1605;                                           % bulk density of the sand in kg/m³   
mu_int = 1.07;                                          % internal friction coefficient of the sand
mu_surf = 0.1806;                                          % intruder-surface interaction coefficient
gravity = 9.81;                                         % gravity in m/s²
xi_n = rho_c * gravity * (894*mu_int^3 - 386*mu_int^2 + 89*mu_int); % initially in N/m³

%% Movement parameters
rotation = true;                                        % true or false
linear_velocity = 0.1;                                  % linear velocity in m/s
direction_angle_xz = -90 * pi / 180;                    % angle between direction and x-z-axis
direction_angle_y = -90 * pi / 180;                     % angle between direction and y-axis
angular_velocity = [0, 0, -8*pi];                       % angular velocity in rad/s
direction_vector = [round(cos(direction_angle_xz), 15) ...
    round(cos(direction_angle_y), 15) round(sin(direction_angle_xz), 15)];

%% Depth parameters
start_depth = 0;
end_depth = 0.12;
step_size = 0.005;

%% Plot options
show_geometry = 1;
show_direction = 0;

show_f_quiver = 0;
show_alpha = 0;

show_f_scatter = 0;
show_f_scatterxyz = 0;

show_linear_f = 0;
show_depth_dependend = 0;

saveFigures = 0;

%% Miscellaneous
unit_test = false;
threshold = 1.0e-12;


%% Read .stl file
TRG = stlread(strcat('./', folder, '/Models/', object, ...
    triangle_size_calculation, '.stl'));                            % mesh for calculation
TRG_visual = stlread(strcat('./', folder, '/Models/', object, ...
    triangle_size_visualization, '.stl'));                          % mesh for force plots


TRG = rotateTriangulationX(TRG, -90);                                % rotate TRG object
TRG_visual = rotateTriangulationX(TRG_visual, -90);                  % rotate Visual


%% Loop over depths
num_steps = (end_depth - start_depth)./step_size;
depths = start_depth:step_size:end_depth;
results = zeros(7, num_steps);

step = 1;

for depth = start_depth:step_size:end_depth

    TRG = moveTriangulationZ(TRG, depth);                         % align align bottom with depth
    TRG_visual = moveTriangulationZ(TRG_visual, depth);           % align bottom with depth (visual)
    
    % step 1: process STL data
    [points, normals, areas, depth_list] = getStlData(TRG, TRG.Points', TRG.ConnectivityList');

    % step 2: define movement vector
    [movement, movement_normalized] = calcVelocity(points, direction_vector, linear_velocity, rotation, angular_velocity, threshold);

    % step 3: RFT conditions
    [include, normals_include, movement_normalized_include, movement_include, areas_include, points_include, depth_list_include] ...
        = checkConditions(points, normals, areas, movement, movement_normalized, depth_list, unit_test, threshold);

    % step 4: find local coordinate frame
    [z_local, r_local, theta_local] = findLocalFrame(normals_include, movement_normalized_include, movement_include, include, unit_test);

    % step 5: RFT characteristic angles
    [beta, gamma, psi] = findAngles(normals_include, movement_normalized_include, r_local, z_local, theta_local);

    % step 6: empirically determined force components
    [f1, f2, f3] = findFit(gamma, beta, psi, z_local, movement_normalized_include, normals_include, threshold, depth_list_include, include, unit_test);

    % step 7: find generic form of alpha, split into normal and tangential
    % components and add soil + interface properties for to get alpha
    [alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = findAlpha(normals_include, movement_normalized_include, beta, gamma, psi, r_local, theta_local, z_local, f1, f2, f3, mu_surf, xi_n);

    % step 8: get forces
    [forces, pressures, force_x, force_y, force_z, resultant] = getForces(depth_list_include, areas_include, alpha);

    % step 9: get torques
    [T, torque_x, torque_y, torque_z] = getTorques(points_include, depth_list_include, forces, unit_test, include);


    results(:,step) = [depth; force_x; force_y; force_z; torque_x; torque_y; torque_z];

    disp(['Depth processed: ', num2str(depth), 'm']);
    step = step + 1;

end


%% Plotting
result_depth_dependend = results(7,:);
plots_publication

%% Finish
varList = evalin('caller', 'who');
matches = regexp(varList, '^show_\w*', 'match');
matches = vertcat(matches{:});
clear(matches{:})

disp("Done!");