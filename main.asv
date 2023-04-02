clear
close

%% Intruder geometry
folder = 'cylinder';                                    % cylinder, simple, robottip
object = 'cylinder';                                    % name of stl
triangle_size_calculation = 'veryrough';                % 'Fine', 'Normal', 'Rough', 'VeryRough'
triangle_size_visualization = 'veryrough';              % 'Fine', 'Normal', 'Rough', 'VeryRough'

%% Physical Properties
rho_c = 1310;                                           % bulk density of the sand in kg/m³   
mu_int = 0.21;                                          % internal friction coefficient of the sand
mu_surf = 0.4;                                          % intruder-surface interaction coefficient
gravity = 9.81;                                         % gravity in m/s²

%% Movement parameters
rotation = true;                                        % true or false
linear_velocity = 0.1;                                  % linear velocity in m/s
direction_angle_xz = -90 * pi / 180;                    % angle between direction and x-z-axis
direction_angle_y = -90 * pi / 180;                     % angle between direction and y-axis
angular_velocity = [0, 0, -1*pi];                       % angular velocity in rad/s
direction_vector = [round(cos(direction_angle_xz), 15) ...
    round(cos(direction_angle_y), 15) round(sin(direction_angle_xz), 15)];

%% Depth parameters
start_depth = 0;
end_depth = 0.125;
step_size = 0.001;

%% Plot options
show_geometry = 0;
show_direction = 0;

show_f_quiver = 0;
show_alpha = 0;

show_f_scatter = 0;
show_f_scatterxyz = 0;

show_linear_f = 0;
show_depth_dependend = 1;

saveFigures = 0;

%% Miscellaneous
unit_test = false;
threshold = 1.0e-12;


%% Read .stl file
TRG = stlread(strcat('./', folder, '/Models/', object, ...
    triangle_size_calculation, '.stl'));                            % mesh for calculation
TRG_visual = stlread(strcat('./', folder, '/Models/', object, ...
    triangle_size_visualization, '.stl'));                          % mesh for force plots


TRG = rotateTriangulationX(TRG, -0);                              % rotate TRG object
TRG_visual = rotateTriangulationX(TRG_visual, -0);                % rotate Visual


%% Loop over depths
num_steps = (end_depth - start_depth)./step_size;
depths = start_depth:step_size:end_depth;
results = zeros(6, num_steps);

step = 1;

for depth = start_depth:step_size:end_depth
    TRG = moveTriangulationZ(TRG, depth);                         % align align bottom with depth
    TRG_visual = moveTriangulationZ(TRG_visual, depth);           % align bottom with depth (visual)
    
    % step 1: read STL data
    [points, normals, areas, depth_list] = readStlData(TRG, TRG.Points', TRG.ConnectivityList');

    % step 2: define movement vector
    [v_vec, v_norm_vec] = calcVelocity(points, direction_vector, linear_velocity, rotation, angular_velocity, threshold);



    
    [c_inc, F, f, forces_x, forces_y, forces_z, forces, T, torque_x, torque_y, torque_z, alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = RFT3Dfunc(points, normals, areas, depth_list, rho_c, mu_int, mu_surf, gravity, unit_test, threshold, v_norm_vec, v_vec);

    results(:,step) = [forces_x; forces_y; forces_z; torque_x; torque_y; torque_z];

    disp(['Depth processed: ', num2str(depth), 'm']);

    step = step + 1;
end


%% Plotting
result_depth_dependend = results(6,:);
plots_publication

