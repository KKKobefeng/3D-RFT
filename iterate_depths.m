% Iterate through depth 
clear
close

start_depth = 0;
end_depth = 0.06;
step_size = 0.005;
num_steps = (end_depth - start_depth)./step_size;

%% Plot options
show_geometry = false;
show_direction = false;

show_f_quiver = false;
show_alpha = false;

show_f_scatter = false;
show_f_scatterxyz = false;

show_linear_f = false;

saveFigures = false;


sim_data_dtos = SimDataTransferObject.empty(0,num_steps);

step = 1;
depths = start_depth:step_size:end_depth;
z = zeros(1, num_steps);

for cur_depth = start_depth:step_size:end_depth
    sim_data_dtos(step) = RFT3D(cur_depth);
    z(step) = sim_data_dtos(step).forces_x;
    step = step + 1;
end

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

RFTPlots

figure
plot(depths, abs(z), 'LineWidth', 1.5)