function [c_inc, F, f, forces_x, forces_y, forces_z, forces, T, torque_x, torque_y, torque_z, alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = RFT3Dfunc(points, normals, areas, depth_list, rho_c, mu_int, mu_surf, gravity, unit_test, threshold, v_norm_vec, v_vec)

%% 2. Calc velocity


%% 3. Check conditions
is_leading_edge = dot(normals, v_norm_vec, 2) >= -threshold;
is_intruding = points(:,3) < 0;
include = is_leading_edge & is_intruding;

n_inc = normals(include,:);
v_inc = v_norm_vec(include,:);
v_nn_inc = v_vec(include,:);
a_inc = areas(include,:);
c_inc = points(include,:);
d_inc = depth_list(include,:);

if unit_test
pointsTest = [90 180];
n_inc = n_inc(pointsTest,:);
v_inc = v_inc(pointsTest,:);
a_inc = a_inc(pointsTest,:);
c_inc = c_inc(pointsTest,:);
d_inc = d_inc(pointsTest,:);
end

%% 4. Find local coordinate frame
if unit_test
z_local = repmat([0,0,1], numel(d_inc), 1);
else
z_local = repmat([0,0,1], sum(include), 1);
end

r_local = zeros(size(n_inc,1),3);
for i = 1:size(n_inc,1)
if (vecnorm(v_inc(i,:) - dot(v_inc(i,:), z_local(i,:), 3) .* z_local(i,:), 2, 2) == 0 && vecnorm(n_inc(i,:) - dot(n_inc(i,:), z_local(i,:), 3) .* z_local(i,:), 2, 2) == 0)
r_local(i,:) = [1; 0; 0];
elseif (vecnorm(v_nn_inc(i,:) - dot(v_nn_inc(i,:), z_local(i,:), 3) .* z_local(i,:), 2, 2) == 0 && vecnorm(n_inc(i,:) - dot(n_inc(i,:), z_local(i,:), 3) .* z_local(i,:), 2, 2) ~= 0)
r_local(i,:) = (n_inc(i,:) - dot(n_inc(i,:), z_local(i,:), 3) .* z_local(i,:)) ./ vecnorm(n_inc(i,:) - dot(n_inc(i,:), z_local(i,:), 3) .* z_local(i,:), 2, 2);
else
r_local(i,:) = (v_nn_inc(i,:) - dot(v_nn_inc(i,:), z_local(i,:), 3) .* z_local(i,:)) ./ vecnorm(v_nn_inc(i,:) - dot(v_nn_inc(i,:), z_local(i,:), 3) .* z_local(i,:), 2, 2);
end
end

theta_local = cross(z_local, r_local, 2);

%% 5. Find RFT angles (beta, gamma, psi)
% beta - surface characteristic angle
beta = zeros(size(n_inc,1),1);
for i = 1:size(n_inc,1)
if (dot(n_inc(i,:),r_local(i,:), 2) >= 0) && (dot(n_inc(i,:),z_local(i,:), 2) >= 0)
    beta(i) = - round(acos(dot(n_inc(i,:),z_local(i,:), 2)), 15);
elseif  (dot(n_inc(i,:),r_local(i,:), 2) >= 0) && (dot(n_inc(i,:),z_local(i,:), 2) < 0)
    beta(i) = +pi - round(acos(dot(n_inc(i,:),z_local(i,:), 2)), 15);
elseif  (dot(n_inc(i,:),r_local(i,:), 2) < 0) && (dot(n_inc(i,:),z_local(i,:), 2) >= 0)
    beta(i) =     + round(acos(dot(n_inc(i,:),z_local(i,:), 2)), 15);
else 
    beta(i) = -pi + round(acos(dot(n_inc(i,:),z_local(i,:), 2)), 15);
end
end

% gamma - velocity characteristic angle
gamma = zeros(size(v_inc,1),1);
for i = 1:size(v_inc,1)
if dot(v_inc(i,:), z_local(i,:), 2) <= 0
gamma(i) = round(acos(dot(v_inc(i,:), r_local(i,:), 2)), 15);
else
gamma(i) = -round(acos(dot(v_inc(i,:), r_local(i,:), 2)), 15);
end
end

% psi - surface characteristic angle
psi = zeros(size(n_inc,1),1);
nr0_inc = zeros(size(n_inc,1),3);
for i = 1:size(n_inc,1)
nr0_inc(i,:) = (n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:))) ./ vecnorm(n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:)),2,2);
if vecnorm(n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:)),2,2) == 0 || dot(nr0_inc(i,:),r_local(i,:),2) == 0
psi(i) = 0;
else
psi(i) = round(atan( dot(nr0_inc(i,:),theta_local(i,:), 2) ./ dot(nr0_inc(i,:),r_local(i,:), 2) ), 15);
end
end

%% 6a. Determine x1, x2, x3

% beta = acos(dot(-z_local, n_inc,2));
% gamma = asin(dot(-z_local, v_inc,2));
% n_dot_v = dot(n_inc, v_inc, 2);
% for i = 1:size(n_inc,1)
%     if (vecnorm(n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:)),2,2) == 0 || dot(nr0_inc(i,:),r_local(i,:),2) == 0)
%         psi(i) = 0;
%     else
%         psi(i) = atan2(sin(beta(i,:))*sin(gamma(i,:)), n_dot_v(i,:) - sin(gamma(i,:)) * cos(beta(i,:)));
%     end
% end

x1 = round(sin(gamma), 15);
x2 = round(cos(beta), 15);
x3 = round(cos(psi), 15) .* round(cos(gamma), 15) .* round(sin(beta), 15) + round(sin(gamma), 15) .* round(cos(beta), 15);

y1 = dot(-z_local, v_inc,2);
y2 = dot(-z_local, n_inc,2);
y3 = dot(n_inc, v_inc,2);

x1(abs(x1) < threshold) = 0;
x2(abs(x2) < threshold) = 0;
x3(abs(x3) < threshold) = 0;

if unit_test
unitx = ones(numel(d_inc), 1);
else
unitx = ones(sum(include), 1);
end

%% 6b. Determine f1, f2, f3
Tk = [unitx x1  x2  x3  x1.^2    x2.^2    x3.^2    x1.*x2  x2.*x3  x3.*x1  x1.^3    x2.^3    x3.^3    x1.*x2.^2    x2.*x1.^2    x2.*x3.^2    x3.*x2.^2    x3.*x1.^2  x1.*x3.^2  x1.*x2.*x3];

c1k = [0.00212; -0.02320; -0.20890; -0.43083; -0.00259; 0.48872; -0.00415; 0.07204; -0.02750; -0.08772; 0.01992; -0.45961; 0.40799; -0.10107; -0.06576; 0.05664; -0.09269; 0.01892; 0.01033; 0.15120]; 
c2k = [-0.06796; -0.10941; 0.04725; -0.06914; -0.05835; -0.65880; -0.11985; -0.25739; -0.26834; 0.02692; -0.00736; 0.63758; 0.08997; 0.21069; 0.04748; 0.20406; 0.18589; 0.04934; 0.13527; -0.33207];
c3k = [-0.02634; -0.03436; 0.45256; 0.00835; 0.02553; -1.31290; -0.05532; 0.06790; -0.16404; 0.02287; 0.02927; 0.95406; -0.00131; -0.11028; 0.01487; -0.02730; 0.10911; -0.04097; 0.07881; -0.27519];

f1 = Tk * c1k;
f2 = Tk * c2k;
f3 = Tk * c3k;

%% 7. Calculate alpha_r_gen, alpha_theta_gen, alpha_z_gen
alpha_r_gen = f1 .* round(sin(beta), 15) .* round(cos(psi), 15) + f2 .* round(cos(gamma), 15);
alpha_theta_gen = f1 .* round(sin(beta), 15) .* round(sin(psi), 15);
alpha_z_gen = -f1 .* round(cos(beta), 15) - f2 .* round(sin(gamma), 15) - f3;

alpha_gen = alpha_r_gen.*r_local + alpha_theta_gen.*theta_local + alpha_z_gen.*z_local;

%% 8. Estimate media specific scaling factor xi_n
%xi_n = 0.082 * 10^6; % Agarwal verification studies - same values of f for xi_n = 0.08
xi_n = rho_c * gravity * (894*mu_int^3 - 386*mu_int^2 + 89*mu_int); % initially in N/m³

%% 9. Calculate the system specific alpha_n and alpha_t in the local coordinate frame
% Correcting minor sign problems
alpha_gen_n = zeros(size(n_inc,1),3);
alpha_gen_t = zeros(size(n_inc,1),3);
for i = 1:size(n_inc,1)
    if dot(alpha_gen(i,:),-n_inc(i,:),2) < 0
    alpha_gen_n(i,:) = -dot(alpha_gen(i,:),-n_inc(i,:),2) .* (-n_inc(i,:));
    alpha_gen_t(i,:) = (alpha_gen(i,:) + alpha_gen_n(i,:));
    else
    alpha_gen_n(i,:) = dot(alpha_gen(i,:),-n_inc(i,:),2) .* (-n_inc(i,:));
    alpha_gen_t(i,:) = (alpha_gen(i,:) - alpha_gen_n(i,:));
    end
end

for i = 1:size(n_inc,1)
    if dot(alpha_gen_t(i,:),-v_inc(i,:),2) < 0
    alpha_gen_t(i,:) = -(alpha_gen_t(i,:));
    end
end

% original
% alpha_gen_n = dot(alpha_gen,n_inc,2) .* (-n_inc);
% alpha_gen_t = (alpha_gen - alpha_gen_n);
 
alpha = xi_n .* (alpha_gen_n + min(mu_surf .* vecnorm(alpha_gen_n,2,2) ./ vecnorm(alpha_gen_t,2,2),1) .* alpha_gen_t);

%% 10. Calculate {alpha_x alpha_y alpha_z} as alpha
% Not needed

%% 11. multiplying up .* alpha * depth * area
F = alpha .* abs(d_inc) .* a_inc; % N
f = F ./ a_inc ./ 1000000; % N/mm²

%% 12. sum all rows discrete intregral
[forces_x] = sum(F(:,1),1);
[forces_y] = sum(F(:,2),1);
[forces_z] = sum(F(:,3),1);
[forces] = sqrt(forces_x^2 + forces_y^2 + forces_z^2);

%% 13. torque calculation
if unit_test
c_null = zeros(numel(d_inc), 1);
else
c_null = zeros(sum(include),1);
end

c_x = [c_null c_inc(:,2) c_inc(:,3)];
c_y = [c_inc(:,1) c_null c_inc(:,3)];
c_z = [c_inc(:,1) c_inc(:,2) c_null];

F_x = [c_null F(:,2) F(:,3)];
F_y = [F(:,1) c_null F(:,3)];
F_z = [F(:,1) F(:,2) c_null];

T_x = cross(c_x,F_x) ./ 1000;
T_y = cross(c_y,F_y) ./ 1000;
T_z = cross(c_z,F_z) ./ 1000;

T = [T_x(:,1) T_y(:,2) T_z(:,3)];

[torque_x] = sum(T_x(:,1), 1);
[torque_y] = sum(T_y(:,2), 1);
[torque_z] = sum(T_z(:,3), 1);
end