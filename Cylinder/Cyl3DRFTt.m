close all
clear all

%% NOTES - TODO
% Something wrong with distinction of force direction on cylinder
% 

%% Define inputs - Agarwal verification studies
linearVelocity = 0.1; % linear velocity in m/s
angularVelocity = pi; % angular velocity in rad/s
rhoC = 1310; % critical density of the sand in kg/m³   
muInt = 0.21; % internal friction coefficient of the sand
muSurf = 0.4; % intruder-surface interaction coefficient
depth = 0.125; % in m
showGeometry = false;
showDirectionV = false;
showFQuiver = false;
showFScatter = false;
showFScatterxyz = true;
showAlpha = false;
saveFigures = false;
unitTest = false;

%% Read .stl file
TRG = stlread('./Cylinder/Models/Cylinder.stl'); % possibilities: CylinderFine, Cylinder, CylinderRough, CylinderVeryRough
TRGVisual = stlread('./Cylinder/Models/CylinderVeryRough.stl');
% Rotate TRG object
TRG = rotateTriangulationX(TRG, 0);
TRGVisual = rotateTriangulationX(TRGVisual, 0);
% Move the triangulation object so that the lowest point regarding the z-axis is at z = depth
TRG = moveTriangulationZ(TRG, depth);
TRGVisual = moveTriangulationZ(TRGVisual, depth);
% Compute the incenter of each triangle in the tip mesh
points = (incenter(TRG)').';
% Compute the normal of each triangle in the tip mesh
normals = (faceNormal(TRG)').';
% Compute the area of each triangle in the tip mesh
area = (generateArea(TRG.Points', TRG.ConnectivityList')).';

% Compute the vector norm of each normal
norms = vecnorm(normals, 2, 2);
% Normalize the normals by dividing each normal by its vector norm
normals = normals ./ norms;

%% Compute forces using 3D-RFT function
[c_inc, vNormVec, F, f, forcesX, forcesY, forcesZ, T, torqueX, torqueY, torqueZ, alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = RFT3Dfunc(points, normals, area, angularVelocity, linearVelocity, rhoC, muInt, muSurf, unitTest);

%% Plots
CylPlots

function [c_inc, vNormVec, F, f, forcesX, forcesY, forcesZ, T, torqueX, torqueY, torqueZ, alpha_gen, alpha_gen_n, alpha_gen_t, alpha] = RFT3Dfunc(points, normals, area, angularVelocity, linearVelocity, rhoC, muInt, muSurf, unitTest)
%% 1. Read Tip Data
pointList = points;
areaList = area/1000000; % mm² to m²
normalList = normals;
depthList = pointList(:,3)/1000; % mm to m

%% 2. Calc velocity
nElements = size(pointList, 1);

% Direction Vector
vcor = [0; 0; -linearVelocity*1000];

rList = sqrt(pointList(:,1) .^ 2 + pointList(:,2) .^ 2);
angleList = atan2(pointList(:,2), pointList(:,1));

vx = sin(angleList) .* rList .* angularVelocity + vcor(1);
vy = -cos(angleList) .* rList .* angularVelocity + vcor(2);
vz = zeros(nElements, 1) + vcor(3);
vVec = [vx vy vz];
vNormVec = vVec ./ vecnorm(vVec, 2, 2);

% Surface normal
nNormVec = normalList ./ vecnorm(normalList, 2, 2);

%% 3. Check conditions
% logical conditions that determine whether or not to count force on element
is_leading_edge = dot(nNormVec, vNormVec, 2) > 0;
is_intruding = pointList(:,3) < 0;
include = is_leading_edge & is_intruding;

% isolate the elements that satisfy this condition for calculation -

n_inc = nNormVec(include,:);
v_inc = vNormVec(include,:);
a_inc = areaList(include,:);
c_inc = pointList(include,:);
d_inc = depthList(include,:);

if unitTest
n_inc = n_inc(212,:); % REMOVE LATER
v_inc = v_inc(212,:); % REMOVE LATER
a_inc = a_inc(212,:); % REMOVE LATER
c_inc = c_inc(212,:); % REMOVE LATER
d_inc = d_inc(212,:); % REMOVE LATER
end

%% 4. Find local coordinate frame
if unitTest
z_local = repmat([0,0,1], sum(1), 1); % JUST FOR TESTING
else
z_local = repmat([0,0,1], sum(include), 1); % reverse g directon (regular z)
end

r_local = (v_inc - dot(v_inc, z_local, 3) .* z_local) ./ vecnorm(v_inc - dot(v_inc, z_local, 3) .* z_local, 2, 2);
theta_local = cross(z_local, r_local, 2);
g_ = -z_local;

%% 5. Find RFT angles (beta, gamma, psi)
% beta - surface characteristic angle
% Initialize beta vector
beta = zeros(size(n_inc,1),1);
% Iterate over each row of the input matrices
for i = 1:size(n_inc,1)
% Check condition for each row
if (dot(n_inc(i,:),r_local(i,:), 2) >= 0) && (dot(n_inc(i,:),z_local(i,:), 2) >= 0)
    beta(i) = - acos(dot(n_inc(i,:),z_local(i,:), 2));
elseif  (dot(n_inc(i,:),r_local(i,:), 2) >= 0) && (dot(n_inc(i,:),z_local(i,:), 2) < 0)
    beta(i) = +pi - acos(dot(n_inc(i,:),z_local(i,:), 2));
elseif  (dot(n_inc(i,:),r_local(i,:), 2) < 0) && (dot(n_inc(i,:),z_local(i,:), 2) >= 0)
    beta(i) =     + acos(dot(n_inc(i,:),z_local(i,:), 2));
else 
    beta(i) = -pi + acos(dot(n_inc(i,:),z_local(i,:), 2));
end
end

% gamma - velocity characteristic angle
% Initialize gamma vector
gamma = zeros(size(v_inc,1),1);
% Iterate over each row of the input matrices
for i = 1:size(v_inc,1)
% Check condition for each row
if dot(v_inc(i,:), z_local(i,:), 2) <= 0
gamma(i) = acos(dot(v_inc(i,:), r_local(i,:), 2));
else
gamma(i) = -acos(dot(v_inc(i,:), r_local(i,:), 2));
end
end

% psi - surface characteristic angle
% Initialize psi vector
psi = zeros(size(n_inc,1),1);
% Iterate over each row of the input matrices
for i = 1:size(n_inc,1)
% Compute nr0_inc for each row
nr0_inc = (n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:))) ./ vecnorm(n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:)),2,2);
% Check condition for each row
if vecnorm(n_inc(i,:) - (dot(n_inc(i,:),z_local(i,:), 2) .* z_local(i,:)),2,2) == 0 || dot(nr0_inc,r_local(i,:),2) == 0
psi(i) = 0;
else
psi(i) = atan( dot(nr0_inc,theta_local(i,:), 2) ./ dot(nr0_inc,r_local(i,:), 2) );
end
end

%% 6a. Determine x1, x2, x3
% Local versions of n and v
% v_inc_local = cos(gamma) .* r_local - sin(gamma) .* z_local;
% n_inc_local = sin(beta) .* cos(psi) .* r_local + sin(beta) .* sin(psi) .* theta_local - cos(beta) .* z_local;

% y1 = dot(g_, v_inc_local, 2);
% y2 = dot(g_, n_inc_local, 2);
% y3 = dot(n_inc_local, v_inc_local, 2);

x1 = sin(gamma);
x2 = cos(beta);
x3 = cos(psi) .* cos(gamma) .* sin(beta) + sin(gamma) .* cos(beta);

if unitTest
unitx = ones(sum(1), 1); % JUST FOR TESTING
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
alpha_r_gen = f1 .* sin(beta) .* cos(psi) + f2 .* cos(gamma);
alpha_theta_gen = f1 .* sin(beta) .* sin(psi);
alpha_z_gen = -f1 .* cos(beta) - f2 .* sin(gamma) - f3;

alpha_gen = alpha_r_gen.*r_local + alpha_theta_gen.*theta_local + alpha_z_gen.*z_local;

%% 8. Estimate media specific scaling factor xi_n
xi_n = 0.12 * 10^6; % Agarwal verification studies
% xi_n = rhoC * 9.81 * (894*muInt^3 - 386*muInt^2 + 89*muInt); % initially in N/m³

%% 9. Calculate the system specific alpha_n and alpha_t in the local coordinate frame
alpha_gen_n = zeros(size(n_inc,1),3);
alpha_gen_t = zeros(size(n_inc,1),3);
for i = 1:size(n_inc,1)
    if dot(alpha_gen(i,:),-n_inc(i,:),2)<0
    alpha_gen_n(i,:) = -dot(alpha_gen(i,:),-n_inc(i,:),2) .* (-n_inc(i,:));
    alpha_gen_t(i,:) = (alpha_gen(i,:) + alpha_gen_n(i,:));
    else
    alpha_gen_n(i,:) = dot(alpha_gen(i,:),-n_inc(i,:),2) .* (-n_inc(i,:));
    alpha_gen_t(i,:) = (alpha_gen(i,:) - alpha_gen_n(i,:));
    end
end

% alpha_gen_n = dot(alpha_gen,-n_inc,2) .* (-n_inc);  % this has some mistake
% alpha_gen_t = (alpha_gen - alpha_gen_n);  % this might be correct

alpha = xi_n .* (alpha_gen_n + min(muSurf .* vecnorm(alpha_gen_n,2,2) ./ vecnorm(alpha_gen_t,2,2),1) .* alpha_gen_t);

%% 10. Calculate {alpha_x alpha_y alpha_z} as alpha
% Not needed when alpha_gen_t is calculated with n_inc instead of n_inc_local
% alpha_global = zeros(size(alpha,1),3); % pre-allocate rotation matrix
% 
% for i = 1:size(alpha,1)
%     % compute rotation axis and angle
%     rotvec = vrrotvec(n_inc_local(i,:), n_inc(i,:));
%     R = vrrotvec2mat(rotvec);
%     alpha_global(i,:) = alpha(i,:) * R;
% end

%% 11. multiplying up .* alpha * depth * area
F = alpha .* abs(d_inc) .* a_inc; % N
f = F ./ a_inc ./ 1000000; % N/mm²
T = cross(F,c_inc,2) ./ 1000; % Nmm to Nm

%% 12. sum all rows discrete intregral
[forcesX] = sum(F(:,1),1);
[forcesY] = sum(F(:,2),1);
[forcesZ] = sum(F(:,3),1);

[torqueX] = sum(T(:,1),1);
[torqueY] = sum(T(:,2),1);
[torqueZ] = sum(T(:,3),1);

end

function areaarray = generateArea(Points,List)
    % Initialize the area array with zeros
    areaarray = zeros(1, size(List, 2));

    % Loop over the triangles
    for i = 1:size(List, 2)  
        % Extract the coordinates of the triangle vertices
        v1 = Points(:, List(1,i));
        v2 = Points(:, List(2,i));
        v3 = Points(:, List(3,i));

        % Compute the lengths of the triangle sides
        a = norm(v1 - v2);
        b = norm(v2 - v3);
        c = norm(v1 - v3);

        % Use Heron's formula to compute the area of the triangle
        s = (a + b + c) / 2;
        AreaTemp = sqrt(s * (s - a) * (s - b) * (s - c));

        % Append the area to the area array
        areaarray(i) = AreaTemp;
    end
end

function TRG = rotateTriangulationX(TRG, theta)
% ROTATETRIANGULATIONX Rotates a triangulation object around the x-axis.
%   TRG = rotateTriangulationX(TRG, THETA) rotates the triangulation object TRG
%   around the x-axis by an angle THETA (in degrees). TRG is a triangulation
%   object containing points and a connectivity list. The function returns the
%   rotated triangulation object.

% Create the rotation matrix
R = [1 0 0; 0 cosd(theta) -sind(theta); 0 sind(theta) cosd(theta)];

% Rotate the points in the triangulation object
Points = TRG.Points * R;

% Create a new triangulation object with the rotated points and the same connectivity list
TRG = triangulation(TRG.ConnectivityList, Points);
end

function TRG = moveTriangulationZ(TRG, depth)

% MOVETRIANGULATIONZ Moves a triangulation object so that the lowest point regarding the z-axis is at z=0.
%   TRG = MOVETRIANGULATIONZ(TRG) moves the triangulation object TRG so that the
%   lowest point regarding the z-axis is at z=0. TRG is a triangulation object
%   containing points and a connectivity list. The function returns the moved
%   triangulation object.

% Find the minimum z-coordinate
minZ = min(TRG.Points(:, 3));

% Shift the points in the z-direction
Points = TRG.Points;
Points(:, 3) = Points(:, 3) - minZ - depth*1000;

% Create a new triangulation object with the shifted points and the same connectivity list
TRG = triangulation(TRG.ConnectivityList, Points);
end