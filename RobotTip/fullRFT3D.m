close all
clear all

%% Define inputs
linearVelocity = 0.01; % linear velocity in m/s
angularVelocity = 2 * pi / 5; % angular velocity in rad/s
rhoC = 1630; % critical density of the sand in kg/m³
thetaInt = 30; % internal friction angle of the sand in °
thetaSurf = 25; % intruder-surface friction angle in °
muInt = tan(thetaInt*pi/180); % internal friction coefficient of the sand
muSurf = tan(thetaSurf*pi/180); % intruder-surface interaction coefficient
depth = 0.2; % in m
showVisualisation = false;
showFVisualisation = true;

%% Read .stl file
TRG = stlread('TipNr8Rough.stl');
TRGVisual = stlread('TipNr8VeryRough.stl');
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

% Show the tip mesh
if showVisualisation
    figure
    hold on
    title ('Tip mesh visualization');
    colormap winter;
    view([45 25])
    trisurf(TRG);
    daspect([1 1 1]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis on;
    hold off;

    % Create a quiver plot with the normals as the arrow directions
    figure
    quiver3(points(:,1), points(:,2), points(:,3), normals(:,1), normals(:,2), normals(:,3), 1.25);
    title ('Normals of subsurfaces');
    colormap summer;
    view([45 25])
    daspect([1 1 1]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on
    axis on
    hold off;

    % Show the points in the mesh
    figure
    title ('Center points of subsurfaces');
    view([45 25])
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    grid on;
    zlim([-inf inf]);
    hold on;
    colormap winter;
    scatter3(points(:,1), points(:,2), points(:,3), 5, 'filled');
    hold off;

else 
end

%% Compute forces using 3D-RFT function
[c_inc, F, f, forces, T, torque] = RFT3Dfunc(points, normals, area, angularVelocity, linearVelocity, rhoC, muInt, muSurf, showVisualisation);

%% Plots
if showFVisualisation
    % Plot forces on each point of the mesh (quiver)
    figure;
    title ('Forces acting on each subsurface (quiver)');
    view([45 25])
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), F(:,1), F(:,2), F(:,3),2, 'LineWidth', 1, 'MaxHeadSize', 0);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);

    % Plot forces on each point of the mesh (scatter)
    figure;
    title ('Forces acting on each subsurface (scatter)');
    view([45 25])
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', vecnorm(f, 2, 2), 'SizeData', 300*abs(vecnorm(f, 2, 2))); %'SizeData', 300*abs(vecnorm(f, 2, 2))
    colormap(jet);
    clim([min(vecnorm(f, 2, 2)) max(vecnorm(f, 2, 2))]);
    colorbar;

    % Plot forces on each point of the mesh (scatter x)
    figure;
    title ('Forces (x) acting on each subsurface');
    view([45 25])
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,1), 'SizeData', 500*abs(f(:,1)));
    colormap(brewermap([],"-RdYlBu"));
    clim([min(f(:,1)) max(f(:,1))]);
    colorbar;

    % Plot forces on each point of the mesh (scatter y)
    figure;
    title ('Forces (y) acting on each subsurface');
    view([45 25])
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,2), 'SizeData', 500*abs(f(:,2)));
    colormap(brewermap([],"-RdYlBu"));
    clim([min(f(:,2)) max(f(:,2))]);
    colorbar;

    % Plot forces on each point of the mesh (scatter z)
    figure;
    title ('Forces (z) acting on each subsurface');
    view([45 25])
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,3), 'SizeData', 500*abs(f(:,3)));
    colormap(brewermap([],"YlOrRd"));
    clim([min(f(:,3)) max(f(:,3))]);
    colorbar;

else
end

function [c_inc, F, f, forces, T, torque] = RFT3Dfunc(points, normals, area, angularVelocity, linearVelocity, rhoC, muInt, muSurf, showVisualisation)
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

if showVisualisation
    % Create a quiver plot with the direction vectors
    figure
    quiver3(points(:,1), points(:,2), points(:,3), vNormVec(:,1), vNormVec(:,2), vNormVec(:,3), 1.25);
    title ('Direction vectors of rotation and translation');
    colormap summer;
    view([45 25])
    daspect([1 1 1]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    axis on;
    hold off;
else
end

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

%% 4. Find local coordinate frame
z_local = repmat([0,0,1], sum(include), 1); % reverse g directon (regular z)
g_ = -z_local;
r_local = (v_inc - dot(v_inc, z_local, 3) .* z_local) ./ vecnorm(v_inc - dot(v_inc, z_local, 3) .* z_local, 2, 2);
theta_local = cross(z_local, r_local, 2);

%% 5. Find RFT angles (beta, gamma, psi)
% beta - surface characteristic angle
% Initialize beta vector
beta = zeros(size(n_inc,1),1);
% Iterate over each row of the input matrices
for i = 1:size(n_inc,1)
% Check condition for each row
if (dot(n_inc(i,:),r_local(i,:), 2) >= 0) & (dot(n_inc(i,:),z_local(i,:), 2) >= 0)
    beta(i) = - acos(dot(n_inc(i,:),z_local(i,:), 2));
elseif  (dot(n_inc(i,:),r_local(i,:), 2) >= 0) & (dot(n_inc(i,:),z_local(i,:), 2) < 0)
    beta(i) = +pi - acos(dot(n_inc(i,:),z_local(i,:), 2));
elseif  (dot(n_inc(i,:),r_local(i,:), 2) < 0) & (dot(n_inc(i,:),z_local(i,:), 2) >= 0)
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
psi(i) = atan( dot(nr0_inc,r_local(i,:), 2) ./ dot(nr0_inc,r_local(i,:), 2) );
end
end

%% 6a. Determine x1, x2, x3
% Local versions of n and v
v_inc_local = cos(gamma) .* r_local - sin(gamma) .* z_local;
n_inc_local = sin(beta) .* cos(psi) .* r_local + sin(beta) .* sin(psi) .* theta_local - cos(beta) .* z_local;

y1 = dot(g_, v_inc_local, 2);
y2 = dot(g_, n_inc_local, 2);
y3 = dot(n_inc_local, v_inc_local, 2);

x1 = sin(gamma);
x2 = cos(beta);
x3 = cos(psi) .* cos(gamma) .* sin(beta) + sin(gamma) .* cos(beta);
unitx = ones(sum(include), 1);

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
alpha_gen_n = (dot(alpha_gen,n_inc,2)) .* (-n_inc);
alpha_gen_t = alpha_gen - alpha_gen_n;

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
[forces] = sum(F,1);
[torque] = sum(T,1);

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

function SetQuiverColor(q,currentColormap,varargin)
%--------------------------------------------------
% function SetQuiverColor(q,currentColormap)
%
% INPUT:
%   q = handle to quiver plot
%   currentColormap = e.g. jet;
% OPTIONAL INPUT ('Field',value):
%   'range' = [min,max]; % Range of the magnitude in the colorbar
%                          (used to possibly saturate or expand the color used compared to the vectors)
%   'mags' = magnitude; % Actual magnitude of the vectors
%
% Example:
%   [x,y] = meshgrid(-2:.2:2,-1:.15:1);
%   z = x .* exp(-x.^2 - y.^2);
%   [u,v,w] = surfnorm(x,y,z);
%   q = quiver3(x,y,z,u,v,w);
%   mag = 1+3.*rand(size(u));   % Creates number between 1 and 4
%   colormap(jet);
%   colorbar;
%   SetQuiverColor(q,jet,'mags',mag,'range',[-2 8]);  % Color range between -2 8 => all colors are not used
%   caxis([-2 8]);
%   set(gca,'Color','k');
%
%--------------------------------------------------
%   Authorship:
%     This code is heavily based from the answer by the user Suever on Stackoverflow forum
%     at: https://stackoverflow.com/questions/29632430/quiver3-arrow-color-corresponding-to-magnitude
%
%     I, Alexandre De Spiegeleer, only added minor changes to the original answer to have a bit more flexibility.
%--------------------------------------------------

%// Set default values
range = [];
mags = [];

%// Read the optional range value
if find(strcmp('range',varargin))
  range = varargin{ find(strcmp('range',varargin))+1 };
end

qU = q.UData(~isnan(q.UData));
qV = q.VData(~isnan(q.VData));
qW = q.WData(~isnan(q.WData));

%// Compute/read the magnitude of the vectors
if find(strcmp('mags',varargin))
  mags = varargin{ find(strcmp('mags',varargin))+1 };
  mags = mags(~isnan(mags)&~isnan(q.UData));  % This reshapes automatically
else
  mags = sqrt(sum(cat(2, qU, qV, ...
             reshape(qW, numel(qU), [])).^2, 2));
end
%// If range is auto, take range as the min and max of mags
if isstr(range) & strcmp(range,'auto')
  range = [min(mags) max(mags)];
end

%// Change value depending on the desired range
if ~isempty(range) & isnumeric(range) & numel(range)==2
  range = sort(range);
  mags(mags>range(2)) = range(2);
  mags(mags<range(1)) = range(1);
end

%// Now determine the color to make each arrow using a colormap
if ~isempty(range) & isnumeric(range) & numel(range)==2
  Edges = linspace(range(1),range(2),size(currentColormap, 1)+1);
  [~, ~, ind] = histcounts(mags, Edges);
else
  [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
end

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// Color data
cd_head = reshape(cmap(1:3,:,:), [], 4).';
cd_tail = reshape(cmap(1:2,:,:), [], 4).';

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, 'ColorBinding', 'interpolated', 'ColorData', cd_head);

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, 'ColorBinding', 'interpolated', 'ColorData', cd_tail);
end