% Show the tip mesh
if showGeometry
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
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/visual_mesh.pdf');
    end
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
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/visual_normals.pdf');
    end
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
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/visual_points.pdf');
    end
    hold off;
end


if showDirectionV
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
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/visual_direction_vectors.pdf');
    end
    hold off;
end

if showFQuiver
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
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), -F(:,1), -F(:,2), -F(:,3),2, 'LineWidth', 1, 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/forces_quiver.pdf');
    end
    hold off;
end

if showFScatter
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
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', vecnorm(f, 2, 2), 'SizeData', 500*abs(vecnorm(f, 2, 2)));
    colormap(jet);
    clim([min(vecnorm(f, 2, 2)) max(vecnorm(f, 2, 2))]);
    colorbar;
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/forces_scatter.pdf');
    end
    hold off;
end

if showFScatterxyz
    % Plot forces on each point of the mesh (scatter x)
    figure;
    %title ('f_X [N/m^2]');
    view([45 25])
    %view([-45 -10])  % bunny drill
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,1), 'SizeData', 500*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,1))); % -f just because of inverted colorbar
    colormap(jet);
    %colormap(brewermap([],"-RdYlBu"));
    clim([-max(max(abs(f))) max(max(abs(f)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
    c.Label.String = '\bf{f_X [N/m²]}';
    c.Location = 'northoutside';
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/forces_scatter_x.pdf');
    end
    hold off;

    % Plot forces on each point of the mesh (scatter y)
    figure;
    %title ('f_Y [N/m^2]');
    view([45 25])
    %view([-45 -10])  % bunny drill
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,2), 'SizeData', 500*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,2)));
    colormap(jet);
    %colormap(brewermap([],"-RdYlBu"));
    clim([-max(max(abs(f))) max(max(abs(f)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
    c.Label.String = '\bf{f_Y [N/m²]}';
    c.Location = 'northoutside';
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/forces_scatter_y.pdf');
    end
    hold off;

    % Plot forces on each point of the mesh (scatter z)
    figure;
    %title ('f_Z [N/m^2]');
    view([45 25])
    %view([-45 -10])  % bunny drill
    daspect([1 1 1]);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    zlim([-inf inf]);
    grid on;
    hold on;
    trimesh(TRGVisual, 'LineWidth', 0.1, 'EdgeColor', '#999999', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,3), 'SizeData', 500*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,3)));
    colormap(jet);
    %colormap(brewermap([],"YlOrRd"));
    clim([-max(max(abs(f))) max(max(abs(f)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
    c.Label.String = '\bf{f_Z [N/m²]}';
    c.Location = 'northoutside';
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/forces_scatter_z.pdf');
    end
    hold off;
end

if showAlpha
    figure;
    title ('forces alpha_{gen} (quiver)');
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
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), alpha(:,1), alpha(:,2), alpha(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/alpha_n_quiver.pdf');
    end
    hold off;

    figure;
    title ('Normal forces alpha_{gen,n} (quiver)');
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
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), alpha_gen_n(:,1), alpha_gen_n(:,2), alpha_gen_n(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/alpha_n_quiver.pdf');
    end
    hold off;

    figure;
    title ('Tangential forces alpha_{gen,t} (quiver)');
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
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), alpha_gen_t(:,1), alpha_gen_t(:,2), alpha_gen_t(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', './Cylinder/Figures/alpha_t_quiver.pdf');
    end
    hold off;

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