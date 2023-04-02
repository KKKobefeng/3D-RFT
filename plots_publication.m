close all
% Show the tip mesh
if show_geometry
    figure
    hold on
    colormap winter;
    view([45 25])
    trisurf(TRG);
    daspect([1 1 1]);
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    axis on;
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/visual_mesh_', object, triangle_size_calculation, '.png'));
    end
    hold off;

    % Create a quiver plot with the normals as the arrow directions
    figure
    quiver3(points(:,1), points(:,2), points(:,3), normals(:,1), normals(:,2), normals(:,3), 1.25);
    colormap summer;
    view([45 25])
    daspect([1 1 1]);
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    axis on;
    grid off;
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/visual_normals_', object, triangle_size_calculation, '.png'));
    end
    hold off;

    % Show the points in the mesh
    figure
    view([45 25])
    daspect([1 1 1]);
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    axis on;
    zlim([-inf inf]);
    hold on;
    colormap winter;
    scatter3(points(:,1), points(:,2), points(:,3), 5, 'filled');
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/visual_points_', object, triangle_size_calculation, '.png'));
    end
    hold off;
end

if show_direction
    % Create a quiver plot with the direction vectors
    figure
    quiver3(points(:,1), points(:,2), points(:,3), v_norm_vec(:,1), v_norm_vec(:,2), v_norm_vec(:,3), 1.25);
    title ('Direction vectors of rotation and translation');
    colormap summer;
    view([45 25])
    daspect([1 1 1]);
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    axis on;
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/visual_direction_vectors_', object, triangle_size_calculation, '.png'));
    end
    hold off;
end

if show_f_quiver
    % Plot forces on each point of the mesh (quiver)
    figure;
    view([45 25])
    daspect([1 1 1]);
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    zlim([-inf inf]);
    hold on;
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), -F(:,1), -F(:,2), -F(:,3),2, 'LineWidth', 1, 'MaxHeadSize', 0);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/forces_quiver_', object, triangle_size_calculation, '.png'));
    end
    hold off;
end

if show_f_scatter
    % Plot forces on each point of the mesh (scatter)
    figure;
    view([45 25])
    daspect([1 1 1]);
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    zlim([-inf inf]);
    hold on;
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', vecnorm(f, 2, 2), 'SizeData', 250*abs(vecnorm(f, 2, 2)));
    colormap(jet);
    clim([min(vecnorm(f, 2, 2)) max(vecnorm(f, 2, 2))]);
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_', object, triangle_size_calculation, '.png'));
    end
    hold off;
end


%     % Plot forces on each point of the mesh (scatter x)
%     figure;
%     %title ('f_X [N/m^2]');
%     view([45 25])
%     %view([-45 -10])  % bunny drill
%     daspect([1 1 1]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     zlim([-inf inf]);
%     grid off;
%     hold on;
%     trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
%     %trisurf(TRG)
%     scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,1), 'SizeData', 500*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,1))); % -f just because of inverted colorbar
%     colormap(jet);
%     %colormap(brewermap([],"-RdYlBu"));
%     clim([-max(max(abs(f))) max(max(abs(f)))]);
%     c = colorbar;
%     c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
%     c.Label.String = '\bf{f_X [N/m²]}';
%     c.Location = 'northoutside';
%     if saveFigures
%     set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_x_', object, triangle_size_calculation, '.pdf'));
%     end
%     hold off;
% 
%     % Plot forces on each point of the mesh (scatter y)
%     figure;
%     %title ('f_Y [N/m^2]');
%     view([45 25])
%     %view([-45 -10])  % bunny drill
%     daspect([1 1 1]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     zlim([-inf inf]);
%     grid off;
%     hold on;
%     trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
%     %trisurf(TRG)
%     scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,2), 'SizeData', 500*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,2)));
%     colormap(jet);
%     %colormap(brewermap([],"-RdYlBu"));
%     clim([-max(max(abs(f))) max(max(abs(f)))]);
%     c = colorbar;
%     c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
%     c.Label.String = '\bf{f_Y [N/m²]}';
%     c.Location = 'northoutside';
%     if saveFigures
%     set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_y_', object, triangle_size_calculation, '.pdf'));
%     end
%     hold off;
% 
%     % Plot forces on each point of the mesh (scatter z)
%     figure;
%     %title ('f_Z [N/m^2]');
%     view([45 25])
%     %view([-45 -10])  % bunny drill
%     daspect([1 1 1]);
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     zlim([-inf inf]);
%     grid off;
%     hold on;
%     trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#999999', 'FaceAlpha', 0);
%     %trisurf(TRG)
%     scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,3), 'SizeData', 500*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,3)));
%     colormap(jet);
%     %colormap(brewermap([],"YlOrRd"));
%     clim([-max(max(abs(f))) max(max(abs(f)))]);
%     c = colorbar;
%     c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
%     c.Label.String = '\bf{f_Z [N/m²]}';
%     c.Location = 'northoutside';
%     if saveFigures
%     set(gcf,'PaperPositionMode','auto')
%     print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_z_', object, triangle_size_calculation, '.pdf'));
%     end
%     hold off;


if show_f_scatterxyz
    % Multiplot
    % Define figure size and aspect ratio
    aspect_ratio = [1, 1, 1]; % aspect ratio for each subplot
    % Create figure
    fig_width = 1500;
    fig_height = 500;
    figure('Units', 'pixels', 'Position', [0, 0, fig_width, fig_height]);
    % Define common plot settings
    view_angle = [45, 25];
    daspect([1 1 1]);
    zlim([-inf inf]);
    grid off;
    % Plot forces on each point of the mesh (scatter x)
    subplot(1, 3, 1);
    hold on;
    title('\bf{f_X [N/mm²]}');
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,1), 'SizeData', 250*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,1))); % -f just because of inverted colorbar
    colormap(jet);
    clim([-max(max(abs(f))) max(max(abs(f)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
    c.Location = 'southoutside';
    view(view_angle);
    set(gca, 'DataAspectRatio', aspect_ratio);
    hold off;
    % Plot forces on each point of the mesh (scatter y)
    subplot(1, 3, 2);
    hold on;
    title('\bf{f_Y [N/mm²]}');
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,2), 'SizeData', 250*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,2)));
    colormap(jet);
    clim([-max(max(abs(f))) max(max(abs(f)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
    c.Location = 'southoutside';
    view(view_angle);
    set(gca, 'DataAspectRatio', aspect_ratio);
    hold off;
    % Plot forces on each point of the mesh (scatter z)
    subplot(1, 3, 3);
    hold on;
    title('\bf{f_Z [N/mm²]}');
    xlabel('X  [mm]');
    ylabel('Y  [mm]');
    zlabel('Z  [mm]');
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#999999', 'FaceAlpha', 0);
    scatter3(c_inc(:,1), c_inc(:,2), c_inc(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', f(:,3), 'SizeData', 250*max(max(abs(f)))/max(abs(f(:,2)))*abs(f(:,3)));
    colormap(jet);
    clim([-max(max(abs(f))) max(max(abs(f)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(f))) 0 max(max(abs(f)))];
    c.Location = 'southoutside';
    view(view_angle);
    set(gca, 'DataAspectRatio', aspect_ratio);
    set(findall(gcf, 'Type', 'Text'), 'FontSize', 14);
    set(findall(gcf, 'Type', 'Colorbar'), 'FontSize', 14);
    hold off;
    % Save figure
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_xyz_', object, triangle_size_calculation, '.png'));
    end
end

if show_alpha
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
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), alpha_gen(:,1), alpha_gen(:,2), alpha_gen(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/alpha_gen_quiver_', object, triangle_size_calculation, '.pdf'));
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
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), alpha_gen_n(:,1), alpha_gen_n(:,2), alpha_gen_n(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/alpha_gen_n_quiver_', object, triangle_size_calculation, '.pdf'));
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
    trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
    %trisurf(TRG)
    q = quiver3(c_inc(:,1), c_inc(:,2), c_inc(:,3), alpha_gen_t(:,1), alpha_gen_t(:,2), alpha_gen_t(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    SetQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/alpha_gen_t_quiver_', object, triangle_size_calculation, '.pdf'));
    end
    hold off;

end

if show_linear_f
    % Plot linear graph forces - depth
    x_depth = [0.1 0.1 0.1];
    y_forces = [forces_x forces_y forces_z];
    y_forces_agarwal = [0 0 42.55];  % pi --> 42.55 ; 0.5pi --> 53.80 ; 0.25pi --> 59.35
    figure;
    for i=1:3
        X = [0 x_depth(i)];
        Y = [0 y_forces_agarwal(i)];
        plot(X, Y, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--') % plot agarwal lines
        hold on
    end
    for i=1:3
        X = [0 x_depth(i)];
        Y = [0 y_forces(i)];
        plot(X, Y, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-')  % plot lines
        hold on
    end
    xlabel('Depth [m]');
    ylabel('Forces [N]');
    grid on;
    xlim([0 0.1]);
    ylim([-10 75]);
    legend('Reference study', '', '', 'Implementation', 'Location', 'northwest')
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/forces_depth_', object, '.pdf'));
    end
    hold off;


    % Plot linear graph torque - depth
    x_depth = [0.1 0.1 0.1];
    y_torque = [torque_x torque_y torque_z];
    y_torque_agarwal = [0 0 0.831];  % pi --> 0.831 ; 0.5pi --> 0.644 ; 0.25pi --> 0.375
    figure;
    for i=1:3
        X = [0 x_depth(i)];
        Y = [0 y_torque_agarwal(i)];
        plot(X, Y, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--') % plot agarwal lines
        hold on
    end
    for i=1:3
        X = [0 x_depth(i)];
        Y = [0 y_torque(i)];
        plot(X, Y, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-')  % plot lines
        hold on
    end
    xlabel('Depth [m]');
    ylabel('Torque [Nm]');
    xlim([0 0.1]);
    ylim([-0.1 1.6]);
    grid on;
    legend('Reference study', '', '', 'Implementation', 'Location', 'northwest')
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/torque_depth_', object, '.pdf'));
    end
    hold off;
end


function SetQuiverColor(q,currentColormap,varargin)
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


