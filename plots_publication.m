% if show_geometry || show_direction || show_f_quiver || show_f_scatter || show_f_scatterxyz || show_alpha
%     figure
%     colormap jet;
%     view([45 25])
%     daspect([1 1 1]);
%     xlabel('X  [mm]');
%     ylabel('Y  [mm]');
%     zlabel('Z  [mm]');
%     axis on;
%     zlim([-inf inf]);


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
    quiver3(points(:,1), points(:,2), points(:,3), movement(:,1), movement(:,2), movement(:,3), 1.25);
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
    q = quiver3(points_include(:,1), points_include(:,2), points_include(:,3), -forces(:,1), -forces(:,2), -forces(:,3),2, 'LineWidth', 1, 'MaxHeadSize', 0);
    currentColormap = jet;
    setQuiverColor(q,currentColormap);
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
    scatter3(points_include(:,1), points_include(:,2), points_include(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', vecnorm(pressures, 2, 2), 'SizeData', 250*abs(vecnorm(pressures, 2, 2)));
    colormap(jet);
    clim([min(vecnorm(pressures, 2, 2)) max(vecnorm(pressures, 2, 2))]);
    set(findall(gcf,'-property','FontSize'),'FontSize',14);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_', object, triangle_size_calculation, '.png'));
    end
    hold off;
end


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
    scatter3(points_include(:,1), points_include(:,2), points_include(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', pressures(:,1), 'SizeData', 250*max(max(abs(pressures)))/max(abs(pressures(:,2)))*abs(pressures(:,1))); % -f just because of inverted colorbar
    colormap(jet);
    clim([-max(max(abs(pressures))) max(max(abs(pressures)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(pressures))) 0 max(max(abs(pressures)))];
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
    scatter3(points_include(:,1), points_include(:,2), points_include(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', pressures(:,2), 'SizeData', 250*max(max(abs(pressures)))/max(abs(pressures(:,2)))*abs(pressures(:,2)));
    colormap(jet);
    clim([-max(max(abs(pressures))) max(max(abs(pressures)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(pressures))) 0 max(max(abs(pressures)))];
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
    scatter3(points_include(:,1), points_include(:,2), points_include(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', pressures(:,3), 'SizeData', 250*max(max(abs(pressures)))/max(abs(pressures(:,2)))*abs(pressures(:,3)));
    colormap(jet);
    clim([-max(max(abs(pressures))) max(max(abs(pressures)))]);
    c = colorbar;
    c.Ticks = [-max(max(abs(pressures))) 0 max(max(abs(pressures)))];
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
    q = quiver3(points_include(:,1), points_include(:,2), points_include(:,3), alpha_gen(:,1), alpha_gen(:,2), alpha_gen(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    setQuiverColor(q,currentColormap);
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
    q = quiver3(points_include(:,1), points_include(:,2), points_include(:,3), alpha_gen_n(:,1), alpha_gen_n(:,2), alpha_gen_n(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    setQuiverColor(q,currentColormap);
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
    q = quiver3(points_include(:,1), points_include(:,2), points_include(:,3), alpha_gen_t(:,1), alpha_gen_t(:,2), alpha_gen_t(:,3),2, 'LineWidth', 2, 'ShowArrowHead','on', 'MaxHeadSize', 5);
    currentColormap = jet;
    setQuiverColor(q,currentColormap);
    if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpdf', '-r300', '-vector', strcat('./', folder, '/Figures/alpha_gen_t_quiver_', object, triangle_size_calculation, '.pdf'));
    end
    hold off;

end


if show_depth_dependend
    figure
    plot(depths, abs(result_depth_dependend), 'LineWidth', 1.5);
end