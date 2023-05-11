% Plot forces on each point of the mesh (scatter)
figure("Position", [200 200 700 600]);
% title ('Pressures on each subsurface');
hold on;
trimesh(TRG_visual, 'LineWidth', 0.5, 'EdgeColor', '#888888', 'FaceAlpha', 0);
scatter3(points_include(:,1), points_include(:,2), points_include(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', vecnorm(pressures, 2, 2), 'SizeData', 100*abs(vecnorm(pressures, 2, 2))); %250*abs(vecnorm(pressures, 2, 2))
colormap(jet);

c = colorbar;
c.Location = "east";
c.Label.String = 'Pressure [N/mm²]';
clim([min(vecnorm(pressures, 2, 2)) max(vecnorm(pressures, 2, 2))]);
setPlotProperties(x_range, y_range, z_range);

set(findall(gcf,'-property','FontSize'),'FontSize',20);
if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_', object, triangle_size_calculation, '.png'));
end