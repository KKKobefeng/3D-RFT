% Plot forces on each point of the mesh (scatter)
figure;
hold on;
trimesh(TRG_visual, 'LineWidth', 0.1, 'EdgeColor', '#888888', 'FaceAlpha', 0);
scatter3(points_include(:,1), points_include(:,2), points_include(:,3), 'filled', 'MarkerEdgeColor', 'none', 'CData', vecnorm(pressures, 2, 2), 'SizeData', 250*abs(vecnorm(pressures, 2, 2)));
colormap(jet);
clim([min(vecnorm(pressures, 2, 2)) max(vecnorm(pressures, 2, 2))]);
setPlotProperties(x_range, y_range, z_range);
if saveFigures
    set(gcf,'PaperPositionMode','auto')
    print(gcf, '-dpng', '-r300', '-vector', strcat('./', folder, '/Figures/forces_scatter_', object, triangle_size_calculation, '.png'));
end