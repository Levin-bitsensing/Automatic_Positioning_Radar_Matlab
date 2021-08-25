function Draw_Image_Grid(grid_x_range, grid_y_range, P, image_handle)

for iy = 1 : length(grid_y_range)
    if iy == 1
        grid_long_wld_pt0 = [ grid_x_range(1), grid_y_range(iy) ];
        grid_long_wld_pt1 = [ grid_x_range(end), grid_y_range(iy) ];
    else
        grid_long_wld_pt0 = [ grid_long_wld_pt0 ; grid_x_range(1), grid_y_range(iy) ];
        grid_long_wld_pt1 = [ grid_long_wld_pt1 ; grid_x_range(end), grid_y_range(iy) ];
    end
end

for ix = 1 : length(grid_x_range)
    if ix == 1
        grid_lat_wld_pt0 = [ grid_x_range(ix), grid_y_range(1) ];
        grid_lat_wld_pt1 = [ grid_x_range(ix), grid_y_range(end) ];
    else
        grid_lat_wld_pt0 = [ grid_lat_wld_pt0 ; grid_x_range(ix), grid_y_range(1) ];
        grid_lat_wld_pt1 = [ grid_lat_wld_pt1 ; grid_x_range(ix), grid_y_range(end) ];
    end
end

% Draw grid
nGrid_lat = length(grid_x_range);
for idx = 1 : nGrid_lat
    grid_lat_img_pt0(idx, :)  = wld2img([grid_lat_wld_pt0(idx, 1:2), 0]',  P);
    grid_lat_img_pt1(idx, :)  = wld2img([grid_lat_wld_pt1(idx, 1:2), 0]',  P);

    plot([grid_lat_img_pt0(idx, 1); grid_lat_img_pt1(idx,1)], ...
         [grid_lat_img_pt0(idx, 2); grid_lat_img_pt1(idx,2)], ...
         'color', [1-idx/nGrid_lat, idx/nGrid_lat, 0], 'LineWidth', 0.5, 'parent', image_handle.Parent);
end

nGrid_long = length(grid_y_range);
for idx = 1 : nGrid_long
    grid_long_img_pt0(idx, :) = wld2img([grid_long_wld_pt0(idx, 1:2), 0]', P);
    grid_long_img_pt1(idx, :) = wld2img([grid_long_wld_pt1(idx, 1:2), 0]', P);

    plot([grid_long_img_pt0(idx, 1); grid_long_img_pt1(idx,1)], ...
         [grid_long_img_pt0(idx, 2); grid_long_img_pt1(idx,2)], ...
         'color', [1-idx/nGrid_long, idx/nGrid_long, 0], 'LineWidth', 0.5, 'parent', image_handle.Parent);
end