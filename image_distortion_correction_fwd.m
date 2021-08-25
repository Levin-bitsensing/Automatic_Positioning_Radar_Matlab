function [dst_image] = image_distortion_correction_fwd(src_image, intrinsic, k, p)

    nx = size(src_image, 2);
    ny = size(src_image, 1);
    dst_image = zeros(ny, nx, 3, 'uint8');

    for y = 1 : ny
        for x = 1 : nx
            
            img_pt = [x; y];
            
            img_pt_cor = image_distortion_correction_point(img_pt, intrinsic, k, p);
            
            xi = img_pt_cor(1);
            yi = img_pt_cor(2);
            
            if (xi >= 1) && (xi <= nx) && (yi >= 1) && (yi <= ny)

                dst_image(yi, xi, :) = src_image(y, x, :);

            end
        end
    end
end
