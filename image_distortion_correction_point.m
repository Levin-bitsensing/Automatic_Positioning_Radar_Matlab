function [img_pt_cor] = image_distortion_correction_point(img_pt, intrinsic, k, p)

    x = img_pt(1,:);
    y = img_pt(2,:);

    fx = intrinsic(1);
    fy = intrinsic(2);
    cx = intrinsic(3);
    cy = intrinsic(4);
    skew = intrinsic(5);

    nIter = 5;
    
    % normalize the points
    yn = (y - cy) / fy;
    xn = (x - cx - yn * skew) / fx;

    xn0 = xn;
    yn0 = yn;

    for i = 1 : nIter

        % compute radial distortion
        r2 = xn .^ 2 + yn .^ 2;
        r4 = r2 .* r2;
        r6 = r2 .* r4;
        alpha = k(1) * r2 + k(2) * r4 + k(3) * r6;

        % compute tangential distortion
        xyProduct = xn .* yn;
        dxTangential = 2 * p(1) * xyProduct + p(2) * (r2 + 2 * xn .^ 2);
        dyTangential = p(1) * (r2 + 2 * yn .^ 2) + 2 * p(2) * xyProduct;

        % apply the distortion to the points
        xnd = xn + xn .* alpha + dxTangential;
        ynd = yn + yn .* alpha + dyTangential;

        dx = xnd - xn0;
        dy = ynd - yn0;

        % undistortion
        xn = xn - dx;
        yn = yn - dy;
    end

    % convert back to pixels
    xi = xn * fx + cx + yn * skew;
    yi = yn * fy + cy;

    xi = round(xi);
    yi = round(yi);

    img_pt_cor = [xi; yi];
end