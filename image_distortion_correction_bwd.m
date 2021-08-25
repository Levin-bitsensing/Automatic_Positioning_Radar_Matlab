function [dst_image] = image_distortion_correction_bwd(src_image, intrinsic, k, p)

    nx = size(src_image, 2);
    ny = size(src_image, 1);
    
    [X, Y] = meshgrid(1:nx, 1:ny);
    points = [X(:) Y(:)]; % remapmex requires singles

    fx = intrinsic(1);
    fy = intrinsic(2);
    cx = intrinsic(3);
    cy = intrinsic(4);
    skew = intrinsic(5);

    % center the points
    center = [cx, cy];
    centeredPoints = bsxfun(@minus, points, center);

    % normalize the points
    yNorm = centeredPoints(:, 2, :) ./ fy;
    xNorm = centeredPoints(:, 1, :) - yNorm * skew;
    xNorm = xNorm ./ fx;

    % compute radial distortion
    r2 = xNorm .^ 2 + yNorm .^ 2;
    r4 = r2 .* r2;
    r6 = r2 .* r4;

    alpha = k(1) * r2 + k(2) * r4 + k(3) * r6;

    % compute tangential distortion
    xyProduct = xNorm .* yNorm;
    dxTangential = 2 * p(1) * xyProduct + p(2) * (r2 + 2 * xNorm .^ 2);
    dyTangential = p(1) * (r2 + 2 * yNorm .^ 2) + 2 * p(2) * xyProduct;

    % apply the distortion to the points
    normalizedPoints = [xNorm, yNorm];
    distortedNormalizedPoints = normalizedPoints ...
                              + normalizedPoints .* [alpha, alpha]  ...
                              + [dxTangential, dyTangential];

    % convert back to pixels
    distortedPoints = [distortedNormalizedPoints(:, 1, :) * fx + cx + ...
                       distortedNormalizedPoints(:, 2, :) * skew, ...
                       distortedNormalizedPoints(:, 2, :) * fy + cy];

    % interpolate image point
    % pad = 0;
    pad = -1;
    XmapSingle = cast(reshape(distortedPoints(:,2),[ny nx]), 'single') + pad;
    YmapSingle = cast(reshape(distortedPoints(:,1),[ny nx]), 'single') + pad;

    method = 'bilinear';
    % method = 'nearest';    
    
    if 1
        dst_image = vision.internal.calibration.interp2d(src_image, YmapSingle, XmapSingle, method, 0);
    else
        MqR = interp2(single(src_image(:,:,1)),YmapSingle,XmapSingle,method);
        MqG = interp2(single(src_image(:,:,2)),YmapSingle,XmapSingle,method);
        MqB = interp2(single(src_image(:,:,3)),YmapSingle,XmapSingle,method);

        Mq = zeros(ny,nx,3); % initialize new image matrix

        Mq(:,:,1) = MqR;
        Mq(:,:,2) = MqG;
        Mq(:,:,3) = MqB;

        dst_image = uint8(Mq);
    end

    % matlab toolbox
    % undistortedImage = undistortImage(src_image, cameraParams);

end