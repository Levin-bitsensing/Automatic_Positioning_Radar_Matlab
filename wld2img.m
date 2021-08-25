function [ pt_img ] = wld2img( pt_wld, P )
    pt_img = P * [pt_wld ; 1];
    pt_img = pt_img(1:2) ./ pt_img(3);
end