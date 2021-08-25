function [ pt_wld ] = img2wld( pt_img, H )
    pt_wld = H * [ pt_img ; 1 ];
    pt_wld(1:2) = pt_wld(1:2) ./ pt_wld(3);
    pt_wld(3) = 0;
end