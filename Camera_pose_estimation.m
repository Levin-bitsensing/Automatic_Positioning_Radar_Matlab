close all
clear
clc

%% Image input and parameters

colors = 'brgkcm';

SHOW_IMAGE = 1;
SHOW_IMAGE_PTS = 1;
SHOW_GRID = 1;
SHOW_GRID_2 = 0;
SHOW_3D_VIEW = 1;
HOMOGRAPHY_ESTIMATION = 0;

nx = 1920;
ny = 1080;

fx = 2119.137451;
fy = 2120.412109;
cx = 925.452271;
cy = 564.02832;
skew = -11.126279;
alpha_c = skew / fx;
intrinsic = [fx, fy, cx, cy, skew]; 

rx_basis = 180;
ry_basis = -90;
rz_basis = -90;

k = [-0.560545; 0.515465; -0.070978];
p = [-0.001732; -1.6e-05];

inst_tl = [-7.0; 0.0; 10.0];            % xc, yc, zc
inst_euld = -1 * [0.0; 9.8; 3.5];       % rx(roll), ry(pitch), rz(yaw)

%% Load data

[pts_img, input_file_name] = load_data();

for idx_cmr = 1 : size(input_file_name,2)
    input_image{idx_cmr} = imread([input_file_name{idx_cmr} '.jpg']);
end


%% Image distortion correction


for idx_cmr = 1 : size(input_image,2)
    
    input_image{idx_cmr} = imread([input_file_name{idx_cmr} '.jpg']);
    
    input_image_undistorted{idx_cmr} = image_distortion_correction_bwd(input_image{idx_cmr}, intrinsic, k, p);
%     input_image_undistorted{idx_cmr} = image_distortion_correction_fwd(input_image{idx_cmr}, intrinsic, k, p);
    
end


%% Image point distortion correction

for idx_cmr = 1 : size(pts_img,2)
    pts_img_cor{idx_cmr} = image_distortion_correction_point(pts_img{idx_cmr}, intrinsic, k, p);
        
end






%% Preparing perspective projection

K = [fx, skew, cx; 0, fy, cy; 0, 0, 1];

R_tf = rotz(rz_basis) * roty(ry_basis) * rotx(rx_basis);

R_inst = rotz(inst_euld(3)) * roty(inst_euld(2)) * rotx(inst_euld(1));
Rc_ref = R_tf * R_inst;
Tc_ref = -1 * Rc_ref * inst_tl;

P_ref = K * [Rc_ref Tc_ref];
P2 = P_ref;

H_ref = [P_ref(:, 1:2) P_ref(:, 4)];
H_ref_inv = inv(H_ref);



%% Project image to world
idx_img_ref = 1;
for idx_pt = 1 : size(pts_img_cor{idx_img_ref}, 2)       

    pts_wld(:, idx_pt) = img2wld(pts_img_cor{idx_img_ref}(:, idx_pt), H_ref_inv);

end


%% data buffering

for idx_cmr = 1 : size(pts_img, 2)
    
    xi{idx_cmr} = [];
    xw{idx_cmr} = [];
    
    for idx_pt = 1 : size(pts_img{idx_cmr}, 2)
        if (pts_img{idx_cmr}(1, idx_pt) ~= 0) && (pts_img{idx_cmr}(2, idx_pt) ~= 0)
            xi{idx_cmr}(:, end+1) = pts_img{idx_cmr}(:, idx_pt);
            xw{idx_cmr}(:, end+1) = pts_wld(:, idx_pt);
        end
    end
end


%% Extrinsic parameter extraction

for idx_cmr = 1 : size(pts_img, 2)

    % Conditioning threshold for view rejection
    thresh_cond = 1e6;
    check_cond = 1;

    % N_points_views(kk) = size(x_kk,2);

    % Extrinsic parameter initialization
    [omc_init,Tc_init,Rckk_init] = compute_extrinsic_init(xi{idx_cmr}, xw{idx_cmr}, [fx; fy], [cx; cy], [k; p], alpha_c);

    % Refine initial extrinsic parameter
    MaxIter = 20;
    [omc,Tc{idx_cmr},Rc{idx_cmr},JJ] = compute_extrinsic_refine(omc_init, Tc_init, xi{idx_cmr}, xw{idx_cmr}, [fx; fy], [cx; cy], [k; p], alpha_c, MaxIter, thresh_cond);

    if check_cond
        if (cond(JJ)> thresh_cond)
            omc = NaN*ones(3,1);
            Tc{idx_cmr} = NaN*ones(3,1);
            fprintf(1,'\nWarning: View is ill-conditioned. This image is now set inactive.\n')
        end
    end

%     omc = omc_init;
%     Tc{idx_cmr} = Tc_init;
%     Rc{idx_cmr} = Rckk_init;

    %% Convert extrinsic matrix to euler angle and translation vector

    
    R_est = rodrigues(omc);

    % R_est = R_tf * R_inst
    R_inst = R_tf' * R_est;

    euld(3:-1:1,idx_cmr) = rad2deg(rotm2eul(R_inst, 'ZYX')');
%     euld(3:-1:1,idx_cmr) = rad2deg(rot2eul(R_inst, 'ZYX'));
    
    
    % rx = atand(R(3,2)/R(3,3));
    % ry = atand(-R(3,1)/sqrt(R(3,2)*R(3,2)+R(3,3)*R(3,3)));
    % rz = atand(R(2,1)/R(1,1));

    R_inst = rotz(euld(3,idx_cmr)) * roty(euld(2,idx_cmr)) * rotx(euld(1,idx_cmr));
    % R2 = eul2rotm(eul, 'ZYX');

    Rc{idx_cmr} = R_tf*R_inst;
    P{idx_cmr} = K * [Rc{idx_cmr} Tc{idx_cmr}];
    tl(:,idx_cmr) = -Rc{idx_cmr}'*Tc{idx_cmr};

    if 1
        fprintf(1,'\n[CAM%d] Extrinsic parameters\n', idx_cmr);
        fprintf(1,'Rotation angle deg (roll):  rx = [ %3.5f ]\n', -euld(1,idx_cmr));
        fprintf(1,'Rotation angle deg (pitch): ry = [ %3.5f ]\n', -euld(2,idx_cmr));
        fprintf(1,'Rotation angle deg (yaw):   rz = [ %3.5f ]\n', -euld(3,idx_cmr));
        fprintf(1,'Translation vector:          t = [ %3.5f   %3.5f   %3.5f ]\n', tl(:,idx_cmr));
    end
    
    ext_param(:,idx_cmr) = [-euld(1,idx_cmr); -euld(2,idx_cmr); -euld(3,idx_cmr); tl(1,idx_cmr); tl(2,idx_cmr); tl(3,idx_cmr) ];
end

ext_param_table = table(ext_param, 'rowNames', {'rx','ry','rz','tx','ty','tz'});




%% Image plot
if SHOW_IMAGE
    
    
    
    for idx_cmr = 1 : size(input_image_undistorted,2)

        figure('Name', ['CAM' num2str(idx_cmr)],'NumberTitle','off')
        imwrite(input_image_undistorted{idx_cmr}, [input_file_name{idx_cmr} '_undistorted.jpg']);
        image_handle = imshow(input_image_undistorted{idx_cmr});
        hold on

        %% %%%%%%%%%% Draw x-y grid on the image %%%%%%%%%%
        if SHOW_GRID

            % Downward Grid definition
            grid_x_range = 10 : 5 : 150;
            grid_y_range = -20 : 1 : 20;
            Draw_Image_Grid(grid_x_range, grid_y_range, P{idx_cmr}, image_handle)


            % Upward grid definition

            if SHOW_GRID_2
                % Grid definition
                grid_x_range = 10 : 5 : 150;
                grid_y_range = -15 : 1 : -1;
                Draw_Image_Grid(grid_x_range, grid_y_range, P{idx_cmr}, image_handle)
            end
        end
        
        
        %% %%%%%%%%%% Draw image correspondence %%%%%%%%%%
        if SHOW_IMAGE_PTS
                
            scatter(pts_img_cor{idx_cmr}(1,:), pts_img_cor{idx_cmr}(2,:), 50, 'filled', 'parent', image_handle.Parent);
                
        end

        
    end
end


%%
if SHOW_3D_VIEW
    
    fHandle = figure('Name', 'World coordinates for downward', 'NumberTitle', 'off');
    hold on;
    axis('equal');
%     title('Extrinsic parameters in world coordinates');
    view(-115,20);
    grid on;
%     axis vis3d;    
    axis([-20 100 -20 20])
    set(fHandle, 'color',[1 1 1]);
    rotate3d on;

    
    % Draw reference marker
    plot3(pts_wld(1,:),pts_wld(2,:),pts_wld(3,:),'.', 'color', colors(rem(idx_pt-1,6)+1));
    
    
    % Draw camera
    for idx_cmr = 1 : size(input_image_undistorted,2)
        
        dX = norm(Tc{idx_cmr}) / 10;
        dY = dX;

        IP = 5*dX*[1 -alpha_c 0;0 1 0;0 0 1]*[1/fx 0 0;0 1/fy 0;0 0 1]*[1 0 -cx;0 1 -cy;0 0 1]*[0 nx-1 nx-1 0 0 ; 0 0 ny-1 ny-1 0;1 1 1 1 1];
        BASE = 5*dX*([0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]);
        IP = reshape([IP;BASE(:,1)*ones(1,5);IP],3,15);

        BASE = Rc{idx_cmr}' * (BASE - Tc{idx_cmr});
        IP = Rc{idx_cmr}' * (IP - Tc{idx_cmr});
        
        if idx_cmr == 1
            cmr_pos_origin = BASE(:,1);
        end

        plot3(BASE(1,1:2),BASE(2,1:2),BASE(3,1:2),'r-','linewidth',2);
        plot3(BASE(1,3:4),BASE(2,3:4),BASE(3,3:4),'g-','linewidth',2);
        plot3(BASE(1,5:6),BASE(2,5:6),BASE(3,5:6),'b-','linewidth',2);
        plot3(IP(1,:),IP(2,:),IP(3,:),'color', [0.5 0.5 0.5],'linewidth',1);
        
        text_pos = BASE(:,1) + 7 * (BASE(:,1) - cmr_pos_origin);
        text_pos(3) = text_pos(3) + 3;
        text(text_pos(1), text_pos(2), text_pos(3), num2str(idx_cmr));
       
    end
end


function [eul] = rot2eul(R, seq)
    if strcmp(seq, 'ZYX')
        sy = sqrt(R(1,1) * R(1,1) + R(2,1) * R(2,1));

        if sy > 1e-6
            eul(1) = atan2(R(3, 2), R(3, 3));
            eul(2) = atan2(-R(3, 1), sy);
            eul(3) = atan2(R(2, 1), R(1, 1));
        else
            eul(1) = atan2(-R(2, 3), R(2, 2));
            eul(2) = atan2(-R(3, 1), sy);
            eul(3) = 0;
        end
    else
        error('unsupported rotation sequence')
    end
end


