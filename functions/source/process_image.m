
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% This code was developed by Li Yuanman 
% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the paper if you use this code:
% Y. M. Li and J. T. Zhou, ¡°Image Copy-Move Forgery Detection Using Hierarchical Feature Point Matching¡±, ASC, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bool_temp, map, inliers1, inliers2] = process_image(imagefile, plotimg, large_resize)
if nargin <= 2
    large_resize = 0;
end
img_rgb     =    imread(imagefile);

%%%%%%%%%%%%%%%set Parameter%%%%%%%%%%%%%%%%
para          =     setPara(img_rgb);
para.plotimg = plotimg;
if large_resize == 1       %for resize factor bigger than 2
   para.scale_seg = -1;
end
%%%%%%%%%%%%%%%do matching%%%%%%%%%%%%%%%%%%%
if para.match_strategy ==1
    [para.p1, para.p2]     =   feature_matching(imagefile, para); 
else
    [para.p1, para.p2]     =   feature_matching2(imagefile, para);
end

%%%%%%%%%%%%%draw figure%%%%%%%%%%%%%%%%%%%%%
if (plotimg==1&&size(para.p1,2)>=1)
   draw_match(img_rgb, para.p1, para.p2)
end

if para.post_strategy==1
    [bool_temp, map, inliers1, inliers2]   =   post_processing(img_rgb, para);
else 
    [bool_temp, map, inliers1, inliers2]   =   post_processing2(img_rgb, para);
end

%%%%%%%%%%%%postprocessing%%%%%%%%%%%%%%%%%%%%



function para = setPara(img)
    [w, h,~] = size(img);
    
    para.max_sup_radius = 20;
    para.min_sup_radius = 16;
    para.max_hole_size = 300;
    para.max_keypoints_octave = 7000;
    para.r_ratio = 16;            
    para.check_orienation  = 1;
    para.match_strategy = 1; 
    para.post_strategy = 1; 
    para.allowed_err = 1;
    para.check_scale = 1; 
    para.check_dis = 0;
    para.check_edges = 0;
    para.check_transition= 1;%for coverage database
    para.illumination=0;%for coverage database
    para.max_scale = 2;%for coverage database
    para.fillholes =0;%for coverage database
    para.allow_error_inliers=1.5;%for coverage database
    if w<840&& h<840
       para.check_orienation=1;
       para.match_strategy = 2; 
       para.post_strategy=2;
       para.bi_direction = 1;
       para.sift_method =2;
       para.match_thr = 0.74;
       para.min_inliers = 4;
       para.min_clone_dis = max(40,25*(sqrt(w*h/(1000*1000))));
       para.scale_seg = [0];
       para.min_size = 250; 
       para.min_total_inliers = 4;
       para.min_neighbors = min(3, round(3*(w*h/(1000*1000))));
       if w<600 && h<600        %for coverage database
           para.sift_method =1;
           para.check_transition=0;
           para.match_thr = 0.65;
           para.matching_rgb_dif = 70;  
           para.min_size = 2000;
           para.min_inliers = 5;
           para.min_total_inliers =6;
           para.allow_error_inliers=1.5;
           para.illumination=15;
           para.fillholes =1;
           para.max_scale = 2.5;
       end
    else 
        para.match_strategy = 2; 
        para.bi_direction = 1;
        para.match_thr = 0.52;
        para.sift_method =1;
        para.min_inliers = 5;
        para.min_clone_dis = 25*(sqrt(w*h/(1000*1000)));  
        para.scale_seg = [-1,0,1,2,3];
        para.min_size = min(500,400 + round(30*w*h/(1000*1000)));
        para.min_neighbors = min(3, round(3*sqrt(w*h/(1000*1000))));
        para.min_total_inliers = 20;
        para.scale_err_ratio = 4; 
        para.match_ratio = 1/5; %the correct match at least 1/5 of total matches
        para.half_dis_add = 50;
        para.min_dis = 25;
        para.check_dis = true;
        para.check_edges = true;
    end
    if w>=1400&& h>=1400
        para.match_strategy = 1; 
        para.half_dis_add = 80;
        para.allowed_err = 8;
        para.bi_direction = 2;
        para.min_inliers = 6;  %5
        para.max_sup_radius = 48;
        para.min_sup_radius = 20;
        para.max_hole_size = 1500;
        para.max_keypoints_octave = 16000;
        para.min_total_inliers = 30;
        para.scale_err_ratio = 30; 
        para.check_scale = 0; 
        para.match_ratio =1/20;
        para.scale_seg = [-1,0,1];
        para.match_thr = 0.6;                  %for F660, set it as 0.5, can achieve much lower FPR
        para.min_neighbors = min(2, round(3*sqrt(w*h/(1000*1000))));
        para.min_dis = 130; 
    end
    if para.bi_direction ==2
        para.min_size = round(para.min_size*1.3);
    end
    para.max_matches   = 50000;
    para.max_err_dis_xy = 100;
   


    
    
