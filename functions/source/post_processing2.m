
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% This code was developed by Li Yuanman 
% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please cite the paper if you use this code:
% Y. M. Li and J. T. Zhou, ¡°Image Copy-Move Forgery Detection Using Hierarchical Feature Point Matching¡±, ASC, 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bool_temp, map, inliers1, inliers2] =   post_processing2(img_rgb, para)
para.r_ratio = 16;
r_ratio = para.r_ratio;
p1 = para.p1; p2 = para.p2;

[w, h,~] = size(img_rgb);

inliers1 = [];
inliers2 = [];
bool_temp=0;
map = zeros(w,h);
if ~(size(p1,2)<4)
    
    img_r = double(img_rgb(:,:,1));
    img_g = double(img_rgb(:,:,2));
    img_b = double(img_rgb(:,:,3));
    
    map_1 = zeros(w,h);
    map_2  = zeros(w,h);
    in_image_idx = p1(1,:)<h & p1(2,:)<w & p2(1,:)<h & p2(2,:)<w;
    p1 = p1(:,in_image_idx,:);
    p2 = p2(:,in_image_idx,:);

    p=[p1(1:2,:) p2(1:2,:)]';
    n_matches = size(p1,2);
    n_ranc = max(15, round(n_matches/20));
    n_ranc = min([n_ranc, n_matches,40]);
    
    half_dis = 25*(sqrt(w*h/(1000*1000)))/2+50;
    dis_map = pdist(p);
    dis_map_q = squareform(dis_map <= half_dis);
    dis_map_q = dis_map_q | eye(size(dis_map_q,1));
%     neighbors = sum(dis_map_q);
%     neighbors_p1 = neighbors(1:size(dis_map_q,2)/2);
%     neighbors_p2 = neighbors(size(dis_map_q,2)/2+1: end);
    
    if n_matches < 4
        return;
    end
    indx = randperm(n_matches);
    indx = indx(1:n_ranc);
    min_size = para.min_size;
    for k = 1 : n_ranc
        cur_p = indx(k);
        neighbors_cur_p = dis_map_q(:,cur_p);
        cur_indx1 = neighbors_cur_p(1:n_matches);
        cur_indx2 = neighbors_cur_p(n_matches+1:end);
        cur_p1 = p1(:, cur_indx1);
        cur_p2 = p2(:, cur_indx1);
        cur_indx_changed = xor( (cur_indx1 | cur_indx2), cur_indx1);
        if sum(cur_indx_changed)>0
            cur_p1_add = p2(:, cur_indx_changed);
            cur_p2_add = p1(:, cur_indx_changed);
            cur_p1 = [cur_p1, cur_p1_add];
            cur_p2 = [cur_p2, cur_p2_add];
        end
        
        %%%%%do ransac%%%%%%%%%%%%%%%%%
        if size(cur_p1,2) < 4                              %4
            %fprintf('find matches less 4 in a local region, continue!\n ');
            continue;
        end  
        
        t  = 0.05;
        [H, inliers1] = ransacfithomography2(cur_p1(1:3,:), cur_p2(1:3,:), t);
        
        if isempty(H) || length(inliers1)<para.min_inliers
            %fprintf('the number of inliers is small, continue!\n ');
            continue;
        end
        
        %%%%anayze the tranformation matrix (scale)%%%%%%%%%%%%%%%%
        allowed_err = 0;
        [U,S,V] = svd(H(1:2,1:2));
        if S(1,1)>para.max_scale||S(1,1)<0.3 || S(2,2)>para.max_scale||S(2,2)<0.3
            %fprintf('estimated H £¨scale£© is not accurate, give up\n');
            continue;
        end
        
        rotation_xy = U*V';
        scale_err = max(abs(S([1 4])-1))*5;
        r_xy = acos(rotation_xy(1))/pi*180;
        rotation_err = r_xy*0.2;
        allowed_err = max(1.5,allowed_err + round(scale_err+ rotation_err)); %50
        if para.check_orienation 
            cur_p1 = cur_p1(:,inliers1); cur_p2 = cur_p2(:,inliers1);
            indx_o  = check_orientation(rotation_xy(1), rotation_xy(2), cur_p1(5,:), cur_p2(5,:), 15);
            if sum(indx_o)<length(indx_o)*0.8
                fprintf('check dominant orientation fail, continue!\n ');
                 continue;
            end
        end    
        
        %%%%%%get all the inliers%%%%%%%%%%%%%%%%%%%%%%
        [inliers] =  get_inliers(H, p1(1:3,:), p2(1:3,:),allowed_err);
        outliers = find(inliers ==0);
        cur_x1 = p1(:, inliers);
        cur_x2 = p2(:, inliers);
        if sum(outliers) > 0
            p1_1 = p2(:, outliers);
            p2_1 = p1(:, outliers);
            [inliers_1] =  get_inliers(H, p1_1(1:3,:), p2_1(1:3,:), allowed_err);
            if sum(inliers_1)>0
                cur_x1 = [cur_x1, p1_1(:, inliers_1)];
                cur_x2 = [cur_x2, p2_1(:,inliers_1)];
            end
        end
              
        if size(cur_x1,2)<max(para.min_total_inliers,round(n_matches/5)) && size(cur_x1,2)<20
            %fprintf('detected random match pairs, give up!\n');
            continue;
        end
        
        %%%%anayze the tranformation matrix (transition)%%%%%%%%%%%%%%%%
        tran_xy = abs(H(1:2,3));
        if para.check_transition && (tran_xy(1)>=h || tran_xy(2)>=w) %|| max(xy_err)>para.max_err_dis_xy
            %fprintf('estimated H £¨transition£© is not accurate, give up\n');
            continue;
        end
        
        
        %%%%%%the pixels difference%%%%%%%%%%%%%%%%%%%%%%
        map_ori = zeros(w,h);
        cur_radius = cur_x1(4,:)*r_ratio;
        cur_radius(cur_radius>25) = 25;
        circles = int32([cur_x1(1,:); cur_x1(2,:); cur_radius]');
        map_ori = insertShape(map_ori, 'FilledCircle',circles);
        map_ori = im2bw(map_ori);
        
        max_error = 10+para.illumination;
        sus_region_1 =  find(map_ori ==1);
        r_indx = floor((sus_region_1-1)/w)+1;
        c_indx = sus_region_1 - (r_indx-1)*w;
        copy_region = round(H * [r_indx'; c_indx'; ones(1, length(r_indx))]);
        sus_region_2 = copy_region(1,:)*w+copy_region(2,:);
        
        
        sus_region_2 = sus_region_2(1<=sus_region_2 & sus_region_2<= w*h);
        sus_region_1 = sus_region_1(1<=sus_region_2 & sus_region_2<= w*h);
        
        dif_r = abs(img_r(sus_region_1) - img_r(sus_region_2'));
        dif_g = abs(img_g(sus_region_1) - img_g(sus_region_2'));
        dif_b = abs(img_b(sus_region_1) - img_b(sus_region_2'));
        
        refind_indx = dif_r< max_error& dif_g<max_error&dif_b<max_error;
        map_ori(sus_region_1(~refind_indx)) =0;
        
        %map_ori(sus_region_2 <= w*h) =0;
        
        map_1 = map_1 | map_ori;
        
        map_copymove = zeros(w,h);
        map_copymove(sus_region_2(refind_indx))=1;
        map_2 = map_2 | map_copymove;
    end
    map_1= imclose(map_1,strel('disk', 2));
    map_2= imclose(map_2,strel('disk', 2));
    map_1 =  bwareaopen(map_1, min_size);
    map_2 =  bwareaopen(map_2, min_size);
    if sum(sum(map_1))>0 && sum(sum(map_2))>0
        map = map_1 | map_2;
        map= imclose(map,strel('disk', 5));
        if para.fillholes
            map = fill_small_holes(map, min_size);
        end
        bool_temp=1;
        iptsetpref('ImshowBorder','loose');
    end
    
end

if(bool_temp)
    fprintf('Tampering detected!\n\n');
else
    fprintf('Image not tampered.\n\n');
end
end

function[inliers] =  get_inliers(H, p1, p2, t)
Hp1 = H*p1;
invHp2 = H\p2;

p1     = hnormalise(p1);
p2     = hnormalise(p2);
Hp1    = hnormalise(Hp1);
invHp2 = hnormalise(invHp2);

d2 = sum((p1-invHp2).^2)  + sum((p2-Hp1).^2);
inliers = abs(d2) < t;
end
function [inliers] = check_orientation(cos_t, sin_t, o1, o2, t)
theta = acos(cos_t);
if cos_t<=0 && sin_t<=0 || cos_t>0 && sin_t<0
    theta = 2*pi - theta;
end
est_o = theta/pi*180;

p1_o =o1./pi.*180;
p2_o = o2./pi.*180;
p1_o(p1_o<0) = p1_o(p1_o<0) + 360;
p2_o(p2_o<0) = p2_o(p2_o<0) + 360;

dif_o = p2_o-p1_o;
dif_o(dif_o<0) = dif_o(dif_o<0) + 360;
dif_o2 = 360-dif_o;
inliers1= abs(dif_o - est_o)<t;
inliers2= abs(dif_o2 - est_o)<t;
inliers = inliers1|inliers2;
end