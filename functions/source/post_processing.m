

function [bool_temp, map, inliers1, inliers2] =   post_processing(img_rgb, para)
img = rgb2gray(img_rgb);
iptsetpref('ImshowBorder','tight');
para.r_ratio = 16;
r_ratio = para.r_ratio;
p1 = para.p1; p2 = para.p2;
imclose_para = 5;
max_sup_radius= para.max_sup_radius;

[w, h,~] = size(img_rgb);

inliers1 = [];
inliers2 = [];
bool_temp=0;
map = zeros(w,h);
if ~(size(p1,2)<5)
    
    img_r = filter2(fspecial('average',3),double(img_rgb(:,:,1)));
    img_g = filter2(fspecial('average',3),double(img_rgb(:,:,2)));
    img_b = filter2(fspecial('average',3),double(img_rgb(:,:,3)));
    
    map_1 = zeros(w,h);
    map_2  = zeros(w,h);
    in_image_idx = p1(1,:)<h & p1(2,:)<w & p2(1,:)<h & p2(2,:)<w;
    p1 = p1(:,in_image_idx,:);
    p2 = p2(:,in_image_idx,:);
    p=[p1(1:2,:) p2(1:2,:)]';
    n_matches = size(p1,2);
    n_ranc = max(15, round(n_matches/30));
    n_ranc = min([n_ranc, max(10,round(n_matches/30)),40, size(p1,2)]);
    
    half_dis = para.min_clone_dis/2+55;
    dis_map = pdist(p);
    dis_map_q = squareform(dis_map <= half_dis);
    dis_map_q = dis_map_q | eye(size(dis_map_q,1));
    neighbors = sum(dis_map_q);
    neighbors_p1 = neighbors(1:size(dis_map_q,2)/2);
    neighbors_p2 = neighbors(size(dis_map_q,2)/2+1: end);
    
    if n_matches < 4
        return;
    end
    indx = randperm(n_matches);
    indx = indx(1:n_ranc);
    min_size = para.min_size;
    min_inliers = para.min_inliers;
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
        if size(cur_p1,2) < 4
            fprintf('find matches less 4 in a local region, continue!\n ');
            continue;
        end
        
        t  = 0.05;
        [H, inliers1] = ransacfithomography2(cur_p1(1:3,:), cur_p2(1:3,:), t);
        if isempty(H) || length(inliers1)<4%para.min_inliers
            fprintf('the number of inliers is small, continue!\n ');
            continue;
        end
        
        %%%%anayze the tranformation matrix (scale)%%%%%%%%%%%%%%%%
        allowed_err = para.allowed_err;
        [U,S,V] = svd(H(1:2,1:2));
        if para.check_scale&&(S(1,1)>8||S(1,1)<0.2 || S(2,2)>8||S(2,2)<0.2)
            %fprintf('estimated H £¨scale£© is not accurate, give up\n');
            continue;
        end
        cur_p1 = cur_p1(:,inliers1); cur_p2 = cur_p2(:,inliers1);
        
        rotation_xy = U*V';
        scale_err = max(abs(S([1 4])-1))*para.scale_err_ratio;
        r_xy = acos(rotation_xy(1))/pi*180;
        rotation_err = max(0,r_xy-5)*0.2;
        allowed_err = max(allowed_err,min(30,allowed_err + floor(scale_err+ rotation_err))); %50
        if ~(r_xy<30 || r_xy>70&&r_xy<105 || r_xy>330)
            min_inliers = max(4, para.min_inliers-1);
        end
        
        if para.check_orienation
            indx_o  = check_orientation(rotation_xy(1), rotation_xy(2), cur_p1(5,:), cur_p2(5,:), 15);
            if sum(indx_o)<length(indx_o)*0.8
                fprintf('check dominant orientation fail, continue!\n ');
                continue;
            end
            %                     if sum(indx_o)<length(indx_o) && sum(indx_o)>4
            %                         cur_p1 = cur_p1(:,indx_o); cur_p2 = cur_p2(:,indx_o);
            %                         [H, inliers1] = ransacfithomography2(cur_p1(1:3,:), cur_p2(1:3,:), t);
            %                     end
        end
        
        %%do noise estimation only for those light attacks£¬ only
        %%work for smooth images
        if S(1,1)<1.2&&S(1,1)>0.8&& S(2,2)<1.2&&S(2,2)>0.8&&(r_xy<30||r_xy>330)
            noise_err = noise_estimation(img, cur_p1(1:2,:),cur_p2(1:2,:));
            allowed_err = allowed_err + noise_err;
            if allowed_err>60
                min_inliers=4;
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
        
        if para.check_orienation && ~isempty(cur_x1)
            indx_o  = check_orientation(rotation_xy(1), rotation_xy(2), cur_x1(5,:), cur_x2(5,:), 20);
            cur_x1 = cur_x1(:, indx_o);
            cur_x2 = cur_x2(:,indx_o);
            %                     if sum(~indx_o)>0
            %                         fprintf('kick out some outliers through dominant orientation!\n ');
            %                     end
        end
        
        if size(cur_x1,2)<max(min_inliers,floor(n_matches*para.match_ratio)) && size(cur_x1,2)<para.min_total_inliers
            %fprintf('detected random match pairs, give up!\n');
            continue;
        end
        
        if size(cur_x1,2)<17&&para.check_dis && ~check_distance(cur_p1(1:2,:), cur_p2(1:2,:), para.min_dis);
            %fprintf('check smaller distance fail, continue!\n ');
            continue;
        end
        
        if size(cur_x1,2)<7&& para.check_edges&&allowed_err==para.allowed_err&&~check_vh_edges(cur_x1(1:2,:),cur_x2(1:2,:), 25)
            %fprintf('suspect vertical or horizontal edge, give up!\n');
            continue;
        end

        if S(1,1)>1.1||S(1,1)<0.9 || S(2,2)>1.1||S(2,2)<0.9
            imclose_para = 7;
            max_sup_radius = para.max_sup_radius+12;
           if S(1,1)>2.8 || S(2,2)>2.8
                imclose_para = 15;
                max_sup_radius = para.max_sup_radius+50;
                r_ratio= r_ratio*2;
            end
        end
        
        for i = 1:para.bi_direction     %bi-direction transform
            H_cur=H;
            if i == 2
                H_cur = inv(H);
                temp = cur_x1;
                cur_x1 = cur_x2;
                cur_x2 = temp;
            end
            %%%%%%the pixels difference%%%%%%%%%%%%%%%%%%%%%%
            map_ori = zeros(w,h);
            cur_radius = cur_x1(4,:)*r_ratio;
            cur_radius(cur_radius>para.max_sup_radius) = max_sup_radius;
            cur_radius(cur_radius<para.min_sup_radius) = para.min_sup_radius;
            circles = int32([cur_x1(1,:); cur_x1(2,:); cur_radius]');
            map_ori = insertShape(map_ori, 'FilledCircle',circles);
            map_ori = im2bw(map_ori);
            
            max_error = 16;%;
            sus_region_1 =  find(map_ori ==1);
            r_indx = floor((sus_region_1-1)/w)+1;
            c_indx = sus_region_1 - (r_indx-1)*w;
            copy_region = round(H_cur * [r_indx'; c_indx'; ones(1, length(r_indx))]);
            sus_region_2 = copy_region(1,:)*w+copy_region(2,:);
            
            
            sus_region_2 = sus_region_2(1<=sus_region_2 & sus_region_2<= w*h);
            sus_region_1 = sus_region_1(1<=sus_region_2 & sus_region_2<= w*h);
            
            dif_r = abs(img_r(sus_region_1) - img_r(sus_region_2'));
            dif_g = abs(img_g(sus_region_1) - img_g(sus_region_2'));
            dif_b = abs(img_b(sus_region_1) - img_b(sus_region_2'));
            
            refind_indx = dif_r< max_error& dif_g<max_error&dif_b<max_error;
            map_ori(sus_region_1(~refind_indx)) =0;
            
            map_copymove = zeros(w,h);
            map_copymove(sus_region_2(refind_indx))=1;
            if i==2
                temp = map_ori;
                map_ori = map_copymove;
                map_copymove = temp;
            end
            map_1 = map_1 | map_ori;
            map_2 = map_2 | map_copymove;
        end
    end
    if sum(sum(map_1))>0 && sum(sum(map_2))>0 && (S(1,1)<4 || S(2,2)<4)             
        map_1 =  bwareaopen(map_1, min_size);
        map_2 =  bwareaopen(map_2, min_size);
    end
    if sum(sum(map_1))>0 && sum(sum(map_2))>0
        map = map_1 | map_2;
        map= imclose(map,strel('disk', imclose_para));
        map =fill_small_holes(map, para.max_hole_size);
        bool_temp=1;
    end
    if sum(sum(map_1))>0 && sum(sum(map_2))>0 && (S(1,1)>=4 || S(2,2)>=4)             
        map =  bwareaopen(map, min_size);
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
   
   function noise_err = noise_estimation(img, p1, p2)
   [w,h] = size(img);
   
   s = 16;
   indx = p1(1,:)>s&p1(2,:)>s&p1(1,:)<w-s&p1(2,:)<h-s...
       &p2(1,:)>s&p2(2,:)>s&p2(1,:)<w-s&p2(2,:)<h-s;
   p1 = p1(:,indx); p2=p2(:,indx);
   N = size(p1,2);
   record_ratio = ones(1,N);
   for i =1:N
       patch_1 = img(p1(1,i)-s:p1(1,i)+s, p1(2,i)-s:p1(2,i)+s);
       patch_2 = img(p2(1,i)-s:p2(1,i)+s, p2(2,i)-s:p2(2,i)+s);
       [grad_tr_x1, grad_tr_y1] = gradient(double(patch_1));
       mag_tr1 = mean(sum(sum(sqrt(grad_tr_x1.^2+grad_tr_y1.^2))));
       [grad_tr_x2, grad_tr_y2] = gradient(double(patch_2));
       mag_tr2 = mean(sum(sum(sqrt(grad_tr_x2.^2+grad_tr_y2.^2))));
       record_ratio(i) = max(mag_tr1, mag_tr2)/max(1,min(mag_tr1, mag_tr2));
   end
   record_ratio = sort(record_ratio);
   noise_ratio = mean(record_ratio(2:end-1));
   noise_err = 0;
   if noise_ratio > 10000
       noise_err = 100;
   elseif noise_ratio > 5000
       noise_err = 60;
   elseif noise_ratio > 1000
       noise_err = 30;
   end
   end
   function check_pass =  check_distance(p1, p2, thr)
   check_pass = true;
   abs_dis = abs(p1-p2);
   dis = sqrt(abs_dis(1,:).^2+abs_dis(2,:).^2);
   if mean(dis) < thr
       check_pass = false;
   end
   end
   function check_pass =  check_vh_edges(p1, p2, thr)
   check_pass = true;
   N = size(p1,2);
   temp1 = min([var(p1(1,:)), var(p1(2,:)), var(p2(1,:)), var(p2(2,:))]); %an vetical or horizontal edge
   abs_p1p2=abs(p1-p2);
   temp2 = max(var(abs_p1p2(1,:)), var(abs_p1p2(2,:)));
   if temp1<thr && temp2<1
       check_pass = false;
   end
   end
