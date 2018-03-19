function [clusters, pixel_cls] =  gray_cluster(im, idx, step, len ,min_pixels)
pixel_vals = im(round(idx));
cls = [0:step:255]';
cls = cls-len;
cls = [cls, cls+step+len];
cls(1) =0; cls(end)=255;
num_cls = size(cls,1);
bin_s  = zeros(1, num_cls);
clusters =  {};
pixel_cls = cls;
re_cls =[];
for i = 1:num_cls
    cur_range          =    cls(i,:);
    cur_samples     =    find(pixel_vals >= cur_range(1) & pixel_vals <= cur_range(2));
    if length(cur_samples) <min_pixels
        re_cls =[re_cls, i];
        continue;
    end
    clusters = [clusters, cur_samples];
end
pixel_cls(re_cls,:)=[];
end


