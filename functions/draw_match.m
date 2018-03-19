function [] = draw_match(img, p1, p2)
figure;
imshow(img);
hold on
[w,h,~] = size(img);
line_width = max(1, int8(w/1024 * h/768));
line_width = min(line_width, 1);
for i = 1: size(p1,2)
    line([p1(1,i)' p2(1,i)'], [p1(2,i)' p2(2,i)'], 'color', 'blue','LineWidth', line_width);
end
scatter(p1(1,:),p1(2,:),'r');
scatter(p2(1,:),p2(2,:),'r');

     