function map_new = fill_small_holes(map, N)
filled_map = bwfill(map,'holes');
holes = filled_map &~map;
bigholes = bwareaopen(holes, N);
smallholes = holes&~bigholes;
map_new = map | smallholes;
end