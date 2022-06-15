%% NV orientations
normalization = 1/(sqrt(3));

o1 = normalization * [1 1 1];
o2 = normalization * [-1 -1 1];
o3 = normalization * [-1 1 -1];
o4 = normalization * [1 -1 -1];
orientations = [o1; o2 ;o3; o4];