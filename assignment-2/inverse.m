% calculations to support in finding the proof for exercise 3.2

A7 = [
1 1 1 0 0 0 0;
0 1 1 1 0 0 0;
0 0 1 1 1 0 0;
0 0 0 1 1 1 0;
0 0 0 0 1 1 1;
1 0 0 0 0 1 1;
1 1 0 0 0 0 1 ]'

A71 = rref([A7 eye(7)])(:,8:14)

A9 = [
1 1 1 1 0 0 0 0 0;
0 1 1 1 1 0 0 0 0;
0 0 1 1 1 1 0 0 0;
0 0 0 1 1 1 1 0 0;
0 0 0 0 1 1 1 1 0;
0 0 0 0 0 1 1 1 1;
1 0 0 0 0 0 1 1 1; 
1 1 0 0 0 0 0 1 1;
1 1 1 0 0 0 0 0 1]'

A91 = rref([A9 eye(9)])(:,10:18)