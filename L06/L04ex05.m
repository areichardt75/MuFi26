% Lecture 04. Example 5.
% Coordinate Metrology
% Sensing machines are used to record the coordinates of points on the
% perimeter of the manufactured products. Determine how close the points
% are to being circular, we can fit a least squares circle to the data and
% check how close the measured points are to the circle. 

% Measured Data
xi = transpose([3 2.6 0.6 -0.5 -1.3 -1.1 0.3 1.6 2.3 2.8]);
yi = transpose([2.1 3.4 4.1 3.8 2.3 1.3 0.2 0.3 0.3 1.5]);

% equation to fit (x-c1)^2+(y-c2)^2 = r^2
% ...
% 
A = [2*xi 2*yi ones(size(xi))];
B = [xi.^2+yi.^2];
Cpars = (transpose(A)*A)\(transpose(A)*B);

% center of fitted circle : (c1, c2) 
% radius : r = sqrt(c3+c1^2+c2^2)
c1 = Cpars(1);
c2 = Cpars(2);
c3 = Cpars(3);
r = sqrt(c3+c1.^2+c2.^2)

% Plotting the results
figure;
  plot(xi,yi,'bo','MarkerFaceColor','b');
  hold on;
  % plot()