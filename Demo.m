clear all

demo_0 = 1;
rad = pi/180;
deg = 180/pi;

pos = [30*pi/180, 120*pi/180, 200]';
vn = [10, 10, 10];


att0 = [10, 10, 10]'*rad;
v0 = [10, 20, 30]';
% ans = euler2dcm(att0)

b = fir1(20, 0.01, 'low');
b = b/sum(b);
% x = repmat([att0; vby]', length(b), 1);

% rv = [10, 2, 1]'*pi/180;
% sqrt(rv'*rv)

wm = [1, 1, 1
      2, 2, 2
      3, 3, 3
      4, 4, 4
      5, 5, 5];
%   
% wmm = wm(1:2, :)
% wm1 = wm(2:3, :)
% cross(wmm,wm1, 2)
% cs = [ [   2,    0,    0,    0,     0]/3
%        [   9,   27,    0,    0,     0]/20
%        [  54,   92,  214,    0,     0]/105
%        [ 250,  525,  650, 1375,     0]/504
%        [2315, 4558, 7296, 7834, 15797]/4620 ];
%    
%     cs(3, 1:3)*wm(1:3, :)
%     wm(1:3, :)

% eth = earth(pos, vn)
% skew(v0)


% u = [1, 2, 3]';
% u = u/norm(u);
% theta = 0.8*pi;
% 
% phi1 = theta;
% phi2 = -(2*pi - theta) ;
% 
% PHI1 = phi1*u;
% PHI2 = phi2*u;
% % norm(PHI)
% qq1 = rv2q(PHI1)'
% qq2 = rv2q(PHI2)'
% 
% PHI1'*PHI1
% rk = [[0.1; 0.1; 0.1]; [10; 10; 10]]';
% Rk = diag(rk)^2;

% R =1 + 0.1.*randn(10000, 1);
% mean(R)
% var(R)
% 
% clearvars -except R
% imu_err = imuerror();
name = 'selfdefine';
exist('name', 'var') && ~strcmp(name, '') && ~strcmp(name, 'zero')

