% oddTriangle
clear; clc
load('bronze.mat')

% oddpoint = [30.3201010000000,-26.3411010000000,-831.630981000000];
oddpoint = [30.32,-26.34,-831.63];
odd_point2 = [-11.0834,-44.414501,-827.2969];
odd_point3 = [-7.5106,-38.3862,-837.9340];
[va, id] = min(sum((vertices - oddpoint).^2,2));
figure(5)
trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3))
axis equal; hold on
v = vertices(id, :);
plot3(v(:,1),v(:,2),v(:,3),'r*');
plot3(odd_point2(:,1),odd_point2(:,2),odd_point2(:,3),'k*');
plot3(odd_point3(:,1),odd_point3(:,2),odd_point3(:,3),'r*');

AZ = 51.7757;
EL = -28.4537;
view(AZ,EL);
