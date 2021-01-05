function [vnew, fnew, meanedge, stdev] = Remesher(V, F, edgelength, iterations)

%clean up patch
[vnew, fnew] = CleanPatch(V, F);
voriginal = vnew;
foriginal = fnew;
[vnew,fnew] = SubdivideLarge( vnew, fnew,edgelength,voriginal,foriginal );
voriginal = vnew;
foriginal = fnew;

for i = 1:iterations
[vnew, fnew,temp] = EdgeCollaps( vnew, fnew, edgelength ,voriginal,foriginal);
[vnew, fnew, temp] = RemoveBadTriangles( vnew, fnew,voriginal,foriginal);
[vnew,fnew] = SubdivideLarge( vnew, fnew,0 ,voriginal,foriginal);
[vnew, fnew, temp] = RemoveBadTriangles( vnew, fnew,voriginal,foriginal);

% disp(['Iteration:' num2str(i) '  Output mesh: ' num2str(size(fnew,1)) ' triangles, ' ... 
%     num2str(size(vnew,1))  ' vertices.']);
end

meanedge = temp(:,1);
stdev=temp(:,2);

