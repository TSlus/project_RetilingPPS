function [projections] = Project(vS,fS,vT,fT)

% global vglobalRP;
% vT = vglobalRP; % size(vglobalRP)
% size(vT)
% global vba;
% vT = vba;
% size(vT)

TRS = triangulation(fS,vS); 
normalsS = vertexNormal(TRS);

[IDXsource,Dsource] = knnsearch(vT,vS);
vector_s_to_t = vT(IDXsource,:) - vS;

projections = vS + [(sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,1) (sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,2) (sum(vector_s_to_t.*normalsS,2)./(norm(normalsS).^2)).*normalsS(:,3)];

end

