if fig1
    figure(nfig); nfig = nfig + 1;
    trimesh(faces, vertices(:,1),vertices(:,2),vertices(:,3));axis equal
    title('初始网格')
    fig1 = 0;
end

if fig2
    figure(nfig); nfig = nfig + 1;
    trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
    hold on
    plot3(vertices_cand(:,1), vertices_cand(:,2), vertices_cand(:,3), 'b*');
    plot3(vertices(CSP_idx,1), vertices(CSP_idx,2), vertices(CSP_idx,3), 'r*');axis equal
    hold off; title('Candidate Points and CSP')
    fig2 = 0;
end

if fig3
    figure(nfig); nfig = nfig + 1;
    trimesh(faces, vertices(:,1), vertices(:,2), vertices(:,3));axis equal;
    hold on
    plot3(vertices_cand(:,1), vertices_cand(:,2), vertices_cand(:,3), 'b*');
    plot3(vertices(CSP_idx,1), vertices(CSP_idx,2), vertices(CSP_idx,3), 'ro');axis equal
    hold off
    title('Candidate Points and CSP after moving')
    fig3 = 0;
end

if fig4
    figure(nfig); nfig = nfig + 1;
    trimesh(faces_Mutual, vertices_Mutual(:,1), vertices_Mutual(:,2), vertices_Mutual(:,3));
    axis equal; title('Mutual Tesselation');
    fig4 = 0;
end

if fig5
    figure(nfig); nfig = nfig + 1;
    trimesh(faces_final, vertices_final(:,1), vertices_final(:,2), vertices_final(:,3));axis equal;
    title('mesh after kicking old points');
    fig5 = 0;
end