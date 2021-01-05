% force_test
aplane = sum(forces(1:nCand, :).*norm_face(nameF_cand3(1:nCand), :),2);
[xdelta,y] = max(abs(aplane));


