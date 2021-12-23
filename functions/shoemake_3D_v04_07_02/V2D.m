function [X, Y, Z, Ux, Uy, Uz] = V2D(tmp)
X  = reshape(tmp(:,1),length(unique(tmp(:,1))),length(unique(tmp(:,2))),length(unique(tmp(:,3))));
Y  = reshape(tmp(:,2),length(unique(tmp(:,1))),length(unique(tmp(:,2))),length(unique(tmp(:,3))));
Z  = reshape(tmp(:,3),length(unique(tmp(:,1))),length(unique(tmp(:,2))),length(unique(tmp(:,3))));
Ux = reshape(tmp(:,4),length(unique(tmp(:,1))),length(unique(tmp(:,2))),length(unique(tmp(:,3))));
Uy = reshape(tmp(:,5),length(unique(tmp(:,1))),length(unique(tmp(:,2))),length(unique(tmp(:,3))));
Uz = reshape(tmp(:,6),length(unique(tmp(:,1))),length(unique(tmp(:,2))),length(unique(tmp(:,3))));