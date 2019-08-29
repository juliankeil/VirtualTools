%% Build SpatDimNeighStruct for findclusters

dummymat = zeros(size(neigh,2),size(neigh,2));

for c = 1:size(neigh,2)
    dummymat(c,c) = 1;
    for n = 1:size(neigh(c).neighblabel,1)
        for ic = 1:size(neigh,2)
            if strcmp(neigh(ic).label,neigh(c).neighblabel(n));
                dummymat(c,ic) = 1;
            end
        end
    end
end