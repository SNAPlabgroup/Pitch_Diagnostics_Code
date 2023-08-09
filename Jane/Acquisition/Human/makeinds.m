function inds = makeinds()
% Make inds matrix for words such that Jane is the target name and such
% that T, M1, and M2 are always different from each other for each category

doneflag = false;

while ~doneflag
    inds = ceil(rand(3, 5)*8);
    inds(1, 1) = 2;
    goodcols = 0;
    for k = 1:5
        if numel(unique(inds(:, k))) == 3
            if ~any( (inds(:) > 8) & (inds(:) < 1))
                goodcols = goodcols + 1;
            end
        end
    end
    if goodcols == 5
        doneflag = true;
    end
end
