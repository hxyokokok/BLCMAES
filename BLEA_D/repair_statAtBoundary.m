function X = repair_stayAtBoundary(X,xrange)
    for i = 1 : size(X,1)
        flag = X(i,:) < xrange(1,:);
        X(i,flag) = xrange(1,flag);
        flag = X(i,:) > xrange(2,:);
        X(i,flag) = xrange(2,flag);
    end
end