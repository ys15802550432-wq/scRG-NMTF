function acc = get_acc(y_true, y_pred)
    y_true = double(y_true); y_pred = double(y_pred); n = max(max(y_true), max(y_pred)); C = zeros(n, n);
    for i = 1:length(y_true), C(y_pred(i), y_true(i)) = C(y_pred(i), y_true(i)) + 1; end
    try, [a,~]=matchpairs(-C, 1e9); acc=sum(C(sub2ind(size(C), a(:,1), a(:,2))))/length(y_true); catch, acc=sum(max(C,[],2))/length(y_true); end
end