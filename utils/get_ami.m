function ami = get_ami(true_labels, pred_labels)
    true_labels = double(true_labels); pred_labels = double(pred_labels); n = length(true_labels); [~,~,u1] = unique(true_labels); [~,~,u2] = unique(pred_labels);
    C = accumarray([u1, u2], 1); a = sum(C, 2); b = sum(C, 1); MI = 0; 
    for i = 1:length(a), for j = 1:length(b), if C(i,j)>0, MI=MI+(C(i,j)/n)*log((C(i,j)*n)/(a(i)*b(j))); end; end; end
    Ha = -sum((a/n).*log(a/n+eps)); Hb = -sum((b/n).*log(b/n+eps)); EMI = 0; 
    for i = 1:length(a), for j = 1:length(b)
            start_val = max(1, a(i)+b(j)-n); end_val = min(a(i), b(j));
            if start_val <= end_val, nij = (start_val:end_val);
                log_num = gammaln(a(i)+1) + gammaln(b(j)+1) + gammaln(n-a(i)+1) + gammaln(n-b(j)+1);
                log_den = gammaln(n+1) + gammaln(nij+1) + gammaln(a(i)-nij+1) + gammaln(b(j)-nij+1) + gammaln(n-a(i)-b(j)+nij+1);
                prob = exp(log_num - log_den); log_term = log( (n .* nij) ./ (a(i) * b(j)) ); EMI = EMI + sum( (nij / n) .* log_term .* prob );
            end
    end; end
    avg_H = (Ha + Hb) / 2; if avg_H - EMI < 1e-10, ami = 0; else, ami = (MI - EMI) / (avg_H - EMI); end
end