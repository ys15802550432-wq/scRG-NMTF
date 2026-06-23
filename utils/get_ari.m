function s = get_ari(t, p)
    n=length(t); [~,~,u1]=unique(t); [~,~,u2]=unique(p); m=accumarray([u1, u2], 1);
    a=sum(m,2); b=sum(m,1); nc2=@(x) x*(x-1)/2; t1=nc2(n); t2=sum(arrayfun(nc2, a)); t3=sum(arrayfun(nc2, b)); nn=sum(arrayfun(nc2, m(:)));
    if t1==(t2+t3)/2, s=0; else, s=(nn-(t2*t3)/t1)/((t2+t3)/2-(t2*t3)/t1); end
end