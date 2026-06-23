function s = get_nmi(A, B)
    t=length(A); idA=unique(A); idB=unique(B); MI=0;
    for i=1:length(idA), for j=1:length(idB), c=sum(A==idA(i)&B==idB(j)); if c>0, p=c/t; MI=MI+p*log2(p/((sum(A==idA(i))/t)*(sum(B==idB(j))/t))); end; end; end
    Ha=0; for i=1:length(idA), p=sum(A==idA(i))/t; if p>0, Ha=Ha-p*log2(p); end; end
    Hb=0; for i=1:length(idB), p=sum(B==idB(i))/t; if p>0, Hb=Hb-p*log2(p); end; end
    s=2*MI/(Ha+Hb+eps);
end