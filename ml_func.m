function f=ml_func(aa,z,n,eps0)
aa=[aa,1,1,1]; a=aa(1); b=aa(2); c=aa(3); q=aa(4);
f=0; k=0; fa=1; aa = aa(1:4); if nargin<4, eps0=eps; end
if nargin<3, n=0; end
if n==0
    while norm(fa,1)>=eps0
        fa=gamma(k*q+c)/gamma(c)/gamma(k+1)/gamma(a*k+b) *z.^k;
        f=f+fa; k=k+1;
    end
    if ~isfinite(f(1))
        if c*q==1
            f=mlf(a,b,z,round(-log10(eps0))); f=reshape(f,size(z));
        else
            error('Error: truncation method failed');
        end
    end
else
    aa(2) = b+n*a; aa(3)=c+q*n;
    f=gamma(q*n+c)/gamma(c)*ml_func(aa,z,0,eps0);
end

end
