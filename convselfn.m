function y = convselfn(u,n,dt)
%CONVSELFN Convolve signal with itself N times

% Convert to column vector
u = colv(u);

if n==0
    y = u;
elseif n>0
    for k=1:n
        if k==1
            y = dt * conv(u,u);
        else
            y = dt * conv(u,y);
        end
        y =  y(1:length(u));
    end
end

end


