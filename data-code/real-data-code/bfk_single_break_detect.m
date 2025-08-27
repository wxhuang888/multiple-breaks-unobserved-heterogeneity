function [est_break_pt] = bfk_single_break_detect(N,T,xmat,ymat,zmat)
z = zmat;
K = size(xmat,3);
p = size(zmat,3);
for k = 1:T-1
    for i = 1:N
        z2          = [zeros(k,p); reshape(z((k+1):T,i,:),T-k,p)];        
        Xcomb       = [reshape(xmat(:,i,:),T,K) z2];        
        bhat        = (Xcomb'*Xcomb)\(Xcomb'*ymat(:,i));
        yhat        = Xcomb*bhat;
        SSR2(i,k)   = (ymat(:,i) - yhat)'*(ymat(:,i) - yhat);        
    end
end
sumSSR              = sum(SSR2); 
[~,est_break_pt]    = min(sumSSR);