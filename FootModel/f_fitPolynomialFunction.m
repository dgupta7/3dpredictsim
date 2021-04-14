function [f_poly,err,varargout] = f_fitPolynomialFunction(x,y,n)

ny = length(y); % number of data elements in vector 
if size(x,1) ~= ny
    ny = min(size(x,1),ny);
    disp(['Input and output size did not match, using first ' num2str(ny) ' elements.']);
end
nx = size(x,2); % number of variables
if length(n) < nx
    for i=1:nx
        n(end+1) = n(end);
    end
end

A = zeros(ny,sum(n)+nx);
b = zeros(ny,1);
ct = 0;

for i=1:nx
    for j=1:ny
        A(j,ct+1) = 1;
        b(j,1) = y(j);
        for k=1:n(i)
            A(j,ct+1+k) = x(j,i)^k;
        end
    end
    ct = ct + 1 + n(i);
end

cff = A\b;
err = sqrt(nanmean(((A*cff-b)).^2));

coeff = zeros(max(n)+1,nx);
ct = 1;
for i=1:nx
    for k=1:n(i)+1
       coeff(k,i) = cff(ct);
       ct = ct + 1;
    end
end

tm=0;

f_poly = @(q) makeFunc(coeff,q);

if nargout >2
    varargout{1} = coeff;
end

    function f = makeFunc(cf,var)
        f = 0;
        for ii=1:size(cf,1)
            for jj=1:size(cf,2)
                f = f + cf(ii,jj)*var(:,jj).^(ii-1);
            end
        end
    end
end