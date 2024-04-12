function [sol,esol] = evaluateNNsolution(w,xm,ym)
[fun,dfun,d2fun,d3fun,d4fun] = ActivationFun();
[v,W,u] = param(w);
[A,~,h,~,~,~,~,~,exact_sol] = setup();

B = h([xm;ym]);
NNfun = zeros(size(xm));
[n1,n2] = size(xm);
for i = 1 : n1
    for j = 1 : n2
        x = [xm(i,j);ym(i,j)];
        NNfun(i,j) = v'*fun(W*x + u) - v'*fun(W*[x(1);1] + u)...
                     - v'*(W(:,2).*dfun(W*[x(1);1] + u));
    end
end
sol = A(xm,ym) + B.*NNfun;
esol = exact_sol(xm,ym);
end
