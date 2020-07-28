function V = getV(Q)
[~, ~, V] = svd(Q');
Qhat = inv(V)*Q;
if Qhat(1,1) > Qhat(1,end) % flip the eiganvectors if the first one is in the negative direction
    V = -V;
end
end