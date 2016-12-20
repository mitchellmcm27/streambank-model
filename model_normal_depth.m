function H_norm = model_normal_depth(Q,H_bkf,B,S,n)

for j=1:numel(Q)
    yn(1) = H_bkf;
    in = 1;
    dyn(1) = 1e-2;
    while (abs(dyn(in))>1e-4)
        An(in) = B*yn(in);
        Tn(in) = B;
        Pn(in) = B + 2*yn(in);
        Rn(in) = An(in)/Pn(in);
        Dn(in) = An(in)/(B);
        fn(in) = sqrt(S)*An(in)*Rn(in)^(2/3)*n^(-1)-Q(j);
        ffn(in) = (sqrt(S)*n^(-1))*...
            ((Rn(in)^(2/3)*Tn(in)) + ...
            (Tn(in)/Pn(in))-((2*yn(in)*Rn(in))/Pn(in)));
        yn(in+1)=yn(in)-fn(in)/ffn(in);
        dyn(in+1)=-fn(in)/ffn(in);
        in=in+1;
    end
    H_norm(j) = yn(in);
end
end