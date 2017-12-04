 function x = traceoutcomN(W, Sys, Dim)
% W is the matrix, sys the system whose complementary will be traced out, dim is the dimension [d_in d_out d_in d_out etc.]
% output is not normalized.

Sys = sort(Sys);
s = 1;
i = 1;
product = 1;
sysT = zeros(1, length(Dim) - length(Sys));
while s <= length(Dim);
    for o = 1:length(Sys);
        if s == Sys(o);
            s = s+1;
        end
    end
    if s <= length(Dim);
        sysT(i) = s;
        i = i+1;
        s = s+1;
        
    end
end
sysT = sort(sysT);
x = TrX(W, sysT, Dim);% TrX just traces out, leaving W unnormalized.