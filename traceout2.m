function x = traceout2(W, sys, dim)
% Traces out subsystem sys of dimensions dim(sys) from matrix W, tensor Idendity of dimension dim(sys) at place where sys was.
% sys is the subsystem to be traced out: the position in dim in which its dimension is stored
% dim is an array with the dimensions of each subsystem


% x1 = syspermute(tensor(eye(dim(sys(1))), TrX(W, sys(1), dim)), [2 1 3], [dim(sys(1)) prod(dim(1:sys(1)-1)) prod(dim(sys(1)+1:length(dim)))])/dim(sys(1));
% x = syspermute(tensor(eye(dim(sys(2))), TrX(x1, sys(2), dim)), [2 1 3], [dim(sys(2)) prod(dim(1:sys(2)-1)) prod(dim(sys(2)+1:length(dim)))])/dim(sys(2));
x1 = traceout(W, sys(1), dim); 
x = traceout(x1, sys(2), dim);
