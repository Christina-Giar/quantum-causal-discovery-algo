function x = traceout(W, sys, dim)
% Traces out subsystem sys of dimensions dim(sys) from matrix W, tensor Idendity of dimension dim(sys) at place where sys was.
% sys is the subsystem to be traced out: the position in dim in which its dimension is stored
% dim is an array with the dimensions of each subsystem


x = syspermute(tensor(eye(dim(sys)), TrX(W, sys, dim)), [2 1 3], [dim(sys) prod(dim(1:sys-1)) prod(dim(sys+1:length(dim)))])/dim(sys);