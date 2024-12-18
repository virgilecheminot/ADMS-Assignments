function [MFF, MFC, MCC] = freefree(M, ndof)

MFF = M(1:ndof, 1:ndof);
MFC = M(1:ndof, ndof+1:end);
MCC = M(ndof+1:end, ndof+1:end);