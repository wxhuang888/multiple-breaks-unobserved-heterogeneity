function [theta,se,sig2]=b_est(y,y1)
NT=size(y,1);
theta=y1'*y1\(y1'*y);
e=y-y1*theta;
% sig2 = sum(e.^2)/NT;
sig2 = sum(e.^2);
% se=sqrt(diag(y1'*y1\(e'*e)/(NT-4)));
se = inv(y1'*y1) * (e'*e)/(NT);