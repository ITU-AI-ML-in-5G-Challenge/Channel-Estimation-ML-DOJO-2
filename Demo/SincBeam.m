function [ p ] = SincBeam( N, w )
%SINCBEAM Summary of this function goes here
%   Detailed explanation goes here

p=sinc(w*((0:N-1)-(N-1)/2)/(2*pi))';
% p=Steer(p, -w/2); % This is for some comparisons
p=bsxfun(@rdivide, p, sqrt(sum(abs(p).^2, 1)));

end