function G = multipletau(data, deltat, varargin)

% Autocorrelation using multiple tau algorithm.
%
% This function calculates the autocorrelation of 
% the vector data on a logarithmic scale using the 
% multiple tau algorithm. 
%
% The algorithm is copied from multipletau.py
% (http://paulmueller.github.io/multipletau/) from Paul Müller.
%
% The original license is:
% Copyright (c) 2014 Paul Müller
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
%  met:
%
%  1. Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   
%  2. Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in
%     the documentation and/or other materials provided with the
%     distribution.
%
%  3. Neither the name of multipletau nor the names of its contributors 
%     may be used to endorse or promote products derived from this 
%     software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL INFRAE OR
%CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
%EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
%PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
%PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
%LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
%NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The license of the matlab version is the same, but Copyright by Patrick
% Happel, 2015
%
% Synopsis:
%   [lags, correlation] = multipletau(data, deltat)
%     Computes the autocorrelation of the vector data on a logarithmic
%     scale using the multiple tau algorithm, each bin width is deltat. 
%     Returns a vector containing the lags and the corresponding
%     correlation. 
%
% [lags, correlation] = multipletau(data, options)
%     As above, but with additional options. Options are a strings with a
%     key followed by a value.
%     
%     Currently, the options are not yet implemented. They will have the
%     following meanings:
%      
%     Key         Value                                      Default
%    ----------------------------------------------------------------------
%     normalize   0 or 1, indicating whether the data        1
%                 should be normailzed                      
%     
%     m           an even number indicating the points       16
%                 on one level    


options.normalize=1;
options.m = 16;
options.deltat = deltat;

% Calulate the maen of the data

avg = mean(data);
m = options.m;
if avg == 0
    error('Mean of the data is zero')
end

fdata = double(data);

if mod(m,2) ~= 0
   oldbin = m;
   options.binsize = 2*(round(m/2)+1);
   warning('Changing binsize from %f to %i', oldbin, m)
end

% Arragne data
if size(fdata,2) > size(fdata,1)
    fdata = fdata';
end

N0 = size(fdata,1);
N = N0;

k = floor(log2(N/m));

lenG = floor(m + k*m/2);

G = zeros(lenG, 2);

normstat = zeros(lenG, 1);
normnump = zeros(lenG, 1);

if options.normalize == 1
    fdata = fdata - avg;
end

if N < 2*m
    error('Size of data must be larger than 2* binsize')
end

for n = 1:m
    G(n, 1) = options.deltat * (n);
    G(n, 2) = sum(fdata(1:end-n).*fdata(n+1:end));
    normstat(n) = N - n;
    normnump(n) = N;
end

if mod(N,2) == 1
    N = N - 1;
end

fdata = (fdata(1:2:N)+fdata(2:2:N+1))/2;

N = N/2;
idx = m;

for step = 1:k
    %step%length(fdata)
    for n = 1:floor(m/2)
        L = floor(n + m/2);
        
        if isempty(fdata(1:N-L))
            
            %size(G)
            G = G(1:idx-1,:);
            
            break;
        else
            
            idx = idx + 1;
            G(idx,1) = options.deltat * (n+m/2)*2^step;
           % length(fdata(1:N-L))
           % length(fdata(L+1:end))

            G(idx,2) = sum(fdata(1:end-L).*fdata(L+1:end));
            
            normstat(idx) = N-(n+m/2);
            normnump(idx) = N;
        end
    
        if mod(N,2) == 1
            N = N - 1;
        end
    end
    fdata = (fdata(1:2:N)+fdata(2:2:N+1))/2;
    N = N/2;
end
%normstat(1:length(G))
if options.normalize==1
    
    G(:,2) = G(:,2) ./ ((avg^2)*normstat(1:length(G)));
else
    G(:,2) = G(:,2) .* (N0/normnump(1:length(G)));
end
    

