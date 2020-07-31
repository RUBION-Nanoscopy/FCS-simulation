function I = gaussint(x,y,z,lambda, NA, varargin)
    faktor=1;
    if nargin > 5 
        faktor = varargin {2}; 
    end

    FWHM = lambda/(2 * NA);
    w0   = sqrt(2 * log(2)) * FWHM;
    wd   = w0/faktor;
    r    = sqrt(x.^2 + y.^2);
    z0   = pi * w0^2/lambda;
    wz   = wd.*sqrt(1 + (z./z0) .^2);
    I    = r<=1.22/2*lambda/NA * (wd./wz).^2 .* exp(-2 .* r .^2./(wz.^2));
    %I = (wd./wz).^2 .* exp(-2.*r.^2./(wz.^2));


