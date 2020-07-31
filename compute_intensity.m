steps  = 1;
Winkel = 0;
lambda = .57; 
phsize = 1; % in Units of Airy disks

fitfunc = @(tau,N,x)(1/N * (1 + x ./ tau) .^-1);
cx      = max(scan01.xdata_lin)/2;
cy      = max(scan01.ydata_lin)/2;

F = griddedInterpolant(scan01.xdata_grid, scan01.ydata_grid, ...
    scan01.zdata_grid);

G_all = [];
tau   = zeros(steps, 1);
f_tau = {};

for i = 1:1
    PHole = Pinhole(1.4, 1.518, lambda+.05, phsize, faktor);
    offx  = cx;
    offy  = cy;
    offz  = (scan01.zdata_lin(((scan01.ypx - 1) / 2) ...
            * scan01.xpx + (((scan01.xpx - 1) / 2) + 1))); 
    K     = size(c,1);
    
    for j = 1:K
        x = c{j}(:,1) - offx;
        y = c{j}(:,2) - offy;
        z = c{j}(:,3) - offz;
            
        r = sqrt(x.^2 + y.^2);
        I_all(:, i) = I_all(:, i) ...
                      +  PHole.get(r, z) ...
                     .* gaussint(x, y, z, lambda, 1.4, 1.518, faktor);
        
    end
end