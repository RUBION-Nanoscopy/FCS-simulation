steps   = 1;
Winkel  = 0;
lambda  = .57; % Wavelength
Dt      = 1e-03;


phsize  = 1; % in Units of Airy disks

fitfunc = @(tau,N,x)(1/N * (1+x./tau).^-1);

cx      = max(scan01.xdata_lin)/2;

cy      = max(scan01.ydata_lin)/2;
F       = griddedInterpolant(scan01.xdata_grid, scan01.ydata_grid, scan01.zdata_grid);

G_all   = [];
tau     = zeros(steps, 1);
fit_tau = {};

for i = 1:1
    G   = FCS.multipletau(I_all(:, i), Dt);
    err = 0;
    [ft, gof, out] = fit(G(:,1),G(:,2),fitfunc,'Start',[1e-4, 1], 'Lower',[1e-6 1e-8]);
    if err < 3
        tau(i)  = ft.tau;
        fit_tau = ft;
        N(i)    = ft.N;
    else
        tau(i)  = NaN;
    end
    G_all(:,i)  = G(:,2)*N(i);
end