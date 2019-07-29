classdef Pinhole < handle
    properties (SetAccess = protected)
        data = [];
        NA = 1;
        n = 1;
        lambda = 1;
        f_airy = 1;
        max_r = 2.5^2;
        max_z = 2.5^2;
        f_res = 1;
        xy_data = [];
        rz_data = [];
    end
    methods
        function self = Pinhole(NA, n, wavelength, varargin)
            self.NA = NA;
            self.n = n;
            self.lambda = wavelength;
            self.f_airy = 1;
            self.f_res = 1;
            if nargin >= 4
                self.f_airy = varargin{1};
            end
            if nargin == 5
                self.f_res = varargin{2};
            end
        end
        
        function generateData(self)
            try
                fid = fopen(self.fname(), 'r');
                self.data = fread(fid, [2048 2048], 'double');
                fclose(fid);
            catch
                warning('Pinhole:DataNotAvailable',[...
                    'Could not find previously calculated data for ' ...
                    'this pinhole.\n\nWill calculate it now. '...
                    'This may take some time. Interrupt with Ctrl+C\n']);
                self.calc_n_save();
            end
        end
        
        function calc(self)
            r = linspace(0, self.max_r,2048);
            z = linspace(0, self.max_z,2048);
            self.data = zeros(length(z),length(r));
            for r_ = 1:length(r)
                for z_=1:length(z)
                    self.data(r_,z_) = pinhole2(r(r_),z(z_),...
                        self.NA, self.n, ...
                        self.lambda, self.f_airy, self.f_res);
                end
            end
            %norm = pinhole(0,0,...
            %            self.NA, self.n, ...
            %            self.lambda, self.f_airy);
            %self.data=self.data/norm;
        end
        function calc_n_save(self)
            self.calc()
            fid = fopen(...
                (self.fname()), 'w');
            fwrite(fid, self.data, 'double');
            fclose(fid);
        end
        
        function fn = fname(self)
            fn = sprintf('pinhole-NA%g-n%g-lambda%g-fairy%g-f-res%g', ...
                    self.NA, self.n, self.lambda, self.f_airy, self.f_res);
        end
        
        function CEF = getIdx(self, i,j)
            if isempty(self.data)
                self.generateData()
            end
            idx = sub2ind(size(self.data), i, j);
            CEF = self.data(idx);
        end
        
        function I = get(self, r, z)
            R = abs(r);
            Z = abs(z);
            if R > self.max_r
                error('Pinhole:OutOFBounds','Requested position of r (%g) is out ouf the data size (max=%g).',r,self.max_r)
            end
            if Z > self.max_z
                error('Pinhole:OutOfBounds','Position of Z is out ouf the data size.')
            end
            
            I = self.getIdx(...
                round((R/self.max_r) * 2047) + 1, ...
                round((Z/self.max_z) * 2047) + 1);
        end
        function plotXY(self)
            if isempty(self.xy_data)
                x = linspace(-1.5, 1.5, 513);
                y = linspace(-1.5, 1.5, 513);
                self.xy_data = zeros(length(y),  length(x));
                P = ProgressBar(length(y));
                for y_ = 1 : length(y)
                    for x_ = 1 : length(x)
                        r = sqrt(x(x_)^2+y(y_)^2);
                        self.xy_data(y_, x_) = ...
                            self.get(r, 0) .* gaussint(x(x_), y(y_), 0,...
                            self.lambda, self.NA, self.f_res);
                        
                    end
                    P.progress();
                end
            end
            imagesc(self.xy_data);
        end
        function plotRZ(self)
            if isempty(self.rz_data)
                r = linspace(-1.5, 1.5, 513);
                z = linspace(-1.5, 1.5, 513);
                self.rz_data = zeros(length(z),  length(r));
                for z_ = 1 : length(z)
                    for r_ = 1 : length(r)
                        self.rz_data(z_, r_) = ...
                            self.get(r(r_), z(z_)) .* gaussint(r(r_), 0,...
                            z(z_), self.lambda, self.NA, self.f_res);
                    end
                end
            end
            imagesc(self.rz_data);
        end
        function cef = plotCEF(self)
            r = linspace(-1.5, 1.5, 513);
            z = linspace(-1.5, 1.5, 513);
            cef = zeros(length(z),  length(r));
            for z_ = 1 : length(z)
                for r_ = 1 : length(r)
                    cef(z_, r_) = ...
                        self.get(r(r_), z(z_));
                end
            end
            imagesc(cef);
        end
    end
end