function [I, comp] = pinhole2( R, Z, NA, n, wavelength, varargin)
% This file generates the effect of a pinhole.


    f_airy = 1;
    if nargin >= 6
        f_airy = varargin{1};
    end
    f_res = 1;
    if nargin == 7
        f_res = varargin{2};
    end
    r_ph = f_airy * 1.22 * wavelength / NA;
    cone_angle = asin( NA / n);
    R0 = 1/f_res * wavelength/(2 * NA);
    
    
    d = abs(R);
    z = abs(Z);
    
    Rz = R0 + z.*(tan(cone_angle));
    
    if r_ph+Rz <= d 
            I = 0;
            comp=0;
            return
    end
    ratio = r_ph/Rz;
    
    
    % ratio >= 1: Pinhole larger than Rz
    x_is = (r_ph^2+d^2-Rz^2)/(2*d);        
    if ratio >= 1
        if d+Rz <= r_ph
            I = 1;
            comp = -1;
        elseif x_is < d
            
            a = l_circle_segment(r_ph, x_is) + l_circle_segment(Rz, d-x_is);
            I = a/(pi*Rz^2);
            comp = -2;
        else
            
            a = pi * Rz^2 + l_circle_segment(r_ph, x_is) - l_circle_segment(Rz, x_is-d);
            I = a/(pi*Rz^2);
            comp=x_is-d;
        end     
        return
    end
    
    % Rz larger than pinhole...
    if d+r_ph <= Rz % Rz covers pinhole completely
        I = ratio^2;
        comp = 1;
        return
    end
    
    if x_is >= 0
        a = l_circle_segment(r_ph, x_is) + l_circle_segment(Rz, d-x_is);
        I = a/(pi * Rz^2);
        comp = 2;
        return
    end
    
    a = pi * r_ph^2 - l_circle_segment(r_ph, x_is)+ l_circle_segment(Rz, d-x_is);
    
    I = a / (pi * Rz^2);
    comp = 3;
end

function a = l_circle_segment(r,x)
    x = abs(x);
    r = abs(r);
    if x > r
        fprintf('Will compute complex number, since r=%g > x%g\n', r, x);
    end
    y = sqrt(r^2-x^2);
    alpha = asin(y/r);
    a = r^2*pi * alpha / pi;
    a = a - y*x;
end
