function [x,y,z] = simulateDiffusion3D_2(scan, D, n, Dt, varargin)

    offsetx = min(scan.xdata_lin);
    offsety = min(scan.ydata_lin);
    
    %% Check for sticky object
    sticky_object = false;
    if nargin > 4
        sticky_object = true;
        sticky_pos = varargin{1};
        sticky_pos(1) = sticky_pos(1) - offsetx;
        sticky_pos(2) = sticky_pos(2) - offsety;
        st_r = sticky_pos(3);
        st_x = sticky_pos(1);
        st_y = sticky_pos(2);
    end
    warn_count = 0;
    
    % Defining the stickyness probability 
    sticky_prob = .1;
    if nargin > 5
        sticky_prob = varargin{2};
    end
    
    %%
    x = zeros(n,1);
    y = zeros(n,1);
      
    try
        F = griddedInterpolant(scan.xdata_grid, scan.ydata_grid, ...
            scan.zdata_grid, 'cubic');
    catch ME
        if strcmp(ME.identifier, ...
                'MATLAB:griddedInterpolant:NdgridNotMeshgrid2DErrId')
            F = griddedInterpolant(scan.xdata_grid', scan.ydata_grid', ...
            scan.zdata_grid, 'cubic');
        else
            rethrow(ME); 
        end
    end
    
    %% Interpolate scan (tmp_scan)
    tmp_scan = scan.interpolate(20, 'cubic');
    tmp_dx = tmp_scan.stepx;
    tmp_dy = tmp_scan.stepy;
            
    sigma = sqrt(2 * D * Dt);
        
    % Random x and y positions:
    pos_rand = normrnd( 0, sigma, n, 2 );
    
    %% Compute start position 
    r = rand(2,1);  
    
    xpos = r(1) * ...
    (max(scan.xdata_lin) - min(scan.xdata_lin) - tmp_scan.stepx)...
    + min(scan.xdata_lin);
    
    ypos = r(2) * ...
    (max(scan.ydata_lin) - min(scan.ydata_lin)- tmp_scan.stepy)...
    + min(scan.ydata_lin);
    
    if sticky_object
        while (xpos-st_x)^2/st_r^2 + (ypos-st_y)^2/st_r^2 <= 1
            r = rand(2,1);
            
            xpos = r(1) * ...
            (max(scan.xdata_lin) - min(scan.xdata_lin) - tmp_scan.stepx)...
            + min(scan.xdata_lin);
            
            ypos = r(2) * ...
            (max(scan.ydata_lin) - min(scan.ydata_lin) - tmp_scan.stepy)...
            + min(scan.ydata_lin);
        end
    end

    x(1) = xpos;
    y(1) = ypos;

    [mgx, mgy] = gradient(tmp_scan.zdata_grid, tmp_dx, tmp_dy);
    
    scan_max_x = max(scan.xdata_lin);
    scan_max_y = max(scan.ydata_lin);
    
    lcase = 0;

    is_sticky = true; 
    
    for i = 1 : n

        row = 1 + round((x(i) - offsetx) / scan_max_x * (tmp_scan.xpx - 1)); 
        col = 1 + round((y(i) - offsety) / scan_max_y * (tmp_scan.ypx - 1)); 
        try
            gx = mgx(row, col);
            gy = mgy(row, col);
        catch E
            fprintf('Row: %g, Col: %g, x: %g, y:%g, last case: %g, i: %g',row,col, x(i), y(i), lcase, i);
            rethrow(E);
        end
        lcase = 0;
        %% Compute x and check whether it is still within the simulated area
        
        x(i+1) = x(i) + cos(gx) * pos_rand(i,1) ;
        if x(i+1) < offsetx             % left the scan area at minimum
            x(i+1) = scan_max_x + (x(i+1) - offsetx);
            lcase = 1;
        elseif x(i+1) >= scan_max_x     % left the scan area at maximum
            x(i+1) = offsetx + (x(i+1) - scan_max_x);
            lcase = 2;
        end
        
        
        %% Compute y and check whether it is still within the simulated area
        
        y(i+1) = y(i) + cos(gy) * pos_rand(i,2);
        
        if y(i+1) < offsety             % left the scan area at minimum
            y(i+1) = scan_max_y + (y(i+1) - offsety);
            lcase = 3;
        elseif y(i+1) > scan_max_y      % left the scan area at maximum
            y(i+1) = offsety + (y(i+1) - scan_max_y);
            lcase = 4;
        end

        if sticky_object
            if is_sticky
                if rand() < sticky_prob
                    is_sticky = false;                    
                else
                    x(i+1) = x(i);
                    y(i+1) = y(i);
                end
            end
            
            if ~is_sticky && lcase == 0
                X0 = x(i) - st_x;
                X1 = x(i+1) - st_x;
                Y0 = y(i) - st_y;
                Y1 = y(i+1) - st_y;
            
                D = st_r^2 * ((X1-X0)^2 + (Y1-Y0)^2)...
                    - ( X0 * (Y1-Y0) - Y0 * (X1-X0) )^2;
    
                st_r2 = st_r;
                while abs(D) <= 1e-3
                    st_r2 = 1.05*st_r2;
                    D = st_r2^2 * ((X1-X0)^2 + (Y1-Y0)^2)...
                        - ( X0 * (Y1-Y0) - Y0 * (X1-X0) )^2;
                end    
                
                if D >= 0
                    t1 = (-X0 * (X1-X0) - Y0 * (Y1-Y0) + sqrt(D)) ...
                        / ((X1-X0)^2 + (Y1-Y0)^2);
                    t2 = (-X0 * (X1-X0) - Y0 * (Y1-Y0) - sqrt(D)) ...
                        / ((X1-X0)^2 + (Y1-Y0)^2);
    
                    if within_or_equal(t1, 0, 1) && within_or_equal(t2, 0, 1)
                        is_sticky = true;
                        t = min([t1, t2]);
                    elseif within_or_equal(t1, 0, 1)
                        t = t1;
                        is_sticky = true;
                    elseif within_or_equal(t2, 0, 1)
                        t = t2;
                        is_sticky = true;
                    end
                
                    if is_sticky
                        x(i+1) = x(i) + t * (x(i+1)-x(i));
                        y(i+1) = y(i) + t * (y(i+1)-y(i));
                    end
                end
                
                if sticky_object && ((x(i+1)-st_x)/st_r)^2 + ((y(i+1)-st_y)/st_r)^2 < 1
                    warn_count = warn_count + 1;
                    is_sticky = true;
                    alpha = atan(Y1/X1);
                    x1 = abs(st_r*cos(alpha)) * sign(X1);
                    y1 = abs(st_r*sin(alpha)) * sign(Y1);
                    x(i+1) = x1 + st_x;
                    y(i+1) = y1 + st_y;
                    if abs(x(i) - x(i+1)) > 0.1 || abs(y(i) - y(i+1)) > 0.1
                        fprintf('Adjusted x-pos from %g to %g and y-pos from %g to %g.\n', X1,x1,Y1,y1);
                        fprintf('Long distance detected (from X: %g to %g, Y %g to %g)\n',x(i),x(i+1), y(i), y(i+1));
                        fprintf('angle was %g\n', alpha*360/(2*pi));
                    end
                end
            end
        end
    end
    
    z = F(x,y);
    x = x * scan_max_x/(scan_max_x - offsetx); 
    y = y * scan_max_y/(scan_max_x - offsetx); 
%     fprintf('There were %g warnings\n.', warn_count);
end

function r = within_or_equal (value, v1, v2)
    low = min([v1,v2]);
    high = max([v1,v2]);
    r = low <= value & value <= high;
end