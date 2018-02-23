%... The MatMol Group (2016)
    function [x y ksi_x ksi_y eta_x eta_y ksi_xx ksi_yy eta_xx eta_yy] = actual_grid(corners,nksi,neta,ksi,eta)
%...
%... This function computes the spatial grid
%...
%...   inputs : 
%...           * corners is a (4 x 2) matrix containg the coordinates of the
%...             vertices of the quadrilateral
%...           * nksi et neta are the number of points along the S1-S2 and
%...             S1-S4, respectively
%...           * ksi and eta are the vectors of points along the ksi and eta
%...             axis, respectively.
%...
%... The corners will be sorted so that S1 denotes the left-most corner, i.e.,
%... the one with the smallest x-coordinate. If two alternative points exist,
%... the lower one, i.e. the one with a smaller y-coordinate, is selected.
%... The next corners are selected according to the slope of the segment
%... S1Si. S2 corresponds to the smallest slope, S4 to the largest.
%... These slopes are between -pi/2 et +pi/2 (this is a direct consequence 
%... of the choice of S1). 
%... The following is an example:
%...
%...
%...           S4 o
%...
%...
%...       S1 o
%...
%...                               S3 o
%...
%...                   S2 o
%...
%...
%...   outputs : 
%...           * x and y are matrices of dimension (neta x nksi) with the x
%...             and y coordinates of the grid. These matrices are useful
%...             for graphing the solution.
%...             S1-S2 corresponds to the ksi-axis and S1-S4 the eta-axis.
%...             Numbering is from left to right, bottom-up, starting from
%...             S1 and ending in S3.
%...
%...           * ksi_x, ksi_y, eta_x et eta_y are matrices of dimension
%...             (nksi x neta) giving the following derivatives evaluated
%...             in the grid points
%...
%...            d(ksi)   d(ksi)   d(eta)   d(eta)
%...            ------ , ------ , ------ , ------      [1]
%...              dx       dy       dx       dy
%...
%...           * ksi_x, ksi_y, eta_x et eta_y are matrices of dimension
%...             (nksi x neta) giving the following derivatives evaluated
%...             in the grid points
%...
%...            d2(ksi)   d2(ksi)   d2(eta)   d2(eta)
%...            ------- , ------- , ------- , -------      [2]
%...              dx2       dy2       dx2       dy2
%...
%...
%... sort the corners
%...
    xmin = min(corners(:,1));
    ind = find(corners(:,1) == xmin);
    if length(ind) == 2
        Sselect = corners(ind,:);
        ymin = min(Sselect(:,2));
        indbis = find(Sselect(:,2) == ymin);
        S1 = Sselect(indbis,:);
    else
        S1 = corners(ind,:);
    end
    S = setdiff(corners,S1,'rows');
    test = sortrows([(S(:,2)-S1(2))./(S(:,1)-S1(1)) S]);
    S2 = test(1,2:3);
    S3 = test(2,2:3);
    S4 = test(3,2:3);
%...
%...compute the actual grid coordinates (x,y)
%...   
    coord = zeros(nksi*neta,2);
    for k = 0:neta-1
        coordl = S1 + k*(S4-S1)/(neta-1);
        coordr = S2 + k*(S3-S2)/(neta-1);
        for j = 0:nksi-1
            coord(1+k*nksi+j,1:2) = coordl + j*(coordr-coordl)/(nksi-1);
        end
    end
    x = zeros(neta,nksi);
    y = zeros(neta,nksi);
    for i = 1:neta
        x(neta+1-i,1:nksi) = coord((i-1)*nksi+1:i*nksi,1);
        y(neta+1-i,1:nksi) = coord((i-1)*nksi+1:i*nksi,2);
    end
%...
%... the change of coordinates can be expressed as
%...
%...       x = a(0) + a(1)*ksi + a(2)*eta + a(3)*ksi*eta
%...
%...       y = b(0) + b(1)*ksi + b(2)*eta + b(3)*ksi*eta
%...
%... where the coefficients a(i) and b(i) depends on the corner definition.
%... In the following, a(0) and b(0) are not evaluated as they are not used
%... for computing the derivatives [1]
%...
    a = zeros(1,3);
    b = zeros(1,3);
    a(1) = S2(1) - S1(1);
    a(2) = S4(1) - S1(1);
    a(3) = S3(1) + S1(1) - S2(1) - S4(1);
    b(1) = S2(2) - S1(2);
    b(2) = S4(2) - S1(2);
    b(3) = S3(2) + S1(2) - S2(2) - S4(2);
%...
%... compute derivatives
%...
    d = zeros(nksi*neta,1);
    ksi_x = zeros(nksi*neta,1);
    ksi_y = zeros(nksi*neta,1);
    eta_x = zeros(nksi*neta,1);
    eta_y = zeros(nksi*neta,1);
    d0 = a(1)*b(2)-a(2)*b(1);
    d1 = a(1)*b(3)-a(3)*b(1);
    d2 = a(3)*b(2)-a(2)*b(3);
%...
    for i = 1:neta
        d((i-1)*nksi+1:(i-1)*nksi+nksi,1) = d0 + d1*ksi + d2*eta(i);
        ksi_x((i-1)*nksi+1:(i-1)*nksi+nksi,1) =  (b(2)+b(3)*ksi)./d((i-1)*nksi+1:(i-1)*nksi+nksi,1);
        ksi_y((i-1)*nksi+1:(i-1)*nksi+nksi,1) = -(a(2)+a(3)*ksi)./d((i-1)*nksi+1:(i-1)*nksi+nksi,1);
        ksi_xx((i-1)*nksi+1:(i-1)*nksi+nksi,1) =  2*(b(2)+b(3)*ksi)*(b(1)+b(3)*eta(i))*(a(3)*b(2)-a(2)*b(3))./(d((i-1)*nksi+1:(i-1)*nksi+nksi,1).^3);
        ksi_yy((i-1)*nksi+1:(i-1)*nksi+nksi,1) =  2*(a(2)+a(3)*ksi)*(a(1)+a(3)*eta(i))*(a(3)*b(2)-a(2)*b(3))./(d((i-1)*nksi+1:(i-1)*nksi+nksi,1).^3);
    end
%...
    for i = 1:nksi
        eta_x(i:nksi:(neta-1)*nksi+i,1) = -(b(1)+b(3)*eta)./d(i:nksi:(neta-1)*nksi+i,1);
        eta_y(i:nksi:(neta-1)*nksi+i,1) =  (a(1)+a(3)*eta)./d(i:nksi:(neta-1)*nksi+i,1);
        eta_xx(i:nksi:(neta-1)*nksi+i,1) = 2*(b(1)+b(3)*eta)*(b(2)+b(3)*ksi(i))*(a(1)*b(3)-a(3)*b(1))./(d(i:nksi:(neta-1)*nksi+i,1).^3);
        eta_yy(i:nksi:(neta-1)*nksi+i,1) = 2*(a(1)+a(3)*eta)*(a(2)+a(3)*ksi(i))*(a(1)*b(3)-a(3)*b(1))./(d(i:nksi:(neta-1)*nksi+i,1).^3);
    end
