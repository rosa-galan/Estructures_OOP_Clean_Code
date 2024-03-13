function [xnod,Tnod,Tmat] = buildSectionGeometry(type,d,h1,h2,t1,t2,t3)
% DO NOT MODIFY THIS FUNCTION
switch type
    case 'closed'
        % Nodal coordinates of the geometry points
        xg = @(d,h1,h2,t1,t2,t3,s) [
            (d+t1)/2	(h1+t3+t1*(h1-h2)/d/2)/2
            -(d+t2)/2	(h2+t3-t2*(h1-h2)/d/2)/2
            -(d+t2)/2	-(h2+t3-t2*(h1-h2)/d/2)/2
            (d+t1)/2	-(h1+t3+t1*(h1-h2)/d/2)/2
            (d-t1)/2	(h1-t3-t1*(h1-h2)/d/2)/2
            -(d-t2)/2	(h2-t3+t2*(h1-h2)/d/2)/2
            -(d-t2)/2	-(h2-t3+t2*(h1-h2)/d/2)/2
            (d-t1)/2	-(h1-t3-t1*(h1-h2)/d/2)/2
            (d-t1)/2	(h1+t3-t1*(h1-h2)/d/2)/2
            -(d-t2)/2	(h2+t3+t2*(h1-h2)/d/2)/2
            -(d+t2)/2	(h2-t3-t2*(h1-h2)/d/2)/2
            -(d+t2)/2	-(h2-t3-t2*(h1-h2)/d/2)/2
            -(d-t2)/2	-(h2+t3+t2*(h1-h2)/d/2)/2
            (d-t1)/2	-(h1+t3-t1*(h1-h2)/d/2)/2
            (d+t1)/2	-(h1-t3+t1*(h1-h2)/d/2)/2
            (d+t1)/2	(h1-t3+t1*(h1-h2)/d/2)/2
        ];
        % Nodal connectivities of geometry subdomains
        Tg = [
            1	9	5	16
            2	11	6	10
            3	13	7	12
            4	15	8	14
            5	9	10	6
            6	11	12	7
            7	13	14	8
            8	15	16	5
        ];
        % Load mesh data
        load('closed_section','xref','Tnod','Tmat','Tcon','Tele');
        %   xref: Nodal coordinates matrix for reference mesh
        %   Tnod: Nodal connectivities matrix
        %   Tmat: Material connectivities matrix
        %   Tcon: Contour elements connectivities
        %   Tele: Cell array of subdomains connectivities
        % Reference parameters (must match the ones used in the mesh)
        d_ = 1;
        h1_ = 1;
        h2_ = 1;
        t1_ = 0.1;
        t2_ = 0.1;
        t3_ = 0.1;
        s_ = 0;
    case 'open'
        % Nodal coordinates of the geometry points
        xg = @(d,h1,h2,t1,t2,t3,s) [
            (d+t1)/2	(h1+t3+t1*(h1-h2)/d/2)/2
            -(d+t2)/2	(h2+t3-t2*(h1-h2)/d/2)/2
            -(d+t2)/2	-(h2+t3-t2*(h1-h2)/d/2)/2
            (d+t1)/2	-(h1+t3+t1*(h1-h2)/d/2)/2
            (d-t1)/2	(h1-t3-t1*(h1-h2)/d/2)/2
            -(d-t2)/2	(h2-t3+t2*(h1-h2)/d/2)/2
            -(d-t2)/2	-(h2-t3+t2*(h1-h2)/d/2)/2
            (d-t1)/2	-(h1-t3-t1*(h1-h2)/d/2)/2
            (d-t1)/2	(h1+t3-t1*(h1-h2)/d/2)/2
            -(d-t2)/2	(h2+t3+t2*(h1-h2)/d/2)/2
            -(d+t2)/2	(h2-t3-t2*(h1-h2)/d/2)/2
            -(d+t2)/2	-(h2-t3-t2*(h1-h2)/d/2)/2
            -(d-t2)/2	-(h2+t3+t2*(h1-h2)/d/2)/2
            (d-t1)/2	-(h1+t3-t1*(h1-h2)/d/2)/2
            (d+t1)/2	-(h1-t3+t1*(h1-h2)/d/2)/2
            (d+t1)/2	(h1-t3+t1*(h1-h2)/d/2)/2
            (d-t1)/2	s
            (d-t1)/2	-s
            (d+t1)/2	-s
            (d+t1)/2	s
        ];
        % Nodal connectivities of geometry subdomains
        Tg = [
            1	9	5	16
            2	11	6	10
            3	13	7	12
            4	15	8	14
            5	9	10	6
            6	11	12	7
            7	13	14	8
            8	15	19	18
            17	20	16	5
        ];
        % Load mesh data
        load('open_section','xref','Tnod','Tmat','Tcon','Tele');
        %   xref: Nodal coordinates matrix for reference mesh
        %   Tnod: Nodal connectivities matrix
        %   Tmat: Material connectivities matrix
        %   Tcon: Contour elements connectivities
        %   Tele: Cell array of subdomains connectivities
        % Reference parameters (must match the ones used in the mesh)
        d_ = 1;
        h1_ = 1;
        h2_ = 1;
        t1_ = 0.1;
        t2_ = 0.1;
        t3_ = 0.1;
        s_ = h1_/10;    
end
% Geometry coordinates for reference mesh
xg_ = xg(d_,h1_,h2_,t1_,t2_,t3_,s_);
% Transform functions
inod_ = cell(size(Tg,1),1);
N_ = cell(size(Tg,1),1);
% Element vertices in natural coordinates
xi_ = [-1,1,1,-1]';
et_ = [-1,-1,1,1]';
% Loop through subdomains
for e = 1:size(Tg,1)
    % Get nodes contained in subdomain
    inod_{e} = unique(Tele{e}(:));
    % Get coordinates of subdomain with respect to subdomain centroid
    x = xref(inod_{e},1)-sum(xg_(Tg(e,:),1))/4;
    y = xref(inod_{e},2)-sum(xg_(Tg(e,:),2))/4;
    % Precomputations
    xh = sum(xi_.*xg_(Tg(e,:),1))/4;
    xv = sum(et_.*xg_(Tg(e,:),1))/4;
    xd = sum(xi_.*et_.*xg_(Tg(e,:),1))/4;
    yh = sum(xi_.*xg_(Tg(e,:),2))/4;
    yv = sum(et_.*xg_(Tg(e,:),2))/4;
    yd = sum(xi_.*et_.*xg_(Tg(e,:),2))/4;
    % Get natural coordinates for all points inside the subdomain 
    X0 = x*yv-y*xv;
    X1 = xh*yv-xv*yh+xd*y-yd*x;
    X2 = xd*yh-yd*xh;
    xi1 = X0./X1;
    if abs(X2)<1e-5
        xi = xi1;
    else
        xi2 = X1./X2;
        xi = xi2/2.*(1+sqrt(1-4*xi1./xi2));
        icond = abs(xi)>1;
        xi(icond) = xi2(icond)/2.*(1-sqrt(1-4*xi1(icond)./xi2(icond)));
    end
    Y0 = y*xh-x*yh;
    Y1 = yv*xh-yh*xv+yd*x-xd*y;
    Y2 = yd*xv-xd*yv;
    et1 = Y0./Y1;
    if abs(Y2)<1e-5
        et = et1;
    else
        et2 = Y1./Y2;
        et = et2/2.*(1+sqrt(1-4*et1./et2));
        icond = abs(et)>1;
        et(icond) = et2(icond)/2.*(1-sqrt(1-4*et1(icond)./et2(icond)));
    end
    % Get the transformation function for all points inside the subdomain
    N_{e} = [
        (1+xi_(1)*xi).*(1+et_(1)*et),...
        (1+xi_(2)*xi).*(1+et_(2)*et),...
        (1+xi_(3)*xi).*(1+et_(3)*et),...
        (1+xi_(4)*xi).*(1+et_(4)*et)
    ]/4;
end
% Update coordinates
l23 = sqrt(d^2+((h1-h2)/2)^2);
xnew = xg(d,h1,h2,t1,t2,t3*d/l23,0);
xnod = zeros(size(xref));
for e = 1:size(Tg,1)
    xnod(inod_{e},1) = N_{e}*xnew(Tg(e,:),1);
    xnod(inod_{e},2) = N_{e}*xnew(Tg(e,:),2);
end
end