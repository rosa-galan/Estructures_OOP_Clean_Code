function [Xp,Yp,sig,tau] = sectionAnalysis(type,d,h1,h2,t1,t2,t3,E,G,xref,Q,Mb,Mt)
% Inputs:
%  - type: 'open' for open section, 'closed' for closed section.
%  - d, h1, h2, t1, t2, t3: geometric parameters of the section.
%  - E, G: Young and shear modulus of the material.
%  - xref: x-coordinate of the point where the force Q is applied/torsion
%          moments Mt are taken.
%  - Q: magnitude of the shear force.
%  - Mb: magnitude of the bending moment.
%  - Mt: magnitude of the torsion moment.
% Outputs:
%  - Xp, Yp: x and y coordinates vectors of the points in the section where
%            stress has been evaluated.
%  - sig: vector with the normal stress at the corresponding points.
%  - tau: vector with the tangential stress at the corresponding points.
%  Note: sig(i) and tau(i) provide the stress components evaluated at point
%        with coordinates Xp(i) and Yp(i).
%  Note: All magnitudes given in international system units.

%% Build geometry
[xnod,Tnod,Tmat] = buildSectionGeometry(type,d,h1,h2,t1,t2,t3);
mat = [E, G];

%% Dimensions
nele = size(Tnod,1);
nnode = size(Tnod,2);
nnod = size(xnod,1);

%% Find boundary
nodCon = [1:nnode;[2:nnode,1]]';
Tbnd = reshape(Tnod(:,nodCon')',size(nodCon,2),[]);
Tmb = repelem(Tmat,nnode,1);
[~,ia,ic] = unique(sort(Tbnd,1)','rows');
Tbnd = Tbnd(:,ia(accumarray(ic,1)==1))';
Tmb = Tmb(ia(accumarray(ic,1)==1),1);
% Sort boundary elements into sets of closed curves
i = 1;
Tnb_{1} = Tbnd(1,:);
Tmb_{1} = Tmb(1,1);
Tbnd = Tbnd(2:end,:);
Tmb = Tmb(2:end,1);
while size(Tbnd,1)>0
    ind_ = find(Tbnd(:,1)==Tnb_{i}(end,2));
    if isempty(ind_)
        i = i+1;
        Tnb_{i} = Tbnd(1,:);
        Tmb_{i} = Tmb(1,1);
        Tbnd = Tbnd(2:end,:);
        Tmb = Tmb(2:end,1);
    else
        Tnb_{i} = [Tnb_{i};Tbnd(ind_,:)];
        Tmb_{i} = [Tmb_{i};Tmb(ind_,1)];
        Tbnd = Tbnd(setdiff(1:size(Tbnd,1),ind_),:);
        Tmb = Tmb(setdiff(1:size(Tmb,1),ind_),:);
    end
end
Tbnd = cat(1,Tnb_{:});
Tmb = cat(1,Tmb_{:});
fdof = 2:nnod;
idof = repmat(Tnod,1,nnode,1);
jdof = repelem(Tnod,1,nnode,1);
idofb = repmat(Tbnd,1,2,1);
jdofb = repelem(Tbnd,1,2,1);

%% Material properties
Emat = mat(Tmat,1);
Gmat = mat(Tmat,2);
Gmatb = mat(Tmb,2);

%% Shape functions and natural coordinates derivatives

% Linear elements (boundary)
xgb = [-1, 1]/sqrt(3);
wgb = [ 1, 1];
Nb(1,:,1) = (1-xgb)/2;
Nb(1,:,2) = (1+xgb)/2;
switch nnode
    case 3 % Triangle elements
        xg = [ 1, 1; 4, 1; 1, 4]/6;
        wg = [    1,    1,    1]/6;
        N(1,:,1,1,1) = 1-xg(:,1)-xg(:,2);
        N(1,:,1,1,2) = xg(:,1);
        N(1,:,1,1,3) = xg(:,2);
        dN(1,1,1,:,1) = -1;
        dN(1,1,1,:,2) =  1;
        dN(1,1,1,:,3) =  0;
        dN(1,2,1,:,1) = -1;
        dN(1,2,1,:,2) =  0;
        dN(1,2,1,:,3) =  1;
    case 4 % Quadrilateral elments
        xg = [-1,-1; 1,-1; 1, 1;-1, 1]/sqrt(3);
        wg = [    1,    1,    1,    1];
        N(1,:,1,1,1) = (1-xg(:,1)).*(1-xg(:,2))/4;
        N(1,:,1,1,2) = (1+xg(:,1)).*(1-xg(:,2))/4;
        N(1,:,1,1,3) = (1+xg(:,1)).*(1+xg(:,2))/4;
        N(1,:,1,1,4) = (1-xg(:,1)).*(1+xg(:,2))/4;
        dN(1,1,1,:,1) = -(1-xg(:,2))/4;
        dN(1,1,1,:,2) =  (1-xg(:,2))/4;
        dN(1,1,1,:,3) =  (1+xg(:,2))/4;
        dN(1,1,1,:,4) = -(1+xg(:,2))/4;
        dN(1,2,1,:,1) = -(1-xg(:,1))/4;
        dN(1,2,1,:,2) = -(1+xg(:,1))/4;
        dN(1,2,1,:,3) =  (1+xg(:,1))/4;
        dN(1,2,1,:,4) =  (1-xg(:,1))/4;
end
Yn = xnod(:,1);
Zn = xnod(:,2);
yel = reshape(Yn(Tnod),nele,1,1,1,nnode); 
zel = reshape(Zn(Tnod),nele,1,1,1,nnode); 
J(:,1,1,:) = sum(yel.*dN(:,1,1,:,:),5);
J(:,1,2,:) = sum(zel.*dN(:,1,1,:,:),5);
J(:,2,1,:) = sum(yel.*dN(:,2,1,:,:),5);
J(:,2,2,:) = sum(zel.*dN(:,2,1,:,:),5);
detJ = squeeze(J(:,1,1,:).*J(:,2,2,:)-J(:,1,2,:).*J(:,2,1,:));

%% Matrices
A1_el = sum(wg.*detJ.*permute(N,[1,2,5,4,3]).*permute(N,[1,2,3,5,4]),2);
A1    = sparse(idof(:),jdof(:),A1_el(:),nnod,nnod);
EA_el = sum(Emat.*wg.*detJ.*permute(N,[1,2,5,4,3]).*permute(N,[1,2,3,5,4]),2);
EA    = sparse(idof(:),jdof(:),EA_el(:),nnod,nnod);
GA_el = sum(Gmat.*wg.*detJ.*permute(N,[1,2,5,4,3]).*permute(N,[1,2,3,5,4]),2);
GA    = sparse(idof(:),jdof(:),GA_el(:),nnod,nnod);

%% Average properties
A = full(sum(A1,'all'));
E = full(sum(EA,'all'))/A;
G = full(sum(GA,'all'))/A;

%% Centroid
yG = full(sum(A1*Yn,1))/A;
zG = full(sum(A1*Zn,1))/A;

%% Neutral centre
y0 = full(sum(EA*Yn,1))/E/A;
z0 = full(sum(EA*Zn,1))/E/A;
Yp = Yn-y0;
Zp = Zn-z0;

%% Principal axes and section inertias
Iyp = full(Zp'*EA*Zp)/E;
Izp = full(Yp'*EA*Yp)/E;
Iyzp = full(Zp'*EA*Yp)/E;
if abs(Iyzp)<1e-8
    sa = 0;
    ca = 1;
else
    ca = cos(atan2(-2*Iyzp,Iyp-Izp)/2);
    sa = sin(atan2(-2*Iyzp,Iyp-Izp)/2);
end
Yloc = ca*Yp + sa*Zp;
Zloc = ca*Zp - sa*Yp;
Iy = ca^2*Iyp + sa^2*Izp - 2*ca*sa*Iyzp;
Iz = ca^2*Izp + sa^2*Iyp + 2*ca*sa*Iyzp;

%% Local derivatives
yloc = reshape(Yloc(Tnod),nele,1,1,1,nnode); 
zloc = reshape(Zloc(Tnod),nele,1,1,1,nnode); 
Jloc(:,1,1,:) = sum(yloc.*dN(:,1,1,:,:),5);
Jloc(:,1,2,:) = sum(zloc.*dN(:,1,1,:,:),5);
Jloc(:,2,1,:) = sum(yloc.*dN(:,2,1,:,:),5);
Jloc(:,2,2,:) = sum(zloc.*dN(:,2,1,:,:),5);
detJloc = Jloc(:,1,1,:).*Jloc(:,2,2,:)-Jloc(:,1,2,:).*Jloc(:,2,1,:);
dNloc(:,1,1,:,:) = (Jloc(:,2,2,:).*dN(:,1,1,:,:) - Jloc(:,1,2,:).*dN(:,2,1,:,:))./detJloc;
dNloc(:,2,1,:,:) = (Jloc(:,1,1,:).*dN(:,2,1,:,:) - Jloc(:,2,1,:).*dN(:,1,1,:,:))./detJloc;
Ly = Yloc(Tbnd(:,2))-Yloc(Tbnd(:,1));
Lz = Zloc(Tbnd(:,2))-Zloc(Tbnd(:,1));

%% Matrices involving local derivatives
dxGAdx_el = sum(Gmat.*wg.*detJ.*sum(permute(dNloc,[1,4,5,3,2]).*permute(dNloc,[1,4,3,5,2]),5),2);
dxGAdx = sparse(idof(:),jdof(:),dxGAdx_el(:),nnod,nnod);
dxGA_el = sum(Gmat.*wg.*detJ.*permute(dNloc(:,1,:,:,:),[1,4,5,3,2]).*permute(N,[1,2,3,5,4]),2);
dxGA = sparse(idof(:),jdof(:),dxGA_el(:),nnod,nnod);
dyGA_el = sum(Gmat.*wg.*detJ.*permute(dNloc(:,2,:,:,:),[1,4,5,3,2]).*permute(N,[1,2,3,5,4]),2);
dyGA = sparse(idof(:),jdof(:),dyGA_el(:),nnod,nnod);
GLx_el = squeeze(Gmatb.*Ly.*sum(wgb/2.*Nb.*permute(Nb,[1,2,4,3]),2));
GLx = sparse(idofb(:),jdofb(:),GLx_el(:),nnod,nnod);
GLy_el = squeeze(Gmatb.*Lz.*sum(wgb/2.*Nb.*permute(Nb,[1,2,4,3]),2));
GLy = sparse(idofb(:),jdofb(:),GLy_el(:),nnod,nnod);

%% Shear distribution

% Shear force in the x-direction
Qy = 1;
phi_y(fdof,1) = dxGAdx(fdof,fdof)\(EA(fdof,:)*Yloc)*Qy/E/Iz;
Mx_y = phi_y'*(dyGA*Yloc-dxGA*Zloc);
ky = full(Qy^2/G/A/(phi_y'*dxGAdx*phi_y));
txy_y = full(Gmat.*squeeze(sum(dNloc(:,1,1,:,:).*permute(phi_y(Tnod),[1,3,4,5,2]),5)));
txz_y = full(Gmat.*squeeze(sum(dNloc(:,2,1,:,:).*permute(phi_y(Tnod),[1,3,4,5,2]),5)));

% Shear force in the y-direction
Qz = 1;
phi_z(fdof,1) = dxGAdx(fdof,fdof)\(EA(fdof,:)*Zloc)*Qz/E/Iy;
Mx_z = phi_z'*(dyGA*Yloc-dxGA*Zloc);
kz = full(Qz^2/G/A/(phi_z'*dxGAdx*phi_z));
txy_z = full(Gmat.*squeeze(sum(dNloc(:,1,1,:,:).*permute(phi_z(Tnod),[1,3,4,5,2]),5)));
txz_z = full(Gmat.*squeeze(sum(dNloc(:,2,1,:,:).*permute(phi_z(Tnod),[1,3,4,5,2]),5)));

% Shear centre
yc = full(Mx_z)/Qz;
zc = full(Mx_y)/Qy;
Yc = Yloc-yc;
Zc = Zloc-zc;

% Pure torsion moment
dtdx = 1;
phi_t(fdof,1) = dxGAdx(fdof,fdof)\(GLx(fdof,:)*Yc+GLy(fdof,:)*Zc);
phi_t = phi_t - sum(A1*phi_t)/A;
Mx_t = (phi_t'*(dyGA*Yc-dxGA*Zc)+Yc'*GA*Yc+Zc'*GA*Zc)*dtdx;
Jc = full(Mx_t)/G/dtdx;
txy_t = full(Gmat/Jc/G.*(squeeze(sum(dNloc(:,1,1,:,:).*permute(phi_t(Tnod),[1,3,4,5,2]),5))-squeeze(sum(N.*permute(Zc(Tnod),[1,3,4,5,2]),5))));
txz_t = full(Gmat/Jc/G.*(squeeze(sum(dNloc(:,2,1,:,:).*permute(phi_t(Tnod),[1,3,4,5,2]),5))+squeeze(sum(N.*permute(Yc(Tnod),[1,3,4,5,2]),5))));

%% Output

% Section area
m.A = A;

% Axial + Bending properties
m.E = E; % Effective Young modulus
m.x0 = [ y0, z0]; m.x0(abs(m.x0)<1e-8) = 0;       % Neutral axis position  
m.alpha = atan2(sa,ca); m.alpha(abs(m.alpha)<1e-8) = 0; % Principal axes orientation
m.Iy = Iy; % Inertia about first principal axis
m.Iz = Iz; % Inertia about secont principal axis

% Shear + Torsion
m.G = G; % Effective shear modulus
m.xc = [ y0+ca*yc-sa*zc, z0+ca*zc+sa*yc]; m.xc(abs(m.xc)<1e-8) = 0; % Shear center position
m.ky = ky; % Shear correction parameter in y
m.kz = kz; % Shear correction parameter in z
m.Jc = Jc; % Torsional inertia about shear center

%% Stress

% Precomputations
y_ge = squeeze(sum(N.*yel,5));
z_ge = squeeze(sum(N.*zel,5));
yloc_el = reshape(Yloc(Tnod),nele,1,1,1,nnode); 
zloc_el = reshape(Zloc(Tnod),nele,1,1,1,nnode); 
yloc_ge = squeeze(sum(N.*yloc_el,5));
zloc_ge = squeeze(sum(N.*zloc_el,5));
E_ge = squeeze(sum(N.*repmat(Emat,1,1,1,1,nnode),5));

% Input loads
Fy = 0;
Fz = Q;
My = Mb;
Mz = 0;
Mx = Mt;

% Local shear forces
Fy_ = ca*Fy+sa*Fz;
Fz_ = ca*Fz-sa*Fy;

% Local moments about the shear centre
My_ = ca*My+sa*Mz;
Mz_ = ca*Mz-sa*My;
Mx_ = Mx+m.xc(2)*Fy-(m.xc(1)-xref)*Fz;

% Bending stress
sig_ge = E_ge(:)/E.*(zloc_ge(:)*My_/Iy-yloc_ge(:)*Mz_/Iz);
% Interpolation function
sig_fun = scatteredInterpolant(y_ge(:),z_ge(:),sig_ge,'linear','linear');

% Shear stress caused by shear load
%   In local coordinates
txyQ_ = txy_y(:)*Fy_+txy_z(:)*Fz_;
txzQ_ = txz_y(:)*Fy_+txz_z(:)*Fz_;
%   In global coordinates
txyQ_ge = ca*txyQ_-sa*txzQ_;
txzQ_ge = ca*txzQ_+sa*txyQ_;
% Interpolation function
tauQ_x_fun = scatteredInterpolant(y_ge(:),z_ge(:),txyQ_ge,'linear','linear');
tauQ_y_fun = scatteredInterpolant(y_ge(:),z_ge(:),txzQ_ge,'linear','linear');

% Shear stress caused by torque
%   In local coordinates
txyT_ = txy_t(:)*Mx_;
txzT_ = txz_t(:)*Mx_;
%   In global coordinates
txyT_ge = ca*txyT_-sa*txzT_;
txzT_ge = ca*txzT_+sa*txyT_;
% Interpolation function
tauT_x_fun = scatteredInterpolant(y_ge(:),z_ge(:),txyT_ge,'linear','linear');
tauT_y_fun = scatteredInterpolant(y_ge(:),z_ge(:),txzT_ge,'linear','linear');

%% Plot results

% Set scale value
scale = 0.2*min([max(yel(:))-min(yel(:)),max(zel(:))-min(zel(:))]);

% Initialize colormap
cmap = jet(255);

% Get evaluation points
ny = 15;
nx = 20;
nt = 6;
xt = (1/nt/2:1/nt:1-1/nt/2)-1/2;
Xp = [];
Yp = [];
for i = 1:nt
    % Vertex 1
    Xp = [Xp,d/2+xt*t1];
    Yp = [Yp,h1/2+xt(i)*t3+xt*t1*(h1-h2)/d/2];
    % Tram 1-2
    x1 =  (d-t1)/2;          
    x2 = -(d-t2)/2;
    Xp = [Xp,x1:(x2-x1)/nx:x2]; 
    y1 = h1/2+xt(i)*t3-t1*(h1-h2)/d/4;
    y2 = h2/2+xt(i)*t3+t2*(h1-h2)/d/4;
    Yp = [Yp,y1:(y2-y1)/nx:y2];
    % Vertex 2
    Xp = [Xp,-d/2-xt*t2];
    Yp = [Yp,h2/2+xt(i)*t3-xt*t2*(h1-h2)/d/2];
    % Tram 2-3
    x1 = -d/2-xt(i)*t2;
    Xp = [Xp,x1*ones(1,ny+1)]; 
    y1 =  (h2-t3)/2; 
    y2 = -(h2-t3)/2;
    Yp = [Yp,y1:(y2-y1)/ny:y2];
    % Vertex 3
    Xp = [Xp,-d/2-xt*t2];
    Yp = [Yp,-h2/2-xt(i)*t3+xt*t2*(h1-h2)/d/2];
    % Tram 3-4
    x1 = -(d-t2)/2;         
    x2 =  (d-t1)/2;         
    Xp = [Xp,x1:(x2-x1)/nx:x2]; 
    y1 = -(h2/2+xt(i)*t3+t2*(h1-h2)/d/4);
    y2 = -(h1/2+xt(i)*t3-t1*(h1-h2)/d/4);
    Yp = [Yp,y1:(y2-y1)/nx:y2];
    % Vertex 4
    Xp = [Xp,d/2+xt*t1];
    Yp = [Yp,-h1/2-xt(i)*t3-xt*t1*(h1-h2)/d/2];
    % Tram 4-1
    x1 =  d/2+xt(i)*t1;
    Xp = [Xp,x1*ones(1,ny+1)]; 
    y1 = -(h1-t3)/2; 
    y2 =  (h1-t3)/2;
    Yp = [Yp,y1:(y2-y1)/ny:y2];
end
Xp = Xp';
Yp = Yp';
sig = sig_fun(Xp,Yp);
tauQ = [tauQ_x_fun(Xp,Yp),tauQ_y_fun(Xp,Yp)];
tauT = [tauT_x_fun(Xp,Yp),tauT_y_fun(Xp,Yp)];
tau = sqrt(sum((tauQ+tauT).^2,2));

% Plot sig
figure
hold on
axis equal
% Plot section
patch(zeros(size(Tnod))',Yn(Tnod)',Zn(Tnod)',ones(size(Tnod))','FaceColor',0.85*[1,1,1],'EdgeColor','none');
patch(zeros(size(Tbnd))',Yn(Tbnd)',Zn(Tbnd)',ones(size(Tbnd))','FaceColor','none','EdgeColor','k');
% Plot shear stress distribution
sig_ = sig;
sig_norm = (sig_-min(sig_))./(max(sig_)-min(sig_));
sig_ind = round(sig_norm*(size(cmap,1)-1))+1;
k = abs(sig_)*h2/5/max(abs(sig_))^2;
for i = 1:numel(Xp)
    quiver3(0,Xp(i),Yp(i),k(i)*sig(i),0,0,0,'filled','color',cmap(sig_ind(i),:));
end
colormap(cmap);
cb = colorbar;
caxis([min(sig_), max(sig_)])
title('Normal stress due to bending');
ylabel(cb,'\sigma (Pa)');
set(gca,'color','white','xcolor','none','ycolor','none','zcolor','none');
set(gcf,'color','white');
view(40,20);

% Plot tau Q
figure
hold on
axis equal
% Plot section
patch(Yn(Tnod)',Zn(Tnod)',ones(size(Tnod))','FaceColor',0.85*[1,1,1],'EdgeColor','none');
patch(Yn(Tbnd)',Zn(Tbnd)',ones(size(Tbnd))','FaceColor','none','EdgeColor','k');
% Plot applied load
quiver(xref,0,0,2*scale,0,'r','linewidth',3);
% Plot reference frame
quiver(0,0,scale,0,0,'k','linewidth',1.5);
text(1.1*scale,0,'x','HorizontalAlignment','center');
quiver(0,0,0,scale,0,'k','linewidth',1.5);
text(0,1.1*scale,'y','HorizontalAlignment','center');
% Plot shear stress distribution
tauQ_ = sqrt(tauQ(:,1).^2+tauQ(:,2).^2);
tauQ_norm = (tauQ_-min(tauQ_))./(max(tauQ_)-min(tauQ_));
tauQ_ind = round(tauQ_norm*(size(cmap,1)-1))+1;
k = tauQ_*h2/ny/max(tauQ_)^2;
for i = 1:numel(Xp)
    quiver(Xp(i),Yp(i),k(i)*tauQ(i,1),k(i)*tauQ(i,2),0,'filled','color',cmap(tauQ_ind(i),:));
end
colormap(cmap);
cb = colorbar;
caxis([min(tauQ_), max(tauQ_)])
title('Tangential stress due to shear load acting on shear centre')
ylabel(cb,'\tau (Pa)');
set(gca,'color','white','xcolor','none','ycolor','none');
set(gcf,'color','white');

% Plot tau T
figure
hold on
axis equal
% Plot section
patch(Yn(Tnod)',Zn(Tnod)',ones(size(Tnod))','FaceColor',0.85*[1,1,1],'EdgeColor','none');
patch(Yn(Tbnd)',Zn(Tbnd)',ones(size(Tbnd))','FaceColor','none','EdgeColor','k');
% Plot applied load
quiver(xref,0,0,2*scale,0,'r','linewidth',3);
% Plot reference frame
quiver(0,0,scale,0,0,'k','linewidth',1.5);
text(1.1*scale,0,'x','HorizontalAlignment','center');
quiver(0,0,0,scale,0,'k','linewidth',1.5);
text(0,1.1*scale,'y','HorizontalAlignment','center');
% Plot shear stress distribution
tauT_ = sqrt(tauT(:,1).^2+tauT(:,2).^2);
tauT_norm = (tauT_-min(tauT_))./(max(tauT_)-min(tauT_));
tauT_ind = round(tauT_norm*(size(cmap,1)-1))+1;
k = tauT_*h2/ny/max(tauT_)^2;
for i = 1:numel(Xp)
    quiver(Xp(i),Yp(i),k(i)*tauT(i,1),k(i)*tauT(i,2),0,'filled','color',cmap(tauT_ind(i),:));
end
colormap(cmap);
cb = colorbar;
caxis([min(tauT_), max(tauT_)])
title('Tangential stress due to torsion about the shear centre')
ylabel(cb,'\tau (Pa)');
set(gca,'color','white','xcolor','none','ycolor','none');
set(gcf,'color','white');

end