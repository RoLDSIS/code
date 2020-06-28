scale = 1;
triangle_size = 1*scale;
triangle_shift = 0.04*scale;
X = [1 0 1
     0 1 1
     1.2 1.2 1]*triangle_size + triangle_shift;
%
% Uncomment to draw the triangle defined by x1, x2 and x3.
%fill3(X(1,:),X(2,:),X(3,:),[1 1 1]*0,'faceAlpha',0.2);
clf
hold on
axis([0 2 0 2 0 2]*scale);
set(gca,'view',[120 30])
set(gca,'box','on')

M = mean(X,2);
%M = X(:,1);
X0 = X - repmat(M,1,3);
Cp = X0'*X0;
[U,S] = svd(Cp);
S(3,3) = 1e-10;
V = X0*U*S^(-1/2);
Z = V'*X0;
Z = Z(1:2,:);
plane_size = 4*scale;
vector_size = 1*scale;
mkr = 4;
% Uncomment to plot z orthonormal base
%arrow3(M',M'+V(:,1)'*vector_size)
%arrow3(M',M'+V(:,2)'*vector_size)
%plot3(M(1),M(2),M(3),'ko--','MarkerFaceColor','k','MarkerSize',mkr)
V1 = V(:,1);
V2 = V(:,2);
plane = M + [-V1+V2, V1+V2, V1-V2, -V1-V2]*plane_size; 
fill3(plane(1,:)',plane(2,:)',plane(3,:)',[1 1 1]*0,'faceAlpha',0.2);
text(X(1,1),X(2,1),X(3,1)+0.1,'$\mathbf{x_1}$','Interpreter','latex');
text(X(1,2),X(2,2),X(3,2)+0.1,'$\mathbf{x_2}$','Interpreter','latex');
text(X(1,3)+0.1,X(2,3)-0.1,X(3,3)-0.02,'$\mathbf{x_3}$','Interpreter','latex');
% Uncomment to plot z orthonormal base vector and point labels
%text(M(1)+0.1,M(2)-0.15,M(3),'$\mathbf{x_0}$','Interpreter','latex');
%text(M(1)+V(1,1)*vector_size,M(2)+V(2,1)*vector_size,M(3)+V(3,1)*vector_size+0.1,'$\mathbf{z_1}$','Interpreter','latex');
%text(M(1)+V(1,2)*vector_size,M(2)+V(2,2)*vector_size,M(3)+V(3,2)*vector_size-0.06,'$\mathbf{z_2}$','Interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1 = [1 9 10]';
A1 = [Z; ones(1,3)]'; 
d1 = inv(A1)*y1;
c1 = d1(3);
a1 = c1 - d1(1:2)'*V(:,1:2)'*M;
b1 = V(:,1:2)*d1(1:2);
d1_norm = d1(1:2)/sqrt(d1(1:2)'*d1(1:2));
b1_norm = V(:,1:2)*d1_norm;
arrow3(M',M'+b1_norm'*vector_size)
text(M(1)+b1_norm(1,1),M(2)+b1_norm(2,1),M(3)+b1_norm(3,1)+0.1,'$\mathbf{\Psi}$','Interpreter','latex');
arrow3(M'-b1_norm'*plane_size,M'+b1_norm'*plane_size,'--')
X1b = b1_norm' * X0(:,1) * b1_norm;
X2b = b1_norm' * X0(:,2) * b1_norm;
X3b = b1_norm' * X0(:,3) * b1_norm;
plot3([X(1,1) X1b(1)+M(1)],[X(2,1) X1b(2)+M(2)],[X(3,1) X1b(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,2) X2b(1)+M(1)],[X(2,2) X2b(2)+M(2)],[X(3,2) X2b(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,3) X3b(1)+M(1)],[X(2,3) X3b(2)+M(2)],[X(3,3) X3b(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
text(M(1)+X1b(1),M(2)+X1b(2),M(3)+X1b(3)+0.1,'$\mathbf{\psi_1}$','Interpreter','latex');
text(M(1)+X2b(1),M(2)+X2b(2)-0.1,M(3)+X2b(3)+0.1,'$\mathbf{\psi_2}$','Interpreter','latex');
text(M(1)+X3b(1),M(2)+X3b(2)+0.1,M(3)+X3b(3)+0.05,'$\mathbf{\psi_3}$','Interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y2 = [1 2 3]';
A2 = [Z; ones(1,3)]'; 
d2 = inv(A2)*y2;
c2 = d2(3);
a2 = c2 - d2(1:2)'*V(:,1:2)'*M;
b2 = V(:,1:2)*d2(1:2);
d2_norm = d2(1:2)/sqrt(d2(1:2)'*d2(1:2));
b2_norm = V(:,1:2)*d2_norm;
arrow3(M',M'+b2_norm'*vector_size)
text(M(1)+b2_norm(1,1),M(2)+b2_norm(2,1),M(3)+b2_norm(3,1)+0.1,'$\mathbf{\Phi}$','Interpreter','latex');
arrow3(M'-b2_norm'*plane_size,M'+b2_norm'*plane_size,'--')
X1b2 = b2_norm' * X0(:,1) * b2_norm;
X2b2 = b2_norm' * X0(:,2) * b2_norm;
X3b2 = b2_norm' * X0(:,3) * b2_norm;
mkr = 4;
plot3(M(1),M(2),M(3),'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,1) X1b2(1)+M(1)],[X(2,1) X1b2(2)+M(2)],[X(3,1) X1b2(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,2) X2b2(1)+M(1)],[X(2,2) X2b2(2)+M(2)],[X(3,2) X2b2(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,3) X3b2(1)+M(1)],[X(2,3) X3b2(2)+M(2)],[X(3,3) X3b2(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
text(M(1)+X1b2(1)+0.1,M(2)+X1b2(2)+0.1,M(3)+X1b2(3)+0.15,'$\mathbf{\phi_1}$','Interpreter','latex');
text(M(1)+X2b2(1)+0.1,M(2)+X2b2(2)-0.1+0.07,M(3)+X2b2(3)-0.1,'$\mathbf{\phi_2}$','Interpreter','latex');
text(M(1)+X3b2(1),M(2)+X3b2(2)+0.07,M(3)+X3b2(3)+0.05,'$\mathbf{\phi_3}$','Interpreter','latex');
set(gca,'xtick',[0 1])
set(gca,'ytick',[0 1])
set(gca,'ztick',[0 1])
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps = [900 300 500 500];
set(gcf, 'Position', ps);
ratio = ps(4)/ps(3);
paperWidth = 12;
paperHeight = paperWidth*ratio;

set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
shift = 1.5;
set(gcf, 'PaperPosition', [-shift -shift   paperWidth+1.75*shift paperHeight+1.75*shift]);

print('../../figures/fig_phy_psy_angle','-dpdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
