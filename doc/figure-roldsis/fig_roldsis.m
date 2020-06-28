%
% Script to plot fig_roldsis.pdf which illustrates the RoLDSIS (Regression on Low
% Dimensional Spanned Input Space
%
scale = 1;
triangle_size = 1*scale;
triangle_shift = 0.04*scale;
% The columns of X contain the position of the 3 3D points that span the 2D
% subspace used to constrain the solutions of EEG High Dimensional Low Sample
% Size (HDLSS) problems.
X = [1 0 1
     0 1 1
     1.2 1.2 1]*triangle_size + triangle_shift;
% Plot the triangle defined by the columns of X.
h=fill3(X(1,:),X(2,:),X(3,:),[1 1 1]*0,'faceAlpha',0.2);
axis([0 2 0 2 0 2]*scale);
set(gca,'view',[120 30])
set(gca,'box','on')

% Take the mean value of X as the reference point in the subspace spanned by X.
% Any linear combination of the columns of X could be used as reference.
M = mean(X,2);
%M = X(:,1);
% The columns of X0 contain the vectors from the reference point M to the points
% defined by the columns of X
X0 = X - repmat(M,1,3);
% Use PCA to calculate an orthonormal base for the subspace spanned by X.
% Note that PCA is computed from SVD of X0'*X0 and not X0*X0' which, in HDLSS
% problems, has a much higher dimension.
Cp = X0'*X0;
[U,S] = svd(Cp);
S(3,3) = 1e-10;
V = X0*U*S^(-1/2);
% The first 2 columns of V contain the eigenvectors corresponding to the
% 2 non zero eigenvalues of X0*X0'. Z contains X0 represented in the orthonormal
% base defined by V.
Z = V'*X0;
Z = Z(1:2,:);
hold on
vector_size = 1*scale;
% Plot the orthonormal base defined by V.
arrow3(M',M'+V(:,1)'*vector_size)
arrow3(M',M'+V(:,2)'*vector_size)
plane_size = 4*scale;
V1 = V(:,1);
V2 = V(:,2);
% Draw the plane spanned by X as a partially transparent rectangle.
plane = M + [-V1+V2, V1+V2, V1-V2, -V1-V2]*plane_size; 
h=fill3(plane(1,:)',plane(2,:)',plane(3,:)',[1 1 1]*0,'faceAlpha',0.2);
text(X(1,1),X(2,1),X(3,1)+0.1,'$\mathbf{x_1}$','Interpreter','latex');
text(X(1,2),X(2,2),X(3,2)+0.1,'$\mathbf{x_2}$','Interpreter','latex');
text(X(1,3),X(2,3)+0.06,X(3,3)-0.02,'$\mathbf{x_3}$','Interpreter','latex');
text(M(1)+0.15,M(2)-0.06,M(3),'$\mathbf{m}$','Interpreter','latex');
text(M(1)+V(1,1)*vector_size,M(2)+V(2,1)*vector_size,M(3)+V(3,1)*vector_size+0.1,'$\mathbf{z_1}$','Interpreter','latex');
text(M(1)+V(1,2)*vector_size,M(2)+V(2,2)*vector_size,M(3)+V(3,2)*vector_size-0.06,'$\mathbf{z_2}$','Interpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the problem y = a + bx constrained to the plane spanned by X: c + dz 
y1 = [1 7 9]';
A1 = [Z; ones(1,3)]'; 
d1 = inv(A1)*y1;
c1 = d1(3);
a1 = c1 - d1(1:2)'*V(:,1:2)'*M;
b1 = V(:,1:2)*d1(1:2);
d1_norm = d1(1:2)/sqrt(d1(1:2)'*d1(1:2));
b1_norm = V(:,1:2)*d1_norm;
% Plot the normalized vector b and, in dashed lines, the direction defined by it.
arrow3(M',M'+b1_norm'*vector_size)
text(M(1)+b1_norm(1,1),M(2)+b1_norm(2,1),M(3)+b1_norm(3,1)+0.1,'$\mathbf{b}$','Interpreter','latex');
arrow3(M'-b1_norm'*plane_size,M'+b1_norm'*plane_size,'--')
X1b = b1_norm' * X0(:,1) * b1_norm;
X2b = b1_norm' * X0(:,2) * b1_norm;
X3b = b1_norm' * X0(:,3) * b1_norm;
% Plot the points contained in X and its projections on the direction defined by b.
mkr = 4;
plot3(M(1),M(2),M(3),'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,1) X1b(1)+M(1)],[X(2,1) X1b(2)+M(2)],[X(3,1) X1b(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,2) X2b(1)+M(1)],[X(2,2) X2b(2)+M(2)],[X(3,2) X2b(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
plot3([X(1,3) X3b(1)+M(1)],[X(2,3) X3b(2)+M(2)],[X(3,3) X3b(3)+M(3)],'ko--','MarkerFaceColor','k','MarkerSize',mkr)
text(M(1)+X1b(1),M(2)+X1b(2),M(3)+X1b(3)+0.1,'$\mathbf{\tilde{x}_1}$','Interpreter','latex');
text(M(1)+X2b(1),M(2)+X2b(2)+0.06,M(3)+0.01+X2b(3),'$\mathbf{\tilde{x}_2}$','Interpreter','latex');
text(M(1)+X3b(1),M(2)+X3b(2)+0.07,M(3)+X3b(3)+0.05,'$\mathbf{\tilde{x}_3}$','Interpreter','latex');
set(gca,'xtick',[0 1])
set(gca,'ytick',[0 1])
set(gca,'ztick',[0 1])
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format the figure to fill the pdf file
ps = [900 300 500 500];
set(gcf, 'Position', ps);
ratio = ps(4)/ps(3);
paperWidth = 12;
paperHeight = paperWidth*ratio;

set(gcf, 'paperunits', 'centimeters');
set(gcf, 'papersize', [paperWidth paperHeight]);
shift = 1.5;
set(gcf, 'PaperPosition', [-shift -shift   paperWidth+1.75*shift paperHeight+1.75*shift]);

print('../../figures/fig_roldsis','-dpdf')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

