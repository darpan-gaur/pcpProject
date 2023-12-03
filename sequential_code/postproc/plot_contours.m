basedir = '../serial_code/';
nxarr = [20 40 80]; nyarr = [10 20 40]; iter = 0;

for it = 1:length(nxarr)
  nx = nxarr(it); ny = nyarr(it);
  fname = strcat(basedir,'T_xy_', num2str(nx,'%03d'), '_', num2str(ny,'%03d'), '_', num2str(iter, '%04d'), '.dat')
  rawdat = dlmread(fname);
  x = rawdat(1:ny:nx*ny,1); y = rawdat(1:ny,2); 
  T   = reshape(rawdat(:,3), [ny,nx])'; 
  Tex = reshape(rawdat(:,4), [ny,nx])';

  figure(1), clf
  contourf(x, y, T', 'LineColor', 'none')
  xlabel('x'), ylabel('y'), pbaspect([2 1 1])
  title(strcat('Numerical :', num2str(nx,'%03d'), 'x', num2str(ny,'%03d')))
  axis([0.0 1.0 0.0 0.5]), caxis([0 0.25])
  colorbar,   colormap('jet')
  set(gca,'fontsize',14)
  figfname = strcat(basedir,'contours_numer_T_xy_', num2str(nx,'%03d'), '_', num2str(ny,'%03d'), '_', num2str(iter, '%04d'), '.png');
  screen2jpeg(figfname);

  figure(1), clf
  contourf(x, y, Tex', 'LineColor', 'none')
  xlabel('x'), ylabel('y'), pbaspect([2 1 1])
  title(strcat('Numerical :', num2str(nx,'%03d'), 'x', num2str(ny,'%03d')))
  axis([0.0 1.0 0.0 0.5]), caxis([0 0.25])
  colorbar,   colormap('jet')
  set(gca,'fontsize',14)
  figfname = strcat(basedir,'contours_exact_T_xy_', num2str(nx,'%03d'), '_', num2str(ny,'%03d'), '_', num2str(iter, '%04d'), '.png');
  screen2jpeg(figfname);
end
