clear;
%usually, I paste the distribution from the clipboard to a spreadsheet, and
%export from there into a text file.  For example, c:\sedfit\cofsff0.txt


%some modification might be necessary if the distribution comes from the
%c(s,f) with 1 discrete component model, in order to take care of the first
%entry in the distribution being from the discrete species.

%one quick fix would be to manually change the s and ff0 value of the first
%species to that of the second (to make it part of the pattern of s and ff0
%values that the following script recognizes), and to set the distribution 
% value zero, and then not to plot the first data point

%I have tried to address this problem in lines 38-40, where the ff0 increment
%and the ff0 range is determined

%this is the revision from 12/08/2007, where a bug was fixed when using
%non-equidistant s-values.


distfile = 'C:\sedfit\cofsff0.txt';




fin = fopen(distfile);
distribution = fscanf(fin,'%f',[6 inf]);
% distribution = distsribution';


ff0vec = distribution(3,:);
svec = distribution(1,:);
Mvec = distribution(2,:); Mvec = Mvec/1000;
Dvec = distribution(4,:); Dvec = Dvec*1e7;
Rvec = distribution(5,:); 
cvec = distribution(6,:);

fclose(fin);

ff0min = min(ff0vec(2:length(ff0vec)));
ff0max = max(ff0vec(2:length(ff0vec)));
ff0inc = max(diff(ff0vec(2:length(ff0vec))));

smin = min(svec);
smax = max(svec);

% sinc = max(diff(svec));

Mmin = min(Mvec); Mmax = max(Mvec);
Dmin = min(Dvec); Dmax = max(Dvec);
Rmin = min(Rvec); Rmax = max(Rvec);
cmin = min(cvec); cmax = max(cvec);

nff0 = 1+round((ff0max-ff0min)/ff0inc);
ns = length(distribution)/nff0;

% svals = linspace(smin,smax,ns);
% ff0vals = linspace(ff0min,ff0max,nff0);
% sff0mat = zeros(ns,nff0);
% for s=1:ns,
%     for f=1:nff0,
%         sff0mat(s,f)=cvec(s+(f-1)*ns);
%     end;
% end;
% 
% 
% [S,F] = meshgrid(svals,ff0vals);
% 
% figure(1); clf;
% h = surfc(S,F,sff0mat');
% % contour(S,F,sff0mat,100);
% ylabel('f/f0');
% xlabel('s-value (S)');
% zlabel('c(s,f/f0)');
% zoffset = -0.3*max(max(sff0mat));
% for i=2:length(h),
%    newz = get(h(i),'Zdata') + zoffset;
%    set(h(i),'Zdata',newz);
% end;
% axis([ smin smax ff0min ff0max zoffset 1.1*max(cvec)]);
% view(30,90);
% % shading interp;
% alpha(0.1);
% rotate3d;


% 
% figure(2); clf;
% h = meshc(S,F,sff0mat');
% ylabel('f/f0');
% xlabel('s-value (S)');
% zlabel('c(s,f/f0)');
% zoffset = -0.4*max(max(sff0mat));
% for i=2:length(h),
%    newz = get(h(i),'Zdata') + zoffset;
%    set(h(i),'Zdata',newz);
% end;
% axis([ smin smax ff0min ff0max zoffset 1.1*max(cvec)]);
% view(0,90);
% % shading interp;
% alpha(0.1);
% rotate3d;
% 
% 
% figure(21); clf;
% h = contour(S,F,sff0mat');
% ylabel('f/f0');
% xlabel('s-value (S)');
% colorbar;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cmat = zeros(ns,nff0);
Smat = zeros(ns,nff0);
Mmat = zeros(ns,nff0);
Dmat = zeros(ns,nff0);
Rmat = zeros(ns,nff0);
Fmat = zeros(ns,nff0);
for s=1:ns,
    for f=1:nff0,
        Cmat(s,f)=cvec(s+(f-1)*ns);
        Smat(s,f) = svec(s+(f-1)*ns);
        Mmat(s,f) = Mvec(s+(f-1)*ns);
        Dmat(s,f) = Dvec(s+(f-1)*ns);
        Rmat(s,f) = Rvec(s+(f-1)*ns);
        Fmat(s,f) = ff0vec(s+(f-1)*ns);
    end;
end;

% %average f/f0
% for s=1:ns,
%     sfcurve(1,s) = svals(s);
%     addc = 0;
%     addf = 0;
%     for f=1:nff0,
%         addc = addc + cvec(s+(f-1)*ns);
%         addf = addf + cvec(s+(f-1)*ns)*ff0vec(s+(f-1)*ns);
%     end;
%     sfcurve(2,s) = addf/addc;
%     sfcurve(3,s) = addc;
% end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot s-f


figure(2); clf;
% h = meshc(Smat,Mmat,Cmat);
h = contour(Smat,Fmat,Cmat);
hold on;
for f=1:nff0,
    for s=1:ns,
        tmps(s) = svec(s+(f-1)*ns);
        tmpf(s) = ff0vec(s+(f-1)*ns);
    end;
    plot(tmps,tmpf,'b:');
end;
hold off;
xlabel('s-value (S)');
ylabel('ff0');
zlabel('c(s,ff0)');
zoffset =0;
axis([ smin smax ff0min ff0max zoffset 1.1*max(cvec)]);
view(0,90);
rotate3d;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot s-M


figure(3); clf;
% h = meshc(Smat,Mmat,Cmat);
h = contour(Smat,Mmat,Cmat);
hold on;
for f=1:nff0,
    for s=1:ns,
        tmps(s) = svec(s+(f-1)*ns);
        tmpM(s) = Mvec(s+(f-1)*ns);
    end;
    plot(tmps,tmpM,'b:');
end;
hold off;
xlabel('s-value (S)');
ylabel('Mw (kDa)');
zlabel('c(s,M)');
zoffset =0;
axis([ smin smax Mmin Mmax zoffset 1.1*max(cvec)]);
view(0,90);
rotate3d;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot s-D


figure(4); clf;
h = contour(Smat,Dmat,Cmat);
hold on;
for f=1:nff0,
    for s=1:ns,
        tmps(s) = svec(s+(f-1)*ns);
        tmpD(s) = Dvec(s+(f-1)*ns);
    end;
    plot(tmps,tmpD,'b:');
end;
hold off;
xlabel('s-value (S)');
ylabel('D (10^-7)');
zlabel('c(s,M)');
zoffset =0;
axis([ smin smax 0 Dmax zoffset 1.1*max(cvec)]);
view(0,90);
rotate3d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot s-R


figure(5); clf;
h = contour(Smat,Rmat,Cmat);
hold on;
for f=1:nff0,
    for s=1:ns,
        tmps(s) = svec(s+(f-1)*ns);
        tmpR(s) = Rvec(s+(f-1)*ns);
    end;
    plot(tmps,tmpR,'b:');
end;
hold off;
xlabel('s-value (S)');
ylabel('RStokes (nm)');
zlabel('c(s,M)');
zoffset =0;
axis([ smin smax 0 Rmax zoffset 1.1*max(cvec)]);
view(0,90);
rotate3d;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot M-ff0


figure(6); clf;
h = contour(Mmat,Fmat,Cmat);
hold on;
% for f=1:nff0,
%     for s=1:ns,
%         tmpM(s) = Mvec(s+(f-1)*ns);
%         tmpF(s) = ff0vec(s+(f-1)*ns);
%     end;
%     plot(tmpM,tmpF,'b:');
% end;
for s=1:ns/10,
    for f=1:nff0,
        tmpM2(f) = Mvec(1+10*s+(f-1)*ns);
        tmpF2(f) = ff0vec(1+10*s+(f-1)*ns);
    end;
    plot(tmpM2,tmpF2,'k:');
end;
hold off;
xlabel('molar mass (kDa)');
ylabel('frictional ratio');
zlabel('c(s,M)');
zoffset =0;
axis([ Mmin Mmax 0.9*ff0min 1.1*ff0max zoffset 1.1*max(cvec)]);
view(0,90);
rotate3d;

mwlineyvals(1) = 1; mwlineyvals(2) = 2.2;
mwx1vals(1) = 6; mwx1vals(2) = 6;
mwx2vals(1) = 30; mwx2vals(2) = 30;
hold on;
plot(mwx1vals,mwlineyvals,'r-');
plot(mwx2vals,mwlineyvals,'r-');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot M-R


figure(7); clf;
h = contour(Mmat,Rmat,Cmat);
hold on;
for f=1:nff0,
    for s=1:ns,
        tmpM(s) = Mvec(s+(f-1)*ns);
        tmpR(s) = Rvec(s+(f-1)*ns);
    end;
    plot(tmpM,tmpR,'b:');
end;
for s=1:ns/10,
    for f=1:nff0,
        tmpM2(f) = Mvec(1+10*s+(f-1)*ns);
        tmpF2(f) = Rvec(1+10*s+(f-1)*ns);
    end;
    plot(tmpM2,tmpF2,'k:');
end;
hold off;
xlabel('molar mass (kDa)');
ylabel('Stokes radius (nm)');
zlabel('c(s,M)');
zoffset =0;
axis([ Mmin Mmax 0 Rmax zoffset 1.1*max(cvec)]);
view(0,90);
rotate3d;

