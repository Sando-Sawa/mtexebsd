%% Grain reconstruction plot created by Sando Sawa (Tohoku University)
% version 1.0.0

%% #1 map direction configuration (necessary)

xdi = 'east';
ydi = 'intoPlane';
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);

figure;
plot(ebsd)
%% #2 Grain boundary setting (all pixels)　（下処理）

% MAD<1.0
ebsdf = ebsd(ebsd.mad<1.0);
ebsd_corrected = ebsdf;

% reconstruct grains with theshold angle 10 degree
% only indexed mineral（埋まる）
%[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected('indexed'),'theshold',10*degree);
% all mineral （埋まらない）
[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected,'angle',10*degree);

% remove small grains(pixel)
ebsd_corrected(grains(grains.grainSize <= 2)) = [];

% reconstruct grains with theshold angle 10 degree
% only indexed mineral（埋まる）
%[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected('indexed'),'theshold',10*degree);
% all mineral （埋まらない）
[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected,'angle',10*degree);

% smooth the grains to avoid the stair casing effect
grains = smooth(grains,5);

%% #3 Plot all (grains)　（結晶相カラー）

% Color change
grains('Forsterite').CS.color=str2rgb('deepskyblue');
%grains('Diopside').CS.color=str2rgb('salmon');
%grains('Antigorite').CS.color=str2rgb('limegreen');
%grains('Talc').CS.color=str2rgb('salmon');
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);
%結晶相カラー
figure;
plot(grains)
hold on
plot(grains.boundary,'lineWidth',0.5)
hold off

%% #4 Pole figure for scatter data
%X-Z direction
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'X','Z'},...
'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);

%which plane?
h = Miller({1,0,0},{0,1,0},{0,0,1},grains('Forsterite').CS);
figure;
plotPDF(grains('Forsterite').meanOrientation,h,'figSize','small','MarkerSize',1,'MarkerColor','black','FontSize',36)

%% #5 Pole figure with contour
%X-Z direction
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'X','Z'},...
'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);

%which plane?
h = Miller({1,0,0},{0,1,0},{0,0,1},grains('Forsterite').CS);
figure;
%mtexColorbar % remove colorbars
%CLim(figure,'equal');
%mtexColorbar % add a single colorbar
%plotPDF(grains('Forsterite').meanOrientation,h,'contourf','figSize','small','MarkerSize',1,'MarkerColor','black','FontSize',36)
plotPDF(grains('Forsterite').meanOrientation,h,'contourf','halfwidth',10*degree,'figSize','small','MarkerSize',1,'MarkerColor','black','FontSize',36)
CLim(gcm,[0 30]);
mtexColorbar % add a single colorbar
%% #6 Pole figure for contour without line
%X-Z direction
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'X','Z'},...
'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);

% compute the ODF with the kernel psi
%odf = calcDensity(ebsd_corrected('Forsterite').orientations,'halfwidth',10*degree)
odf = calcDensity(grains('Forsterite').meanOrientation,'halfwidth',10*degree)

h = [Miller(1,0,0,odf.CS),Miller(0,1,0,odf.CS),Miller(0,0,1,odf.CS)];

figure;
plotPDF(odf,h,'antipodal','silent')

%% #7 IPF

setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);
%IPF
ipfKey = ipfColorKey(grains('Forsterite'));
ipfKey.inversePoleFigureDirection = vector3d.Z;
% this is the colored fundamental sector
figure;
plot(ipfKey)
csFo = grains('Forsterite').CS;
colorKey = BungeColorKey(csFo);
% this computes the colors for each orientation specified as input
colors = ipfKey.orientation2color(grains('Forsterite').meanOrientation);

% this plots the grains colorized according to the RGB values stored in colors
figure;
plot(grains('Forsterite'),colors)

% plot the grain boundaries on top of the ipf map
hold on
plot(grains.boundary,'lineWidth',0.5)
%plot(grains('notIndexed'),'FaceColor','white')
hold off

%% #8 misorientation (mtex version >= 5.5.2)

% MAD<1.0
ebsdg = ebsd(ebsd.mad<1.0);

%[grains,ebsd.grainId] = calcGrains(ebsd('indexed'));
[grainsg,ebsdg.grainId,ebsdg.mis2mean] = calcGrains(ebsdg);
% remove one-three pixel grains
ebsdg(grainsg(grainsg.grainSize <= 3)) = [];
%[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',[1*degree, 10*degree]);
[grainsg,ebsdg.grainId,ebsdg.mis2mean] = calcGrains(ebsdg,'threshold',[1*degree, 10*degree]);

% smooth grain boundaries
grainsg = smooth(grainsg,5);

% choose forsterite
ebsdg=ebsdg('Forsterite')
grainsh=grainsg('Forsterite')
% denoise the orientations
F = halfQuadraticFilter;
ebsdg = smooth(ebsdg,F,grainsg,'fill');

%plot(ebsd('indexed'),ebsd('indexed').orientations)
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);
figure;
plot(ebsdg,ebsdg.orientations)
hold on
plot(grainsg.boundary,'lineWidth',2)
hold off

% compute the grain reference orientation deviation
grod = ebsdg.calcGROD(grainsg);

% plot the misorientation angle of the GROD
figure;
%plot(ebsdg,grod.angle./degree,'micronbar','off','colorrange',[0 4])
plot(ebsdg,grod.angle./degree,'micronbar','off')
CLim(gcm,[0 4]);
mtexColorbar('title','misorientation angle to meanorientation in degree')
mtexColorMap LaboTeX

% overlay grain and subgrain boundaries
hold on
plot(grainsg.boundary,'lineWidth',1.5)
plot(grainsh.innerBoundary,'edgeAlpha',grainsh.innerBoundary.misorientation.angle / (5*degree))
hold off

%% #9 misorientation histgram at one point

pointebsd=ebsdg(grainsh(42,42));

% get the misorientations to mean
mori = pointebsd.mis2mean

% plot a histogram
figure;
plotAngleDistribution(mori)
xlabel('Misorientation angles in degree')
set(gca, 'Yscale', 'log');

%% #10 dislocation density

% reconstruct grains
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',5*degree);

% remove small grains
ebsd(grains(grains.grainSize<=5)) = [];

% redo grain reconstruction
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',2.5*degree);

% smooth grain boundaries
grains = smooth(grains,5);

hold on
plot(grains.boundary,'linewidth',2)
hold off

%%

ebsd=ebsd('Forsterite');

% a key the colorizes according to misorientation angle and axis
ipfKey = axisAngleColorKey(ebsd);

% set the grain mean orientations as reference orinetations
ipfKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;

% plot the data
plot(ebsd,ipfKey.orientation2color(ebsd('indexed').orientations),'micronBar','off','figSize','medium')

hold on
plot(grains.boundary,'linewidth',2)
hold off

% denoise orientation data
F = halfQuadraticFilter;

ebsd = smooth(ebsd('indexed'),F,'fill',grains);

% plot the denoised data
ipfKey.oriRef = grains(ebsd('indexed').grainId).meanOrientation;
plot(ebsd('indexed'),ipfKey.orientation2color(ebsd('indexed').orientations),'micronBar','off','figSize','medium')

hold on
plot(grains.boundary,'linewidth',2)
hold off

%%

% consider only the Forsterite phase
ebsd_curv = ebsd('indexed').gridify;

% compute the curvature tensor
kappa = ebsd_curv.curvature

% one can index the curvature tensors in the same way as the EBSD data.
% E.g. the curvature in pixel (2,3) is
kappa(2,3)

% According to Kroener the curvature tensor is directly related to the dislocation density tensor.
alpha = kappa.dislocationDensity

figure;
newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsd,kappa{i,j},'micronBar','off')
    hold on; plot(grains.boundary,'linewidth',2); hold off

  end
end

% unify the color rage  - you may also use setColoRange equal
setColorRange([-0.005,0.005])
drawNow(gcm,'figSize','large')

%%

figure;
newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsd,alpha{i,j},'micronBar','off')
    hold on; plot(grains.boundary,'linewidth',2); hold off

  end
end

mtexColorMap('hot')
mtexColorbar
% unify the color rage  - you may also use setColoRange equal
%setColorRange([-0.005,0.005])
drawNow(gcm,'figSize','large')


%% test dislocation

cs=crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]);
dS = [dislocationSystem(Miller(1,0,0,cs,'uvw'),Miller(0,1,0,cs,'uvw')),...
    dislocationSystem(Miller(1,0,0,cs,'uvw'),Miller(0,0,1,cs,'uvw')),...
    dislocationSystem(Miller(0,0,1,cs,'uvw'),Miller(1,0,0,cs,'uvw')),...
    dislocationSystem(Miller(0,0,1,cs,'uvw'),Miller(0,1,0,cs,'uvw')),...
    dislocationSystem(Miller(1,0,0,cs,'uvw'),Miller(1,0,0,cs,'uvw')),...
    dislocationSystem(Miller(0,0,1,cs,'uvw'),Miller(0,0,1,cs,'uvw'))]

%%

%dS = dislocationSystem.fcc(ebsd.CS)

% size of the unit cell
%a = norm(ebsd.CS.aAxis);

% in bcc and fcc the norm of the burgers vector is sqrt(3)/2 * a
%[norm(dS(1).b), 25/12 * a, norm(dS(end).b)]

%
%R =
%r_0 =
%U = norm(dS.b).^2

%nu = 0.3;


% energy of the edge dislocations
%dS(dS.isEdge).u = 1;

% energy of the screw dislocations
%dS(dS.isScrew).u = 1 - 0.3;

% Question to verybody: what is the best way to set the enegry? I found
% different formulae
%
% E = 1 - poisson ratio
% E = c * G * |b|^2,  - G - Schubmodul / Shear Modulus Energy per (unit length)^2

dSRot = ebsd.orientations * dS

[rho,factor] = fitDislocationSystems(kappa,dSRot);

% the restored dislocation density tensors
alpha = sum(dSRot.tensor .* rho,2);

% we have to set the unit manualy since it is not stored in rho
alpha.opt.unit = '1/um';

% the restored dislocation density tensor for pixel 2
%alpha(2)

% the dislocation density dervied from the curvature in pixel 2
%kappa(2).dislocationDensity

kappa = alpha.curvature

figure;
newMtexFigure('nrows',3,'ncols',3);

% cycle through all components of the tensor
for i = 1:3
  for j = 1:3

    nextAxis(i,j)
    plot(ebsd,kappa{i,j},'micronBar','off')
    hold on; plot(grains.boundary,'linewidth',2); hold off

  end
end

setColorRange([-0.005,0.005])
drawNow(gcm,'figSize','large');

%close all
figure;
plot(ebsd,factor*sum(abs(rho .* dSRot.u),2),'micronbar','off')
mtexColorMap('hot')
mtexColorbar

set(gca,'ColorScale','log'); % this works only starting with Matlab 2018a
set(gca,'CLim',[1e11 5e14]);

hold on
plot(grains.boundary,'linewidth',2)
hold off