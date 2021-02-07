%% Grain reconstruction plot created by Sando
%% #1 map direction configuration (necessary)

xdi = 'east';
ydi = 'intoPlane';
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);

figure;
plot(ebsd)

%% rectangle
plot(ebsd);
region = [120 60 100 60];
rectangle('position',region,'edgecolor','k','linewidth',2)

%% cut
figure;
ebsd = ebsd(inpolygon(ebsd,region))
plot(ebsd)
%% #2 Grain boundary setting (all pixels)　（下処理）

% MAD<1.0
ebsdf = ebsd(ebsd.mad<1.0);
ebsd_corrected = ebsdf;

% reconstruct grains with theshold angle 10 degree
% only indexed mineral（埋まる）
[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected('indexed'),'theshold',15*degree);
% all mineral （埋まらない）
%[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected,'angle',15*degree);

% remove small grains(pixel)
ebsd_corrected(grains(grains.grainSize <= 2)) = [];

% reconstruct grains with theshold angle 10 degree
% only indexed mineral（埋まる）
[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected('indexed'),'theshold',15*degree);
% all mineral （埋まらない）
%[grains,ebsd_corrected.grainId] = calcGrains(ebsd_corrected,'angle',15*degree);

% smooth the grains to avoid the stair casing effect
grains = smooth(grains,5);

%% Plot all (grains)　（結晶相カラー）

% Color change
grains('germanate olivine').CS.color=str2rgb('blue');
%grains('Diopside').CS.color=str2rgb('salmon');
%grains('Antigorite').CS.color=str2rgb('limegreen');
%grains('Talc').CS.color=str2rgb('gold');

setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);
%結晶相カラー
figure;
plot(grains)
hold on
plot(grains.boundary,'lineWidth',0.5)
hold off

% grain diameter
d = diameter(grains('indexed'));
ave = mean(d)


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
plotPDF(grains('germanate olivine').meanOrientation,h,'figSize','small','MarkerSize',1,'MarkerColor','black','FontSize',36)

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
plotPDF(grains('germanate olivine').meanOrientation,h,'contourf','halfwidth',10*degree,'figSize','small','MarkerSize',1,'MarkerColor','black','FontSize',36)
%% #6 Pole figure for contour without line
%X-Z direction
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'X','Z'},...
'BackgroundColor','w','tag','axesLabels',varargin{:});
setMTEXpref('pfAnnotations',pfAnnotations);
setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);

% compute the ODF with the kernel psi
%odf = calcDensity(ebsd_corrected('Forsterite').orientations,'halfwidth',10*degree)
odf = calcDensity(grains('germanate olivine').meanOrientation,'halfwidth',10*degree)

h = [Miller(1,0,0,odf.CS),Miller(0,1,0,odf.CS),Miller(0,0,1,odf.CS)];

figure;
plotPDF(odf,h,'antipodal','silent')

%% #7 IPF

setMTEXpref('xAxisDirection',xdi);
setMTEXpref('zAxisDirection',ydi);
%IPF
ipfKey = ipfColorKey(grains('germanate olivine'));
ipfKey.inversePoleFigureDirection = vector3d.Z;
% this is the colored fundamental sector
figure;
plot(ipfKey)
csFo = grains('germanate olivine').CS;
colorKey = BungeColorKey(csFo);
% this computes the colors for each orientation specified as input
colors = ipfKey.orientation2color(grains('germanate olivine').meanOrientation);

% this plots the grains colorized according to the RGB values stored in colors
figure;
plot(grains('germanate olivine'),colors)

% plot the grain boundaries on top of the ipf map
hold on
plot(grains.boundary,'lineWidth',0.5)
%plot(grains('notIndexed'),'FaceColor','white')
hold off

%% mesh

G=gmshGeo(grains);
mesh(G,'hip2.msh')