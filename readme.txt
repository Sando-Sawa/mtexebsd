Release Note of ebsd_filled_plot.m
製作者：Sando Sawa
バグ報告はsando.cuso4@mail.comまで

v1.0.0
初期リリース版


使い方

まずはじめに
1. EBSDデータをインポート
2. #1（のセクション）を実行し、図の向きを調整する。
CHANNEL 5で得られたマップと見比べながら、figureタブのMTEX→x axis directionとz axis directionを調整。
一致するaxis directionがわかったら、スクリプトのxdi = 'east';とydi = 'intoPlane';を書き換える。
例：x axis direction: north, z axis direction: Out of planeのときにCHANNEL 5で得られたマップと一致する場合
→xdi = 'north';
　ydi = 'OutOfPlane';

結晶相のマップを描きたいとき
#2→#3

極図を描きたいとき
#2（まだやっていなければ）→#4（点）
　　　　　　　　　　　　　→#5（線付きコンター）
                          →#6（線なしコンター）

IPFを描きたいとき
#2（まだやっていなければ）→#7

misorientationのマップを描きたいとき
#8

構造地質学分野なので極図はX-Zプロットになっています。
X-Yに直したいときは
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'X','Z'},...
'BackgroundColor','w','tag','axesLabels',varargin{:});
を
pfAnnotations = @(varargin) text([vector3d.X,vector3d.Y],{'X','Y'},...
'BackgroundColor','w','tag','axesLabels',varargin{:});
に直してください

極図はone point per grain（1粒子につき1点）になっています。
1点ずつに変えたい場合は
grains('Forsterite').meanOrientation
を
ebsd_corrected('Forsterite').orientations
に変えてください

