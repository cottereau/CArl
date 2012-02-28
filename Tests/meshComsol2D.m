function out = model
%
% meshComsol2D.m
%
% Model exported on Feb 21 2012, 15:34 by COMSOL 4.2.1.166.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/Users/cottereau/matlab/CArl/Tests');

model.modelNode.create('mod1');

model.geom.create('geom1', 2);
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature.create('r2', 'Rectangle');
model.geom('geom1').feature('r1').set('pos', {'0' '0'});
model.geom('geom1').feature('r1').set('size', {'1' '1'});
model.geom('geom1').feature('r2').set('pos', {'.25' '.25'});
model.geom('geom1').feature('r2').set('size', {'.5' '.5'});
model.geom('geom1').run;

model.physics.create('hteq', 'HeatEquation', 'geom1');

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');

% model.view('view1').axis.set('xmin', '-0.05000000074505806');
% model.view('view1').axis.set('xmax', '1.0499999523162842');
% model.view('view1').axis.set('ymin', '-0.3482394218444824');
% model.view('view1').axis.set('ymax', '1.3482393026351929');

model.physics('hteq').feature('hteq1').set('f', '0');
model.physics('hteq').feature('hteq1').set('da', '0');

model.mesh('mesh1').run;

model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');

out = model;
