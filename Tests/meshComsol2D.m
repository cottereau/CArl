function out = meshComsol2D
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
model.geom('geom1').feature('r2').set('pos', {'.25' '.25'});
model.geom('geom1').feature('r2').set('size', {'.5' '.5'});
model.geom('geom1').run;

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
model.mesh('mesh1').run;

out = model;
