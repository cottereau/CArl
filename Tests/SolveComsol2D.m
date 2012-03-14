function out = SolveComsol2D( model )
%
% SolveComsol.m
%
% Model exported on Mar 13 2012, 17:37 by COMSOL 4.2.1.110.
model.physics.create('ht', 'HeatTransfer', 'geom1');

model.physics('ht').feature('solid1').set('k_mat', 'userdef');
model.physics('ht').feature('solid1').set('k', {'1'; '0'; '0'; '0'; '1'; '0'; '0'; '0'; '1'});

model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
%model.sol('sol1').feature('v1').set('scalemethod','none');

out = model;
