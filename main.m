clear; close all; clc;

format shortEng

addpath("classes\")
addpath("FilesForAgglo\")

load('trainDat')
load('testDat')


% initilize dmd model
mod = dmdmodel();

% train dmd model
mod = mod.train(trainDat,4,0.95,0.95,false); 

% model predictions on training data
[mod, trainDatPred] = mod.test(trainDat);

% model predictions on test data
[mod, testDatPred] = mod.test(testDat);

%% Simulation
% gather parameter for simulink simulation
pmod = mod.p;
pagg = trainDat.p;

% generate input for simulink
Tmax = 3*3600;
inputFunction = generateInputFunction(mod.uop,mod.ts,Tmax/mod.ts);
temperature.time = 0:1:Tmax;
temperature.signals.values = inputFunction(temperature.time);
% temperature.signals.dimensions = [1];
r0 = 0.6;
referenceFunction = generateReferenceFunction(r0,mod.ts,Tmax/mod.ts);
reference.time = 0:1:Tmax;
reference.signals.values = referenceFunction(reference.time);
writematrix([reference.time'/3600 reference.signals.values],['FilesForAgglo\results\data\d32ref','control','.csv'])

disturbanceFunction = generateDisturbanceFunction(0.2,mod.ts,Tmax/mod.ts);
disturbance.time = 0:1:Tmax;
disturbance.signals.values = disturbanceFunction(disturbance.time);
writematrix([disturbance.time'/3600 disturbance.signals.values],['FilesForAgglo\results\data\dfeed','controlDist','.csv'])

pmod.U = mod.U;
pmod.N_x = mod.N_x;
pmod.N_x_transformed = mod.N_x_transformed;
pmod.N_del = mod.N_del;
pmod.xop = mod.xop;
pmod.uop = mod.uop;
pmod.stateScale = mod.stateScale;
pmod.inputScale = mod.inputScale;
pmod.ext = mod.ext;

% MPC object
mpcobj = mpc(mod.stateSpace);
mpcobj.ManipulatedVariables.Min = (353 - pmod.uop)/pmod.inputScale;
mpcobj.ManipulatedVariables.Max = (373 - pmod.uop)/pmod.inputScale;
mpcobj.PredictionHorizon = 30;
mpcobj.ControlHorizon = 15;

out = sim("ClosedLoopSim");

%% after simulation
% gather datasets
label = 'controlDist';
col = dataCollection();
col = col.addData(trainDat,trainDatPred,testDat,testDatPred);
col.exportData([1 2 3 4])
controlDat = data();
controlDat = controlDat.importData([],[],squeeze(out.x),out.u(1:end-1)',trainDat.p,0:mod.ts:Tmax,mod.ts,label);
col = col.addData(controlDat);
col.exportData(5)

writematrix([(0:mod.ts:Tmax)'/3600 out.dfeed],['FilesForAgglo\results\data\dfeed',label,'.csv'])
save('col','col')