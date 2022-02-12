function main(file_name)
%MAIN function
clearvars -except file_name
clc;

sim = Engine();
sim.start(file_name)
sim.solver()
sim.finish()
end