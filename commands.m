%% Normal Simulation
addpath('Model','-end') % adds the path of the C++ code

% initialise parameters 
p0 = 0.2;
psc = 1e-5;
dmax = 20;
gage = 270;
page = 2;
EC50 = 0.01;

p0 = 0.2;

tumour_volume_initial = 201; %this is the size of the tumour we want to grow before we simulate treatment (or in this case no treatment)
p = clib.Model.SeedAndGrowToStartVolumeM(p0, psc, dmax, gage, page, EC50, tumour_volume_initial); % simulates the growth of a tumour starting at 1 cell and then stopping when the tumour volume reaches tumour_volume_initial

% initialise drug injection location and concentration
xinj1 = 0;
yinj1 = 0;
C0 = 1000; 

for jj = 1:1 % does 250 different simulations of tumour growth
   
    p = clib.Model.SeedAndGrowToStartVolumeM(p0, psc, dmax, gage, page, EC50, tumour_volume_initial); % simulates the growth of a tumour starting at 1 cell and then stopping when the tumour volume reaches tumour_volume_initial
    psim = clib.Model.CreateNewParticle(p0, psc, dmax, gage, page, EC50, p); %sets the initial tumour size for the simulation as the size of the tumour "p", i.e. tumour_volume_initial
   % psim.InjectFibre(xinj1, yinj1, C0*2000/(10+1));% injects fibre at position xinj1 yinj1 of concentration C0
    
    % initialising the vectors

    TotalAout(1) = psim.ReturnDrugConcentrationAout;
    [cellpositions celltype]= cellposextractor(psim);
    cellpositions_mat{1} = cellpositions;
    celltype_mat{1} = celltype;
        
    for ii = 1:33 %simulates 33 days of tumour growth one day at a time - if possible can you also simulate for 50 and 100 days so we can see the difference between the three
        
        Tvol(ii) = psim.SimulateOneDay(1); %simulates the growth of the tumour for one day and returns its volume
        NumberTcells(ii) = psim.ReturnTotalNumberTumourCells;% returns total number of tumour cells
        NumberDeadcells(ii) = psim.ReturnTotalNumberDeadCells;% returned total number of dead cells
        NumberPSCcells(ii) = psim.ReturnTotalNumberPSCCells;% returned total number of PSC cells
        NumberHealthycells(ii) = psim.ReturnTotalNumberHealthyCells;
        Totaldrugconc(ii) = psim.ReturnDrugConcentrationDomain;% returned total number of healthy cells
        Totaldrugconcfibre(ii) = psim.ReturnDrugConcentrationinFibre; % returned total Drug conc in domain
        TotalAout(ii+1) = psim.ReturnDrugConcentrationAout; % returned total concentration released from fibre
        
        [cellpositions celltype]= cellposextractor(psim);
        cellpositions_mat{ii+1} = cellpositions;
        celltype_mat{ii+1} = celltype;
    end
    
    Tvol_mat(jj,:) = Tvol; %record the volume of the tumour over 33 days in a matrix for all the tumour growths of the 10 particles
    NumberTcells_mat(jj,:) = NumberTcells;
    NumberDeadcells_mat(jj,:) = NumberDeadcells;
    TotalDeadcells_mat(jj) = sum(NumberDeadcells);
    NumberPSCcells_mat(jj,:) = NumberPSCcells;
    NumberHealthycells_mat(jj,:) = NumberHealthycells;
    Totaldrugconc_mat(jj,:) = Totaldrugconc;
    TotalAout_mat(jj,:) = TotalAout;
    
save('Longinj.mat', 'Tvol_mat');
    jj
end
save('Longinj.mat', 'Tvol_mat');

% % plots the tumour volume of 10 simulations of the model for the same parameter values
figure
hold on 
yyaxis left
plot(1:33,Totaldrugconc_mat')%,'Color',[0.5 0.5 0.5], 'LineWidth',1)
ylabel('Drug outside fibre')
yyaxis right
plot(0:33,TotalAout_mat')%,'Color',[0.5 0.5 0.5], 'LineWidth',1)
xlabel('Time (days)')
ylabel('Fibre drug release curve')
set(gca,'FontSize',18)
title('Total drug concentration')

% % plots the tumour volume of 10 simulations of the model for the same parameter values
figure
hold on 
plot(1:33,Tvol_mat',':','Color',[0.5 0.5 0.5], 'LineWidth',1)
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)
ylabel('Cells')
title('Tumour Volume')

% % plots the tumour volume of 10 simulations of the model for the same parameter values
figure
hold on 
yyaxis left
plot(1:33,NumberTcells_mat',':','Color',[0.5 0.5 0.5], 'LineWidth',1)
ylabel('Tumour cells')
yyaxis right
plot(NumberHealthycells_mat',':','Color',[0.5 0.5 0.5], 'LineWidth',1)
xlabel('Time (days)')
ylabel('Healthy cells')
set(gca,'FontSize',18)
title('Tumour vs healthy cells')

% % plots the tumour volume of 10 simulations of the model for the same parameter values
figure
hold on 
plot(1:33,NumberDeadcells_mat',':','Color',[0.5 0.5 0.5], 'LineWidth',1)
xlabel('Time (days)')
ylabel('Cells')
set(gca,'FontSize',18)

% % plots the tumour volume of 10 simulations of the model for the same parameter values
figure
hold on 
histogram(TotalDeadcells_mat)
xlabel('Cells')
ylabel('Frequency')
set(gca,'FontSize',18)
title('Dead cells')
% 

% % plots the tumour volume of 10 simulations of the model for the same parameter values
figure
hold on 
plot(1:33,NumberPSCcells_mat',':','Color',[0.5 0.5 0.5], 'LineWidth',1)
xlabel('Time (days)')
ylabel('Cells')
set(gca,'FontSize',18)
title('PSC cells')
% 
% 
% %plots the mean and std of the 10 simulations of tumour volume
fig1 = figure
hold on
options.handle = fig1;
options.color_area = [141,160,203]/255;
options.color_line =  [141,160,203]/255;
options.alpha = 0.3;
options.error = 'std';
options.line_width = 2;
options.x_axis = 1:33;
plot_areaerrorbar(Tvol_mat,options)
xlabel('Time (days)')
ylabel('Tumour volume (mm^3)')
set(gca,'FontSize',18)

%% Plots 4 frames of tumour VCBM-PDE simulation 

figure
subplot(2,2,1)
domain_plotter(cellpositions_mat{1},celltype_mat{1})
title('Day 1')
subplot(2,2,2)
domain_plotter(cellpositions_mat{5},celltype_mat{5})
title('Day 5')
subplot(2,2,3)
domain_plotter(cellpositions_mat{10},celltype_mat{10})
title('Day 10')
subplot(2,2,4)
domain_plotter(cellpositions_mat{15},celltype_mat{15})
title('Day 15')

