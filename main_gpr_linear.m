clear; close all; clc
% Parameters
group = 1;
run_abaqus = 1; % 1 if you want to run abaqus, 0 if no
pcNum = 5; % number of pcs to use in model

% Initialize variables
numRun = 500;
iter = 1;
initMicro = 100;
totMicro = 3125;
numTest = 1000;
remMicro = totMicro - initMicro - numTest;
numMicro = initMicro;

% Load PCA & Vf
load('eigenvectors.mat','eigenvectors');
load('pcs.mat','pcs');
load('eigenvalues.mat','eigenvalues');
load('micro_param.mat','micro_param');
microParamArray = [(1:totMicro)',micro_param(1:totMicro,:)];
pcArray = [(1:totMicro)',pcs(1:totMicro,1:pcNum)];
gpTrainArray = [];
gpTestArray = [];
availArray = pcArray;

% Initialization of sparse PC space
ptInit = pcArray(20,:);   % Initialize as point 1386 first
ptInit = pcArray(1000,:);   % Initialize as point 1386 first
gpTrainArray(end+1,:) = ptInit;

availArrayRem = availArray(:,1) == ptInit(1,1);    % remove initializing point from array
availArray(availArrayRem,:) = [];

for i = 1:(initMicro-1)
    for j = 1:size(gpTrainArray,1)
        pcDistInt = 0;
        for k = 1:pcNum
            pcDistInt = pcDistInt + (availArray(:,(k+1)).*var(k)).^2;
        end
        fprintf('pcDistInt: %4.0f\n',pcDistInt);
        pcDistTot(:,j) = pcDistInt.^0.5;
        pcDist = mean(pcDistTot,2);
    end
    pcDist = [availArray(:,1),pcDist];
    [~,idx] = sort(pcDist(:,2),'descend');
    pcDist = pcDist(idx,:);
    pcArrayTemp = pcArray(:,1) == pcDist(i,1);
    gpTrainArray(end+1,:) = pcArray(pcArrayTemp,:); 
    availArrayRem = availArray(:,1) == pcDist(i,1);
    availArray(availArrayRem,:) = [];
    if i < (initMicro-1)
        clear pcDistTot
    end
    j = 1;
end

gpTrainArrayInit = gpTrainArray;    % Set asside initial array for plotting

 scatter3(gpTrainArrayInit(:,2),gpTrainArrayInit(:,3),gpTrainArrayInit(:,4),'filled');
 xlabel('PC1');
 ylabel('PC2');
 zlabel('PC3');

%% Launch initial subset and read k results
for i = 1:numMicro
    numLabel{i} = num2str(gpTrainArray(i,1),'%04.f');   
    abqFile = strcat(numLabel{i},'_1550_JobQ');
    abqCmd = strcat('abaqus job=',abqFile,' input=',abqFile,'.inp cpus=12 int double');
    %abaqus job=${j}_1550_JobQ.inp cpus=24 int double interactive
    abqCmdPy = strcat(['abaqus cae noGUI=avg_hflux.py --',' ',numLabel{i}]);
    lckfilename = strcat(abqFile,'.lck');
    resultfile = strcat(abqFile,'_k.txt');
    
    if (run_abaqus == 1)
        if(isfile(resultfile) == 0)
            % launch ABAQUS analysis
            system(abqCmd);
            pause(20);
            
            while isfile(lckfilename)   % While ABAQUS lck file exists just wait
                pause(60);
                fprintf(strcat(['ABAQUS Run: ',numLabel{i},'\n']));
            end
            % post processing python script
            system(abqCmdPy);
            
            % Delete additional ABAQUS files
            delete(sprintf('%s.odb',abqFile))
            delete(sprintf('%s.com',abqFile))
            delete(sprintf('%s.dat',abqFile))
            delete(sprintf('%s.prt',abqFile))
            delete(sprintf('%s.sim',abqFile))
            delete(sprintf('%s.sta',abqFile))
        end
    end
    
    outputResults = strcat(numLabel{i},'_1550_JobQ_k.txt');
    
    fid = fopen(outputResults,'r');
    dataResults = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace','');
    fclose(fid);
    
    k11(i) = str2double(cell2mat(dataResults{1}(2)));
    k22(i) = str2double(cell2mat(dataResults{1}(3)));
    k33(i) = str2double(cell2mat(dataResults{1}(4)));
    kplane(i) = (k11(i)+k22(i))/2;
    Vtow(i) = microParamArray(str2double(numLabel{i}),2);
    Vmat(i) = microParamArray(str2double(numLabel{i}),4);
    Vpore(i) = microParamArray(str2double(numLabel{i}),3);
    
    fprintf('****************************************\n');
    fprintf('Initial Count: %4.0f\n',i);
end

%% Read in test dataset
load('k11t.mat','k11t');
load('k22t.mat','k22t');
load('k33t.mat','k33t');
load('gpTestArray.mat','gpTestArray');

for i = 1:numTest
    availArrayRem = availArray(:,1) == gpTestArray(i,1);
    availArray(availArrayRem,:) = [];
    
    numLabelTest{i} = num2str(gpTestArray(i,1),'%04.f');   
    outputResults = strcat(numLabelTest{i},'_1550_JobQ_k.txt');
    
    fid = fopen(outputResults,'r');
    dataResults = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace','');
    fclose(fid);
    
    VtowTest(i) = microParamArray(str2double(numLabelTest{i}),2);
    VmatTest(i) = microParamArray(str2double(numLabelTest{i}),4);
    VporeTest(i) = microParamArray(str2double(numLabelTest{i}),3);
    
    gpTestArray(i,5) = pcs(gpTestArray(i,1),4);
    gpTestArray(i,6) = pcs(gpTestArray(i,1),5);
    
    fprintf('****************************************\n');
    fprintf('Test Count: %4.0f\n',i);
    
end

%% Kernel Parameters
close all

gprMdl1 = fitrgp(gpTrainArray(1:numMicro,2:(pcNum+1)),k11(:),'BasisFunction','linear',...
    'KernelFunction','ardsquaredexponential');

gprMdl2 = fitrgp(gpTrainArray(1:numMicro,2:(pcNum+1)),k22(:),'BasisFunction','linear',...
    'KernelFunction','ardsquaredexponential');

gprMdl3 = fitrgp(gpTrainArray(1:numMicro,2:(pcNum+1)),k33(:),'BasisFunction','linear',...
    'KernelFunction','ardsquaredexponential');

sigma = [gprMdl1.Sigma,gprMdl2.Sigma,gprMdl3.Sigma];
kparams = [gprMdl1.KernelInformation.KernelParameters';...
    gprMdl2.KernelInformation.KernelParameters';gprMdl3.KernelInformation.KernelParameters'];
beta = [gprMdl1.Beta';gprMdl2.Beta';gprMdl3.Beta'];

%% Running Loop
while (iter <= numRun)
% Gaussian Process Regression
        
        gprMdl1 = fitrgp(gpTrainArray(1:numMicro,2:(pcNum+1)),k11(:),'BasisFunction','linear',...
            'KernelFunction','ardsquaredexponential','ConstantSigma',true,'Sigma',0.24,...
            'Standardize',1,'Optimizer','quasinewton');
    
        gprMdl2 = fitrgp(gpTrainArray(1:numMicro,2:(pcNum+1)),k22(:),'BasisFunction','linear',...
            'KernelFunction','ardsquaredexponential','ConstantSigma',true,'Sigma',0.23,...
            'Standardize',1,'Optimizer','quasinewton');
    
        gprMdl3 = fitrgp(gpTrainArray(1:numMicro,2:(pcNum+1)),k33(:),'BasisFunction','linear',...
            'KernelFunction','ardsquaredexponential','ConstantSigma',true,'Sigma',0.30,...
            'Standardize',1,'Optimizer','quasinewton');    
    
    [k11pred,~,ci1] = resubPredict(gprMdl1);
    [k22pred,~,ci2] = resubPredict(gprMdl2);
    [k33pred,~,ci3] = resubPredict(gprMdl3);
    
    %% Check Error on Test Set
    [k11test,~,ci1Test] = predict(gprMdl1, gpTestArray(:,2:(pcNum+1)));
    [k22test,~,ci2Test] = predict(gprMdl2, gpTestArray(:,2:(pcNum+1)));
    [k33test,~,ci3Test] = predict(gprMdl3, gpTestArray(:,2:(pcNum+1)));
    mseTest(iter,1) = loss(gprMdl1,gpTestArray(:,2:(pcNum+1)),k11t);
    mseTest(iter,2) = loss(gprMdl2,gpTestArray(:,2:(pcNum+1)),k22t);
    mseTest(iter,3) = loss(gprMdl3,gpTestArray(:,2:(pcNum+1)),k33t);
    maeTest(iter,1) = mean(abs(k11test - k11t'))/mean(k11t);
    maeTest(iter,2) = mean(abs(k22test - k22t'))/mean(k22t);
    maeTest(iter,3) = mean(abs(k33test - k33t'))/mean(k33t);    
    
    %% Record kernel function parameters and sigma
    kparams(1,:) = gprMdl1.KernelInformation.KernelParameters';
    kparams(2,:) = gprMdl2.KernelInformation.KernelParameters';
    kparams(3,:) = gprMdl3.KernelInformation.KernelParameters'; 
    kernel(iter,1:(pcNum+1),1) = gprMdl1.KernelInformation.KernelParameters';
    kernel(iter,1:(pcNum+1),2) = gprMdl2.KernelInformation.KernelParameters';
    kernel(iter,1:(pcNum+1),3) = gprMdl3.KernelInformation.KernelParameters';    
    beta(1,:) = gprMdl1.Beta';
    beta(2,:) = gprMdl2.Beta';
    beta(3,:) = gprMdl3.Beta';      
    betap(iter,1:(pcNum+1),1) = gprMdl1.Beta'; 
    betap(iter,1:(pcNum+1),2) = gprMdl2.Beta'; 
    betap(iter,1:(pcNum+1),3) = gprMdl3.Beta';  
    sigma(iter,1) = gprMdl1.Sigma;
    sigma(iter,2) = gprMdl2.Sigma;
    sigma(iter,3) = gprMdl3.Sigma;
    
    %% Find highest variance location and launch next run
    gpResultArray = [str2double(numLabel)', k11pred, ci1(:,1), ci1(:,2), (ci1(:,2) - ci1(:,1)),...
        k22pred, ci2(:,1), ci2(:,2), (ci2(:,2) - ci2(:,1)), k33pred, ci3(:,1), ci3(:,2), (ci3(:,2) - ci3(:,1))];
    
    [k11var,~,cipred1] = predict(gprMdl1, availArray(1:remMicro,2:(pcNum+1)));
    [k22var,~,cipred2] = predict(gprMdl2, availArray(1:remMicro,2:(pcNum+1)));
    [k33var,~,cipred3] = predict(gprMdl3, availArray(1:remMicro,2:(pcNum+1)));
    
    ciRange = [abs(cipred1(1:remMicro,1)-cipred1(1:remMicro,2)),abs(cipred2(1:remMicro,1)-cipred2(1:remMicro,2)),abs(cipred3(1:remMicro,1)-cipred3(1:remMicro,2))];
    ciRange = mean(ciRange,2);
%    ciRange = [abs(cipred3(1:remMicro,1)-cipred3(1:remMicro,2))];
    ciRange = [availArray(1:remMicro,1),ciRange];
    
    [~,idx] = sort(ciRange(:,2),'descend');
    ciRange = ciRange(idx,:);
    
    % Print out next micros to run and remove them from availArray
    fid = fopen('next_simulation_ind.txt','w');
    for i = 1:group
        fprintf(fid,'%04.f\n',ciRange(i,1));
        numLabel{group*(iter-1)+initMicro+i} = num2str(ciRange(i,1),'%04.f');  % add microstructures to read results to list
        availArrayRem = availArray(:,1) == ciRange(i,1);    % remove top 5 closest PC score micros from array
        availArray(availArrayRem,:) = [];
        
        % look up top micro and add to gpTrainArray
        pcArrayTemp = pcArray(:,1) == ciRange(i,1);
        gpTrainArray(end+1,:) = pcArray(pcArrayTemp,:);    
        fprintf('Next Point Index: %4.0f\n',ciRange(i,1));
    end
    fclose(fid);
    fprintf('****************************************\n');
    
    %% Wait until analysis is complete and read k results
    for i = 1:group
        j = size(numLabel,2) - group + i;
        abqFile = strcat(numLabel{j},'_1550_JobQ');
        abqCmd = strcat('abaqus job=',abqFile,' input=',abqFile,'.inp cpus=12 int double');
        %abaqus job=${j}_1550_JobQ.inp cpus=24 int double interactive
        abqCmdPy = strcat(['abaqus cae noGUI=avg_hflux.py --',' ',numLabel{j}]);
        lckfilename = strcat(abqFile,'.lck');
        resultfile = strcat(abqFile,'_k.txt');
        
        if (run_abaqus == 1)
            if(isfile(resultfile) == 0)
                % launch ABAQUS analysis
                system(abqCmd);
                pause(20);
                
                while isfile(lckfilename)   % While ABAQUS lck file exists just wait
                    pause(60);
                    fprintf(strcat(['ABAQUS Run: ',numLabel{j},'\n']));
                end
                % post processing python script
                system(abqCmdPy);
                
                % Delete additional ABAQUS files
                delete(sprintf('%s.odb',abqFile))
                delete(sprintf('%s.com',abqFile))
                delete(sprintf('%s.dat',abqFile))
                delete(sprintf('%s.prt',abqFile))
                delete(sprintf('%s.sim',abqFile))
                delete(sprintf('%s.sta',abqFile))              
            end
        end
        
        outputResults = strcat(numLabel{j},'_1550_JobQ_k.txt');
        
        fid = fopen(outputResults,'r');
        dataResults = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace','');
        fclose(fid);
        
        k11(j) = str2double(cell2mat(dataResults{1}(2)));
        k22(j) = str2double(cell2mat(dataResults{1}(3)));
        k33(j) = str2double(cell2mat(dataResults{1}(4)));
        kplane(j) = (k11(j)+k22(j))/2;
        Vtow(j) = microParamArray(str2double(numLabel{j}),2);
        Vmat(j) = microParamArray(str2double(numLabel{j}),4);
        Vpore(j) = microParamArray(str2double(numLabel{j}),3);
%        [kmax_11(j), kh2l_11(j), khj_11(j), khjme_11(j)] = analyticalK1(Vpore(j),Vmat(j),Vtow(j),km,(km*kf)^0.5);
%        [kmax_33(j), kh2l_33(j), khj_33(j), khjme_33(j)] = analyticalK3(Vpore(j),Vmat(j),Vtow(j),km,kft);
    end
    
    %% Incrememnt counter
    numMicro = numMicro + group;
    remMicro = remMicro - group;
    iter = iter + 1;
    fprintf('****************************************\n');
    fprintf('Count: %4.0f\n',iter);
    close all
    if (remMicro == 0)
        break
    end
%     %% Plot hyperparameters
%     %figure('Renderer', 'painters', 'Position', [200 200 900 800])
%     for i = 1:(pcNum+1)
%         subplot(pcNum+1,2,i);
%         hold on
%         box on
%         %plot(1:size(kernel,1),kernel(:,i,1));
%         %plot(1:size(kernel,1),kernel(:,i,2));
%         plot(1:size(kernel,1),kernel(:,i,3));
%         xlabel('Iteration');
%         ylabel('\lambda_1');
%         hold off
%         if i <= pcNum
%             ylabel(strcat(['\lambda_',num2str(i)]));
%         else
%             ylabel('\sigma_f');
%         end
%     end
%     
%     % Beta
%     for i = 1:size(gprMdl1.Beta,1)
%         subplot(pcNum+1,2,i+pcNum+1);
%         hold on
%         box on
%         %plot(1:size(betap,1),betap(:,i,1));
%         %plot(1:size(betap,1),betap(:,i,2));
%         plot(1:size(betap,1),betap(:,i,3));
%         xlabel('Iteration');
%         hold off
%         if i == 1
%             ylabel('\beta');
%         else
%             ylabel(strcat(['\beta_',num2str(i-1)]));
%         end
%     end
%     hL = legend({'k_1_1','k_2_2','k_3_3'});
%     % Programatically move the Legend
%     newPosition = [0.835 0.85 0.05 0.05];
%     newUnits = 'normalized';
%     set(hL,'Position', newPosition,'Units', newUnits);
    
end

%% Parity Plot
figure('Renderer', 'painters', 'Position', [300 300 1500 400])
subplot(1,3,1);
hold on
errorbar(k11(1:size(k11pred,1)),k11pred,(gpResultArray(:,5)/2),'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
errorbar(k11t,k11test,((ci1Test(:,2)-ci1Test(:,1))/2),'o','MarkerSize',4,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
plot([min([k11pred,k11(1:size(k11pred,1))']),max([k11pred,k11(1:size(k11pred,1))'])],...
    [min([k11pred,k11(1:size(k11pred,1))']),max([k11pred,k11(1:size(k11pred,1))'])],'k','LineWidth',2)
xlabel('k_1_1 Actual');
ylabel('k_1_1 Predicted');
box on
axis tight
hold off

subplot(1,3,2);
hold on
errorbar(k22(1:size(k22pred,1)),k22pred,(gpResultArray(:,9)/2),'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
errorbar(k22t,k22test,((ci2Test(:,2)-ci2Test(:,1))/2),'o','MarkerSize',4,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
plot([min([k22pred,k22(1:size(k22pred,1))']),max([k22pred,k22(1:size(k22pred,1))'])],...
    [min([k22pred,k22(1:size(k22pred,1))']),max([k22pred,k22(1:size(k22pred,1))'])],'k','LineWidth',2)
xlabel('k_2_2 Actual');
ylabel('k_2_2 Predicted');
box on
axis tight
hold off

subplot(1,3,3);
hold on
errorbar(k33(1:size(k33pred,1)),k33pred,(gpResultArray(:,13)/2),'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
errorbar(k33t,k33test,((ci3Test(:,2)-ci3Test(:,1))/2),'o','MarkerSize',4,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
plot([min([k33pred,k33(1:size(k33pred,1))']),max([k33pred,k33(1:size(k33pred,1))'])],...
    [min([k33pred,k33(1:size(k33pred,1))']),max([k33pred,k33(1:size(k33pred,1))'])],'k','LineWidth',2)
xlabel('k_3_3 Actual');
ylabel('k_3_3 Predicted');
box on
axis tight
hold off
hL = legend({'Train','Test'});
% Programatically move the Legend
newPosition = [0.92 0.86 0.05 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
%savefig('parity_plot.fig')
%close

%% Analytical Results
towWidth = micro_param(:,4);
towHeight = micro_param(:,5);

km = 65;
kf = 32.58;
kft = 5.1;
kp = 0.001;
%
% % Calculate for available
% for i = 1:size(availArray,1)   
%     numLabelA{i} = num2str(availArray(i,1),'%04.f');   
%     outputResults = strcat(numLabelA{i},'_1550_JobQ_k.txt');
%     
%     fid = fopen(outputResults,'r');
%     dataResults = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace','');
%     fclose(fid);
%     
%     k11A(i) = str2double(cell2mat(dataResults{1}(2)));
%     k22A(i) = str2double(cell2mat(dataResults{1}(3)));
%     k33A(i) = str2double(cell2mat(dataResults{1}(4)));
%     VtowAvail(i) = microParamArray(str2double(numLabelA{i}),2);
%     VmatAvail(i) = microParamArray(str2double(numLabelA{i}),4);
%     VporeAvail(i) = microParamArray(str2double(numLabelA{i}),3);
%     
%     fprintf('****************************************\n');
%     fprintf('Avail Count: %4.0f\n',i);
%     
% end
% 
% k11t = [k11t k11A];
% k22t = [k22t k22A];
% k33t = [k33t k33A];
% VporeTest = [VporeTest VporeAvail];
% VmatTest = [VmatTest VmatAvail];
% VtowTest = [VtowTest VtowAvail];
%%
% Find best beta fit
beta = 0:0.01:10;
Pp = 0:0.01:1;
for i = 1:length(beta)
    for j = 1:length(Pp)
        [khj_33, khjme_33, kh2l_33] = analyticalK3(Vpore,Vmat,Vtow,km,kf,kft,beta(i),Pp(j));
        k33h2lError = (kh2l_33(1:length(k33pred)) - k33(1:length(k33pred)))./k33(1:length(k33pred));
        err_norm(i,j) = norm(k33h2lError);
    end
end

val = min(min(err_norm));
[xind,yind]=find(err_norm==val);

[khj_33, khjme_33, kh2l_33] = analyticalK3(Vpore,Vmat,Vtow,km,kf,kft,beta(xind),Pp(yind));
[khj_33T, khjme_33T, kh2l_33T] = analyticalK3(VporeTest,VmatTest,VtowTest,km,kf,kft,beta(xind),Pp(yind));
h2l_MAE = mean(abs(kh2l_33T - k33t))/mean(k33t);
fprintf('h2l_MAE: %4.8f\n',h2l_MAE);

kEffs = Vmat.*km + Vtow.*kft + Vpore.*kp;
kEffp = 1./(Vmat./km + Vtow./kft + Vpore./(kp+0.001));

% Plot analytical results vs model
k33predError = (k33pred' - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33predLError = (gpResultArray(:,11)' - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33predUError = (gpResultArray(:,12)' - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33predCIError = k33predUError - k33predLError;

k33predErrorT = (k33test' - k33t(1:length(k33test)))./k33t(1:length(k33test));

k33hjError = (khj_33(1:length(k33pred)) - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33hjmeError = (khjme_33(1:length(k33pred)) - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33sError = (kEffs(1:length(k33pred)) - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33pError = (kEffp(1:length(k33pred)) - k33(1:length(k33pred)))./k33(1:length(k33pred));
k33h2lError = (kh2l_33(1:length(k33pred)) - k33(1:length(k33pred)))./k33(1:length(k33pred));

k33h2lErrorT = (kh2l_33T(1:length(k33t)) - k33t(1:length(k33t)))./k33t(1:length(k33t));

%figure('Renderer', 'painters', 'Position', [50 50 1500 400])
%subplot(1,2,1);
figure
hold on
box on
grid on
scatter(VporeTest(1:length(k33predErrorT)),100.*k33predErrorT,16,'o','filled');
%errorbar(Vpore(1:length(k33predError)),k33predError.*100,(k33predCIError/2).*100,'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
%scatter(Vpore(1:length(k33predError)),100.*k33hjError,'o','filled');
%scatter(Vpore(1:length(k33predError)),100.*k33hjmeError,'o','filled');
%scatter(Vpore(1:length(k33predError)),100.*k33sError,'o','filled');
scatter(VporeTest,100.*k33h2lErrorT,16,'o','filled',...
    'MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerEdgeColor',[0.25, 0.25, 0.25]);
plot([0 0.4],[0 0],'k--');
ylim([-25 25]);
xlabel('Vpore');
ylabel('Error (%)');
hold off
legend('GPR','H2L','Location','SouthEast');

% subplot(1,2,2);
% hold on
% box on
% scatter(Vpore(1:length(k33pred)),k33pred,8,'o','filled');
% scatter(Vpore(1:length(k33pred)),k33(1:length(k33pred))',8,'o','filled');
% %scatter(Vpore(1:length(k33pred)),khj_33(1:length(k33pred)),'o','filled');
% %scatter(Vpore(1:length(k33pred)),khjme_33(1:length(k33pred)),'o','filled');
% scatter(Vpore(1:length(k33pred)),kh2l_33(1:length(k33pred)),8,'o','filled');
% xlabel('Vpore');
% ylabel('k_3_3 (W/mK)');
% legend('Predicted','Actual','H2L','Location','NorthEastOutside');
% hold off

%savefig('pore_error.fig')
%close

%% Parity Plot with Analytical
figure('Renderer', 'painters', 'Position', [300 300 1500 400])
subplot(1,3,1);
hold on
errorbar(k11(1:size(k11pred,1)),k11pred,(gpResultArray(:,5)/2),'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
errorbar(k11t,k11test,((ci1Test(:,2)-ci1Test(:,1))/2),'o','MarkerSize',4,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
%scatter(k11,kEffp,'filled');
%scatter(k11,kEffs,'filled');
plot([min([k11pred,k11(1:size(k11pred,1))']),max([k11pred,k11(1:size(k11pred,1))'])],...
    [min([k11pred,k11(1:size(k11pred,1))']),max([k11pred,k11(1:size(k11pred,1))'])],'k','LineWidth',2)
xlabel('k_1_1 Actual');
ylabel('k_1_1 Predicted');
box on
axis tight
hold off
legend('Train','Test','Location','NorthWest')

subplot(1,3,2);
hold on
errorbar(k22(1:size(k22pred,1)),k22pred,(gpResultArray(:,9)/2),'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
errorbar(k22t,k22test,((ci2Test(:,2)-ci2Test(:,1))/2),'o','MarkerSize',4,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
%scatter(k22,kEffp,'filled');
%scatter(k22,kEffs,'filled');
plot([min([k22pred,k22(1:size(k22pred,1))']),max([k22pred,k22(1:size(k22pred,1))'])],...
    [min([k22pred,k22(1:size(k22pred,1))']),max([k22pred,k22(1:size(k22pred,1))'])],'k','LineWidth',2)
xlabel('k_2_2 Actual');
ylabel('k_2_2 Predicted');
box on
axis tight
hold off

subplot(1,3,3);
hold on
errorbar(k33(1:size(k33pred,1)),k33pred,(gpResultArray(:,13)/2),'o','MarkerSize',4,'MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410]);
errorbar(k33t,k33test,((ci3Test(:,2)-ci3Test(:,1))/2),'o','MarkerSize',4,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
%scatter(k33,kEffp,'filled');
scatter(k33,kh2l_33,10,'filled','MarkerFaceColor',[0.9290, 0.6940, 0.1250],'MarkerEdgeColor',[0.25, 0.25, 0.25]);
plot([min([k33pred,k33(1:size(k33pred,1))']),max([k33pred,k33(1:size(k33pred,1))'])],...
    [min([k33pred,k33(1:size(k33pred,1))']),max([k33pred,k33(1:size(k33pred,1))'])],'k','LineWidth',2)
xlabel('k_3_3 Actual');
ylabel('k_3_3 Predicted');
%legend('Train','Test','Location','NorthWest');
box on
axis tight
hold off
%legend('Train','Test','H2L','Location','NorthWest')

hL = legend({'Train','Test','H2L'});
% Programatically move the Legend
newPosition = [0.92 0.86 0.05 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

%% MSE/MAE/Max
RSME = [sqrt(mseTest(:,1))./mean(k11),sqrt(mseTest(:,2))./mean(k22),sqrt(mseTest(:,3))./mean(k33)];
MAE = [maeTest(:,1),maeTest(:,2),maeTest(:,3)];

figure
hold on
box on
plot((1:size(mseTest,1))+initMicro,sqrt(mseTest(:,1))./mean(k11),'--','Color',[0, 0.4470, 0.7410]);
plot((1:size(mseTest,1))+initMicro,sqrt(mseTest(:,2))./mean(k22),'--','Color',[0.8500, 0.3250, 0.0980]);
plot((1:size(mseTest,1))+initMicro,sqrt(mseTest(:,3))./mean(k33),'--','Color',[0.4660, 0.6740, 0.1880]);
plot((1:size(mseTest,1))+initMicro,maeTest(:,1),'-','Color',[0, 0.4470, 0.7410]);
plot((1:size(mseTest,1))+initMicro,maeTest(:,2),'-','Color',[0.8500, 0.3250, 0.0980]);
plot((1:size(mseTest,1))+initMicro,maeTest(:,3),'-','Color',[0.4660, 0.6740, 0.1880]);
ylim([0 0.1]);
xlabel('Iteration');
ylabel('Normalized Error');
legend('RMSE k_1_1','RMSE k_2_2','RMSE k_3_3','MAE k_1_1','MAE k_2_2','MAE k_3_3');
hold off
%savefig('RMSE_MAE.fig')
close
figure
hold on
box on
plot((1:size(mseTest,1))+initMicro,sqrt(mseTest(:,1))./mean(k11),'-','Color',[0, 0.4470, 0.7410]);
plot((1:size(mseTest,1))+initMicro,sqrt(mseTest(:,2))./mean(k22),'-','Color',[0.8500, 0.3250, 0.0980]);
plot((1:size(mseTest,1))+initMicro,sqrt(mseTest(:,3))./mean(k33),'-','Color',[0.4660, 0.6740, 0.1880]);
ylim([0 0.1]);
xlabel('GPR Iteration');
ylabel('RMSE');
legend('k_1_1','k_2_2','k_3_3');
hold off
%savefig('RMSE.fig')
close

figure
hold on
box on
grid on
plot((1:size(mseTest,1))+initMicro,maeTest(:,1).*100,'-','Color',[0, 0.4470, 0.7410],'LineWidth',2);
plot((1:size(mseTest,1))+initMicro,maeTest(:,2).*100,'-','Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
plot((1:size(mseTest,1))+initMicro,maeTest(:,3).*100,'-','Color',[0.4660, 0.6740, 0.1880],'LineWidth',2);
ylim([0 10]);
xlabel('Number of Training Microstructures');
ylabel('Normalized Test MAE (%)');
legend('k_1_1','k_2_2','k_3_3');
hold off
%savefig('MAE.fig')
%close
% 
%% Plot hyperparameters Norm and angle
for i = 1:size(kernel,1)
    norm_length1(i) = norm([kernel(i,:,1),betap(i,:,1),sigma(i,1)]);
    norm_length2(i) = norm([kernel(i,:,2),betap(i,:,2),sigma(i,1)]);
    norm_length3(i) = norm([kernel(i,:,3),betap(i,:,3),sigma(i,1)]);
end

% for i = 1:size(kernel,1)
%     norm_length1(i) = norm([betap(i,:,1),sigma(i,1)]);
%     norm_length2(i) = norm([betap(i,:,2),sigma(i,1)]);
%     norm_length3(i) = norm([betap(i,:,3),sigma(i,1)]);
% end

for i = 1:size(kernel,1)-1
    vec1 = [kernel(i,:,1),betap(i,:,1),sigma(i,1)];
    vec2 = [kernel(i+1,:,1),betap(i+1,:,1),sigma(i+1,1)];
    angleDiff1(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
    vec1 = [kernel(i,:,2),betap(i,:,2),sigma(i,2)];
    vec2 = [kernel(i+1,:,2),betap(i+1,:,2),sigma(i+1,2)];
    angleDiff2(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
    vec1 = [kernel(i,:,3),betap(i,:,3),sigma(i,3)];
    vec2 = [kernel(i+1,:,3),betap(i+1,:,3),sigma(i+1,3)];
    angleDiff3(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);    
end

% for i = 1:size(kernel,1)-1
%     vec1 = [betap(i,:,1),sigma(i,1)];
%     vec2 = [betap(i+1,:,1),sigma(i+1,1)];
%     angleDiff1(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
%     vec1 = [betap(i,:,2),sigma(i,2)];
%     vec2 = [betap(i+1,:,2),sigma(i+1,2)];
%     angleDiff2(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
%     vec1 = [betap(i,:,3),sigma(i,3)];
%     vec2 = [betap(i+1,:,3),sigma(i+1,3)];
%     angleDiff3(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);    
% end

figure('Renderer', 'painters', 'Position', [50 50 1600 800])
subplot(2,3,1)
hold on
box on
%plot((1:length(norm_length1))+initMicro,norm_length1,'o')
scatter((1:length(norm_length1))+initMicro,norm_length1,6,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
%semilogy((1:length(norm_length1))+initMicro,norm_length1,'o')
hold off
ylabel('L2 Norm of Hyperparameters');
%xlabel('Number of Training Microstructures');
%title('k_1_1');

subplot(2,3,2)
hold on
box on
%plot((1:numRun)+initMicro,norm_length2,'o')
scatter((1:length(norm_length1))+initMicro,norm_length2,6,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
hold off
%xlabel('Number of Training Microstructures');
%title('k_2_2');

subplot(2,3,3)
hold on
box on
%plot((1:numRun)+initMicro,norm_length3,'o')
scatter((1:length(norm_length1))+initMicro,norm_length3,6,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
%xlabel('Number of Training Microstructures');
%title('k_3_3');

subplot(2,3,4)
hold on
box on
%plot((1:numRun-1)+initMicro,angleDiff1)
scatter((1:numRun-1)+initMicro,angleDiff1,6,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
ylim([0 90]);
hold off
%ylim([0 14]);
ylabel('\theta (Degrees)');
xlabel('Number of Training Microstructures');
%title('k_1_1');

subplot(2,3,5)
hold on
box on
%plot((1:numRun-1)+initMicro,angleDiff2,'o')
scatter((1:numRun-1)+initMicro,angleDiff2,6,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
ylim([0 90]);
hold off
xlabel('Number of Training Microstructures');
%title('k_2_2');

subplot(2,3,6)
hold on
box on
%plot((1:numRun-1)+initMicro,angleDiff3,'o')
scatter((1:numRun-1)+initMicro,angleDiff3,6,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
ylim([0 90]);
hold off
xlabel('Number of Training Microstructures');
%title('k_3_3');

savefig('hyp_stabilization.fig')
%close;
%%
%figure('Renderer', 'painters', 'Position', [100 100 500 800])
%subplot(2,1,1)
figure
hold on
box on
grid on
scatter((1:length(norm_length1))+initMicro,norm_length1,8,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
scatter((1:length(norm_length1))+initMicro,norm_length2,8,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
scatter((1:length(norm_length1))+initMicro,norm_length3,8,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
ylim([0 1000]);
ylabel('L2 Norm of Hyperparameters');
xlabel('Number of Training Microstructures');
legend('k_1_1','k_2_2','k_3_3');

%subplot(2,1,2)
figure
grid on
hold on
box on
scatter((1:numRun-1)+initMicro,angleDiff1,8,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
scatter((1:numRun-1)+initMicro,angleDiff2,8,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
scatter((1:numRun-1)+initMicro,angleDiff3,8,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
ylim([0 90]);
ylabel('\theta (Degrees)');
xlabel('Number of Training Microstructures');
legend('k_1_1','k_2_2','k_3_3');

% hL = legend({'k_1_1','k_2_2','k_3_3'});
% % Programatically move the Legend
% newPosition = [0.92 0.83 0.05 0.05];
% newUnits = 'normalized';
% set(hL,'Position', newPosition,'Units', newUnits);
%%
for i = 1:size(kernel,1)
    norm_length1(i) = norm([kernel(i,:,1),betap(i,:,1),sigma(i,1)]);
    norm_length2(i) = norm([kernel(i,:,2),betap(i,:,2),sigma(i,1)]);
    norm_length3(i) = norm([kernel(i,:,3),betap(i,:,3),sigma(i,1)]);
end

for i = 1:size(kernel,1)-1
    vec1 = [kernel(i,:,1),betap(i,:,1),sigma(i,1)];
    vec2 = [kernel(i+1,:,1),betap(i+1,:,1),sigma(i+1,1)];
    angleDiff1(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
    vec1 = [kernel(i,:,2),betap(i,:,2),sigma(i,2)];
    vec2 = [kernel(i+1,:,2),betap(i+1,:,2),sigma(i+1,2)];
    angleDiff2(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
    vec1 = [kernel(i,:,3),betap(i,:,3),sigma(i,3)];
    vec2 = [kernel(i+1,:,3),betap(i+1,:,3),sigma(i+1,3)];
    angleDiff3(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);    
end
% Hyperparameters
figure
hold on
box on
grid on
scatter((1:length(norm_length1))+initMicro,norm_length1,8,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
scatter((1:length(norm_length1))+initMicro,norm_length2,8,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
scatter((1:length(norm_length1))+initMicro,norm_length3,8,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
ylim([0 1000]);
ylabel('L2 Norm of Hyperparameters');
xlabel('Number of Training Microstructures');
legend('k_1_1','k_2_2','k_3_3');

figure
grid on
hold on
box on
scatter((1:numRun-1)+initMicro,angleDiff1,8,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
scatter((1:numRun-1)+initMicro,angleDiff2,8,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
scatter((1:numRun-1)+initMicro,angleDiff3,8,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
ylim([0 90]);
ylabel('\theta Hyperparameters (Degrees)');
xlabel('Number of Training Microstructures');
legend('k_1_1','k_2_2','k_3_3');

for i = 1:size(kernel,1)
    norm_length1(i) = norm([betap(i,:,1)]);
    norm_length2(i) = norm([betap(i,:,2)]);
    norm_length3(i) = norm([betap(i,:,3)]);
end

for i = 1:size(kernel,1)-1
    vec1 = [betap(i,:,1)];
    vec2 = [betap(i+1,:,1)];
    angleDiff1(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
    vec1 = [betap(i,:,2)];
    vec2 = [betap(i+1,:,2)];
    angleDiff2(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);
    vec1 = [betap(i,:,3)];
    vec2 = [betap(i+1,:,3)];
    angleDiff3(i) = acos((vec1*vec2')/(norm(vec1)*norm(vec2)))*(180/pi);    
end
% Parameters
figure
hold on
box on
grid on
scatter((1:length(norm_length1))+initMicro,norm_length1,8,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
scatter((1:length(norm_length1))+initMicro,norm_length2,8,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
scatter((1:length(norm_length1))+initMicro,norm_length3,8,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
ylim([0 45]);
ylabel('L2 Norm of Parameters');
xlabel('Number of Training Microstructures');
legend('k_1_1','k_2_2','k_3_3');

figure
grid on
hold on
box on
scatter((1:numRun-1)+initMicro,angleDiff1,8,'filled',...
    'MarkerFaceColor',[0, 0.4470, 0.7410])
scatter((1:numRun-1)+initMicro,angleDiff2,8,'filled',...
    'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
scatter((1:numRun-1)+initMicro,angleDiff3,8,'filled',...
    'MarkerFaceColor',[0.4660, 0.6740, 0.1880])
hold off
%ylim([0 90]);
ylabel('\theta Parameters (Degrees)');
xlabel('Number of Training Microstructures');
legend('k_1_1','k_2_2','k_3_3');
%% Plot hyperparameters

figure('Renderer', 'painters', 'Position', [200 200 900 800])
for i = 1:(pcNum+1)
    subplot(pcNum,3,i);
    hold on
    box on
    plot(1:size(kernel,1),kernel(:,i,1));
    plot(1:size(kernel,1),kernel(:,i,2));
    plot(1:size(kernel,1),kernel(:,i,3));
    xlabel('Iteration');
    ylabel('\lambda_1');
    hold off
    %if i <= 3 && i > 1
    %    ylim([0 30]);
    %end
    if i <= pcNum
        ylabel(strcat(['\lambda_',num2str(i)]));
    else
        ylabel('\sigma_f');
    end
end

% Beta
for i = 1:size(gprMdl1.Beta,1)
    subplot(pcNum,3,i+pcNum+1);
    hold on
    box on
    plot(1:size(betap,1),betap(:,i,1));
    plot(1:size(betap,1),betap(:,i,2));
    plot(1:size(betap,1),betap(:,i,3));
    xlabel('Iteration');
    hold off
    ylabel(strcat(['\beta_',num2str(i)]));
end

hL = legend({'k_1_1','k_2_2','k_3_3'});
% Programatically move the Legend
newPosition = [0.93 0.86 0.05 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

subplot(pcNum,3,9)
hold on
box on
plot(1:size(sigma,1),sigma(:,1));
plot(1:size(sigma,1),sigma(:,2));
plot(1:size(sigma,1),sigma(:,3));
hold off
ylabel('\sigma_n');
%savefig('hyp_legend.fig')

%% Plot PC space of test/train split
figure
hold on
box on
grid on
scatter3(pcs(:,1),pcs(:,2),pcs(:,3),10,'filled','k');
scatter3(gpTestArray(:,2),gpTestArray(:,3),gpTestArray(:,4),14,'filled');
hold off
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
legend('Avail','Test');

%% Save evolution of hyperparameters and RSME

seqDOE = {RSME,MAE,sigma,betap,kernel};
save('seqDOE.mat','seqDOE');
