%% Question 1: Repeat the same analysis for the other 3 lithologies: Generate 
% 3 further masks and generate GRIDobj of (1)elevation, (2)gradient, (3)aspect 
% and (4)curvature for each of the three masks. Generate four plots that show 
% the 4 sub-DEMs only. (40 points)
% add path 
addpath('/Users/kingsleynmoye/Desktop/School/Analysis of DEM/lab04/topotoolbox-master_2021');
addpath('/Users/kingsleynmoye/Desktop/School/Analysis of DEM/lab04/topotoolbox-master_2021/topotoolbox-master');
addpath('/Users/kingsleynmoye/Desktop/School/Analysis of DEM/lab04/DEMs_Lithology');

% Load and display the DEM
DEM = GRIDobj('Pozo_DTM_noveg_UTM11_NAD83.tif');
%Gradient
G = gradient8(DEM,'deg'); 
%Aspect
A = aspect(DEM);
%Curvature
C = curvature(DEM,'profc');

%Load Lithology Map
Lithology = GRIDobj('lithology.tif');  

%Resample DEM
L_res = resample(Lithology,DEM);

%generating a new GRIDobj equal to our resampled lithologicl GRIDobj
L1 = L_res;
L2 = L_res;
L3 = L_res;
L4 = L_res;
% The lithology values for 1
L1.Z(L1.Z >1) = 0;
L1.Z(isnan(L1.Z)) = 0;
% The lithology values for 2
L2.Z(L2.Z >2) = 0;
L2.Z(L2.Z <2) = 0;
L2.Z(isnan(L2.Z)) = 0;
% The lithology values for 3
L3.Z(L3.Z >3) = 0;
L3.Z(L3.Z <3) = 0;
L3.Z(isnan(L3.Z)) = 0;

% The lithology values for 4
L4.Z(L4.Z <4) = 0;
L4.Z(isnan(L4.Z)) = 0;

DEM_L1 = crop(DEM,L1,NaN); % Extract L1 from Elevation GRIDobj
G_L1 = crop(G,L1,NaN);  % Extract L1 from Gradient GRIDobj
A_L1 = crop(A,L1,NaN);  % Extract L1 from Aspect GRIDobj
C_L1 = crop(C,L1,NaN);  % Extract L1 from Curvature GRIDobj

DEM_L2 = crop(DEM,L2,NaN); 
G_L2 = crop(G,L2,NaN);  
A_L2 = crop(A,L2,NaN);  
C_L2 = crop(C,L2,NaN);  

DEM_L3 = crop(DEM,L3,NaN); 
G_L3 = crop(G,L3,NaN);  
A_L3 = crop(A,L3,NaN);  
C_L3 = crop(C,L3,NaN);  

DEM_L4 = crop(DEM,L4,NaN); 
G_L4 = crop(G,L4,NaN);  
A_L4 = crop(A,L4,NaN);  
C_L4 = crop(C,L4,NaN);  

figure
imageschs(DEM_L1,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L1 (Sub DEM)', 'fontsize',14,'fontweight','bold')

figure
imageschs(DEM_L2,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L2 (Sub DEM)', 'fontsize',14,'fontweight','bold')

figure
imageschs(DEM_L3,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L3 (Sub DEM)', 'fontsize',14,'fontweight','bold')


figure
imageschs(DEM_L4,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L4 (Sub DEM)', 'fontsize',14,'fontweight','bold')

figure 
subplot (2,2,1)
imageschs(DEM_L1,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L1 (Sub DEM)')
subplot (2,2,2)
imageschs(DEM_L2,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L2 (Sub DEM)')
subplot (2,2,3)
imageschs(DEM_L3,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L3 (Sub DEM)')
subplot (2,2,4)
imageschs(DEM_L4,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
hold on
title('Lithology L4 (Sub DEM)')

sgtitle('Lithology L1,L2,L3,L4 (sub DEM)','fontsize',14,'fontweight','bold')

figure
imageschs(A_L4,[],'colorbarylabel','Elevation [m]', 'ticklabels','none');
%%
subplot(4,4,1)
imageschs(DEM_L1,DEM_L1,'caxis',[0 350],'ticklabels','none', 'colorbar', false);
hold on 
title('L1','FontSize',14, 'fontweight','bold')
ylabel('Elevation','FontSize',14, 'fontweight','bold')


subplot(4,4,2)
imageschs(DEM_L2,DEM_L2,'caxis',[0 350],'ticklabels','none', 'colorbar', false);
hold on 
title('L2','FontSize',14, 'fontweight','bold')


subplot(4,4,3)
imageschs(DEM_L3,DEM_L3,'caxis',[0 350],'ticklabels','none', 'colorbar', false);
hold on 
title('L3','FontSize',14, 'fontweight','bold')


subplot(4,4,4)
imageschs(DEM_L4,DEM_L4,'caxis',[0 350],'colorbarylabel','Elevation [m]','ticklabels','none');
hold on 
title('L4','FontSize',14, 'fontweight','bold')


subplot(4,4,5)
imageschs(DEM_L1,G_L1,'caxis',[0 50],'ticklabels','none', 'colorbar', false);
hold on
ylabel('Gradient','FontSize',14, 'fontweight','bold')

subplot(4,4,6)
imageschs(DEM_L2,G_L2,'caxis',[0 50],'ticklabels','none', 'colorbar', false);


subplot(4,4,7)
imageschs(DEM_L3,G_L3,'caxis',[0 50],'ticklabels','none', 'colorbar', false);


subplot(4,4,8)
imageschs(DEM_L4,G_L4,'caxis',[0 50],'colorbarylabel','Gradient [deg]','ticklabels','none');


subplot(4,4,9)
imageschs(DEM_L1,A_L1,'caxis',[0 360],'ticklabels','none', 'colorbar', false);
hold on
ylabel('Aspect','FontSize',14, 'fontweight','bold')


subplot(4,4,10)
imageschs(DEM_L2,A_L2,'caxis',[0 360],'ticklabels','none', 'colorbar', false);


subplot(4,4,11)
imageschs(DEM_L3,A_L3,'caxis',[0 360],'ticklabels','none', 'colorbar', false);


subplot(4,4,12)
imageschs(DEM_L4,A_L4,'caxis',[0 360],'colorbarylabel','Aspect [deg]','ticklabels','none');


subplot(4,4,13)
imageschs(DEM_L1,C_L1,'caxis',[-0.1 0.1],'ticklabels','none', 'colorbar', false);
hold on
ylabel('Curvature','FontSize',14, 'fontweight','bold')


subplot(4,4,14)
imageschs(DEM_L2,C_L2,'caxis',[-0.1 0.1],'ticklabels','none', 'colorbar', false);

 
subplot(4,4,15)
imageschs(DEM_L3,C_L3,'caxis',[-0.1 0.1],'ticklabels','none', 'colorbar', false);


subplot(4,4,16)
imageschs(DEM_L4,C_L4,'caxis',[-0.1 0.1],'colorbarylabel','Curvature [-]','ticklabels','none');



sgtitle('Elevation, Gradient, Aspect and Curvature of Lithologies L1, L2, L3 and L4','FontSize',14, 'fontweight','bold')
%%
subplot(4,4,1)
imageschs(DEM_L1,[],'caxis',[0 350],'ticklabels','none', 'colorbar', false);
hold on 
title('L1','FontSize',14, 'fontweight','bold')
ylabel('Elevation','FontSize',14, 'fontweight','bold')


subplot(4,4,2)
imageschs(DEM_L2,[],'caxis',[0 350],'ticklabels','none', 'colorbar', false);
hold on 
title('L2','FontSize',14, 'fontweight','bold')


subplot(4,4,3)
imageschs(DEM_L3,[],'caxis',[0 350],'ticklabels','none', 'colorbar', false);
hold on 
title('L3','FontSize',14, 'fontweight','bold')


subplot(4,4,4)
imageschs(DEM_L4,[],'caxis',[0 350],'colorbarylabel','Elevation [m]','ticklabels','none');
hold on 
title('L4','FontSize',14, 'fontweight','bold')


subplot(4,4,5)
imageschs(G_L1,[],'caxis',[0 50],'ticklabels','none', 'colorbar', false);
hold on
ylabel('Gradient','FontSize',14, 'fontweight','bold')

subplot(4,4,6)
imageschs(G_L2,[],'caxis',[0 50],'ticklabels','none', 'colorbar', false);


subplot(4,4,7)
imageschs(G_L3,[],'caxis',[0 50],'ticklabels','none', 'colorbar', false);


subplot(4,4,8)
imageschs(G_L4,[],'caxis',[0 50],'colorbarylabel','Gradient [deg]','ticklabels','none');


subplot(4,4,9)
imageschs(A_L1,[],'caxis',[0 360],'ticklabels','none', 'colorbar', false);
hold on
ylabel('Aspect','FontSize',14, 'fontweight','bold')


subplot(4,4,10)
imageschs(A_L2,[],'caxis',[0 360],'ticklabels','none', 'colorbar', false);


subplot(4,4,11)
imageschs(A_L3,[],'caxis',[0 360],'ticklabels','none', 'colorbar', false);


subplot(4,4,12)
imageschs(A_L4,[],'caxis',[0 360],'colorbarylabel','Aspect [deg]','ticklabels','none');


subplot(4,4,13)
imageschs(C_L1,[],'caxis',[-0.1 0.1],'ticklabels','none', 'colorbar', false);
hold on
ylabel('Curvature','FontSize',14, 'fontweight','bold')


subplot(4,4,14)
imageschs(C_L2,[],'caxis',[-0.1 0.1],'ticklabels','none', 'colorbar', false);

 
subplot(4,4,15)
imageschs(C_L3,[],'caxis',[-0.1 0.1],'ticklabels','none', 'colorbar', false);


subplot(4,4,16)
imageschs(C_L4,[],'caxis',[-0.1 0.1],'colorbarylabel','Curvature [-]','ticklabels','none');



sgtitle('Elevation, Gradient, Aspect and Curvature of Lithologies L1, L2, L3 and L4','FontSize',14, 'fontweight','bold')


%% Question 2: Generate one figure with 4 subplots. In the four subplots, you
% plot the distributions of (1) elevation, (2) slope, (3) aspect and
% 4(curvature) for all four lithologies, respectively. (40 points)

%convert to vector to plot histogram
DEM_L1_vector = DEM_L1.Z(:); 
G_L1_vector = G_L1.Z(:);
A_L1_vector = A_L1.Z(:); 
C_L1_vector = C_L1.Z(:);

DEM_L2_vector = DEM_L2.Z(:); 
G_L2_vector = G_L2.Z(:);
A_L2_vector = A_L2.Z(:); 
C_L2_vector = C_L2.Z(:);

DEM_L3_vector = DEM_L3.Z(:); 
G_L3_vector = G_L3.Z(:);
A_L3_vector = A_L3.Z(:); 
C_L3_vector = C_L3.Z(:);

DEM_L4_vector = DEM_L4.Z(:); 
G_L4_vector = G_L4.Z(:);
A_L4_vector = A_L4.Z(:); 
C_L4_vector = C_L4.Z(:);

figure
subplot(4,4,1)
histfit(DEM_L1_vector,90);
hold on 
xlim([0 400]);
xlabel('Elevation [m]');
ylabel('Frequency');
title('L1','FontSize',14, 'fontweight','bold')


subplot(4,4,2)
histfit(DEM_L2_vector,90);
hold on 
xlim([0 400]);
xlabel('Elevation [m]');
title('L2','FontSize',14, 'fontweight','bold')


subplot(4,4,3)
histfit(DEM_L3_vector,90);
hold on 
xlim([0 400]);
xlabel('Elevation [m]');
title('L3','FontSize',14, 'fontweight','bold')


subplot(4,4,4)
histfit(DEM_L4_vector,90);
hold on 
xlim([0 400]);
xlabel('Elevation [m]');
title('L4','FontSize',14, 'fontweight','bold')


subplot(4,4,5)
histfit(G_L1_vector,90);
hold on 
xlim([0 80]);
xlabel('Gradient [deg]');
ylabel('Frequency');


subplot(4,4,6)
histfit(G_L2_vector,90);
hold on 
xlim([0 80]);
xlabel('Gradient [deg]');

subplot(4,4,7)
histfit(G_L3_vector,90);
hold on 
xlim([0 80]);
xlabel('Gradient [deg]');

subplot(4,4,8)
histfit(G_L4_vector,90);
hold on 
xlim([0 80]);
xlabel('Gradient [deg]');

subplot(4,4,9)
histfit(A_L1_vector,90);
hold on 
xlim([0 350]);
xlabel('Aspect [deg]');
ylabel('Frequency');

subplot(4,4,10)
histfit(A_L2_vector,90);
hold on 
xlim([0 350]);
xlabel('Aspect [deg]');

subplot(4,4,11)
histfit(A_L3_vector,90);
hold on 
xlim([0 350]);
xlabel('Aspect [deg]');

subplot(4,4,12)
histfit(A_L4_vector,90);
hold on 
xlim([0 350]);
xlabel('Aspect [deg]');

subplot(4,4,13)
histfit(C_L1_vector,90);
hold on 
xlim([-0.5 0.5]);
xlabel('Curvature [-]');
ylabel('Frequency');


subplot(4,4,14)
histfit(C_L2_vector,90);
hold on 
xlim([-0.5 0.5]);
xlabel('Curvature [-]');

subplot(4,4,15)
histfit(C_L3_vector,90);
hold on 
xlim([-0.5 0.5]);
xlabel('Curvature [-]');


subplot(4,4,16)
histfit(C_L4_vector,90);
hold on 
xlim([-0.5 0.5]);
xlabel('Curvature [-]');



sgtitle('Elevation, Slope, Aspect and Curvature Distribution (Bin Size=90) of Lithologies L1, L2, L3 and L4','FontSize',14, 'fontweight','bold')
