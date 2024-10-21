close all
clear
clc


doplot = false;

%set up pathway to where data for fig is stored
%rootpath = '/Users/scherler/Nextcloud_GFZ/work/Paper_finished/2023_ESD_Emma_Nahuelbuta/mfiles/';

rootpath = 'C:\Users\emmal\Dropbox\GFZ\Corrections Esurf paper\';
fn = '10Be_data_Esurf.xlsx';

T = readtable(fn);

% Units used for calculations: cm, year, gram, atom

% Constants
L_n = 160; %attenuation length by spallation - well known
L_muf = 4320; %attenuation length by fast muons from Braucher 2011
L_mun = 1500; %attenuation length by negative muons from Braucher 2011
rho = 2.6; %density
lambda10 = 4.9975e-07;  % decay constant - radioactive decay of 10be


figure('Position',[109 302 1250 420],'Color',[1 1 1])
for site = 1 : 3

    switch (site)
        case 1 % 'NA'
            IX = [1,3,5,4,2];
            maxE1 = 50; maxE2 = 20;
        case 2 % 'LC'
            IX = [10,9];
            maxE1 = 200; maxE2 = 200;
        case 3 % 'SG'
            IX = [13,15,14,16,18,17];
            maxE1 = 30; maxE2 = 30;
%         case 3 % 'SG JGR'
%             IX = [29,30,31,32];
%             maxE1 = 30; maxE2 = 30;
     
    end


    locname = T.SpecificLocation(IX);

    P = T.ProductionRate_spallation__atoms_g_yr_(IX);  %Production rate of 10Be by spallation
    Pmu_fast = T.ProductionRateFASTMuons_atoms_g_yr_(IX);  %Production rate of 10Be by FAST muons
    Pmu_neg = T.ProductionRateNEGMuons_atoms_g_yr_(IX);  %Production rate of 10Be by NEGATIVE muons

    N = T.ConcentrationInSample_N10_(IX);   %N10 or 10Be concentration of the corestones
    Nsig = T.Uncertainty10BeConc__N10_(IX);   %uncertainty in 10Be concentration of the corestones
    z = T.AverageProtrusion_cm_(IX); %protrusion heights of 6 corestones in cm 
    %zm = T.AverageProtrusion_m_(IX);  %protrusion heights of the 6 corestones in m
    %z = zm.*100; %height in cm for the figure
    %E1 = 0.002; % cm/yr - this is the estimated soil erosion rate we are using.
    %E1_real = T.CorrespondingSoilER(1:6);

    %set up vectors for E1 and E2 - these will be the axes of the final figure

    %the concentration of NA15 (76cm) is so low that, since we are adding the
    %concentrations of E1 and E2 together, you need to set the input erosion
    %rates to up to 50 for E1 and E2 so that added together it produces a
    %concentration low enough to overlap with the concentration of NA15.

    all_E1= linspace(0.000,maxE1/10e3,1000);   %range of possible soil erosion rates (0.003 cm/yr = 30m/Myr)
    all_E2= linspace(0.000,maxE2/10e3,1000);   %range of possible corestone erosion rates in cm/yr

    all_E1_mMyr = all_E1.*10e3; %same thing just in m/Myr, for the purpose of plotting to be consistent with units used elsewhere in paper
    all_E2_mMyr = all_E2.*10e3; %same thing just in m/Myr, for the purpose of plotting to be consistent with units used elsewhere in paper

    %t0 = [1000000000 1000000000 1000000000 1000000000 1000000000 1000000000];
    %t = z./E1_real; % exposure time of corestone in years. E1_real uses the measured soil erosion rates surrounding the corestones
    % t = z./E1;   %exposure time of corestone in years. E1 uses a single set soil erosion rate that you can set above.


    C = cell(1,numel(IX));

    for k = 1 : numel(IX)

        P10spall = P(k);
        P10muf = Pmu_fast(k);
        P10mun = Pmu_neg(k);
        this_z = z(k);
        t1 = 1e7;
        %t2 = t(k);

        %below the middle term (exp(-z.*rho/L_n)) is taken out because it is equal
        %to 1 and you multiply by it so it is unnecessary


        %[P10spall P10muf P10mun this_z N(k)]

%         %OLD
%         N10calc = @(E1, E2) (exp(-this_z./(E1-E2).*lambda10).* ...
%             (P10spall./(E1.*rho./L_n + lambda10)) .* (1-exp(-t1.*(E1.*rho./L_n + lambda10))) + ...
%             (P10muf./(E1.*rho./L_muf + lambda10)) .* (1-exp(-t1.*(E1.*rho./L_muf + lambda10))) + ...
%             (P10mun./(E1.*rho./L_mun + lambda10)) .* (1-exp(-t1.*(E1.*rho/L_mun + lambda10))))+ ...
%             ((P10spall./(E2.*rho./L_n + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_n + lambda10))) + ...
%             (P10muf./(E2.*rho./L_muf + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_muf + lambda10))) + ...
%             (P10mun./(E2.*rho./L_mun + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_mun + lambda10))));

%         % NEW
%         N10calc = @(E1, E2) (exp(-this_z./(E1-E2).*lambda10).* ...
%             (P10spall./(E1.*rho./L_n + lambda10)) .* exp(-(this_z./(E1-E2).*E2).*rho./L_n) + ...
%             (P10muf./(E1.*rho./L_muf + lambda10)) .* exp(-(this_z./(E1-E2).*E2).*rho./L_muf) + ...
%             (P10mun./(E1.*rho./L_mun + lambda10)) .* exp(-(this_z./(E1-E2).*E2).*rho./L_mun) ) + ...
%             ((P10spall.*1.2./(E2.*rho./L_n + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_n + lambda10))) + ...
%             (P10muf./(E2.*rho./L_muf + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_muf + lambda10))) + ...
%             (P10mun./(E2.*rho./L_mun + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_mun + lambda10))));

% NEW 2
        N10calc = @(E1, E2) exp(-this_z./(E1-E2).*lambda10).* ...
          ( (P10spall./(E1.*rho./L_n + lambda10)) .* exp(-(this_z./(E1-E2).*E2).*rho./L_n) + ...
            (P10muf./(E1.*rho./L_muf + lambda10)) .* exp(-(this_z./(E1-E2).*E2).*rho./L_muf) + ...
            (P10mun./(E1.*rho./L_mun + lambda10)) .* exp(-(this_z./(E1-E2).*E2).*rho./L_mun) ) + ...
          ( (P10spall.*1.2./(E2.*rho./L_n + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_n + lambda10))) + ...
            (P10muf./(E2.*rho./L_muf + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_muf + lambda10))) + ...
            (P10mun./(E2.*rho./L_mun + lambda10)) .* (1-exp(-this_z./(E1-E2).*(E2.*rho./L_mun + lambda10))));

        % Calculate expected 10Be concentrations based on all possible soil erosion
        % rates (all_E1) and all possible corestone erosion rates (all_E2)
        M = nan(length(all_E1),length(all_E2));
        A = nan(length(all_E1),length(all_E2));
        B = nan(length(all_E1),length(all_E2));

        for i = 1 : length(all_E1)
            this_E1 = all_E1(i);
            for j = 1 : length(all_E2)
                this_E2 = all_E2(j);
                if this_E1>this_E2
                    this_N = N10calc(this_E1,this_E2);
                    A(i,j) = this_N;
                    M(i,j) = abs(this_N - (N(k)));  %write down the absolute value of the difference between the expected and observed concentrations
                end
            end
        end

        v = M(:);
        minv = min(v);
        ix = find(v<2*Nsig(k));  %tell me for which combination of E1 and E2 are the differences between expected and observed concentration < 10Be uncertainty
        [r,c] = ind2sub(size(M),ix);   %setting up the cells used for the fig


        BW = M<2*Nsig(k);
        B = bwboundaries(BW);

        C{k} = B{1};
        %imagesc(M), hold on
        %plot(B{1}(:,2),B{1}(:,1),'r-')


        if (doplot)
            nexttile
            imagesc(all_E2_mMyr,all_E1_mMyr,log10(M))
            xlabel('E2 (m/Myr)')
            ylabel('E1 (m/Myr)')
            hold on
            hc = colorbar;
            hc.Label.String = '(Expected - Observed) Concentrations';
            hold on
            plot(all_E2_mMyr(c),all_E1_mMyr(r), 'r.')
            title(sprintf('Height = %1.2f m',zm(k)))
            axis square
        end
    end



    cols = {'r','b','m','c','g','y'};
    %cols = {'#0072BD','r','b','m','c','g','y', 'k'};

    subplot(1,3,site)
    clear hp
    legstr = cell(1,numel(IX));
    for i = 1 : length(C)
        hp(i) = patch(all_E2(C{i}(:,2)).*1e4,all_E1(C{i}(:,1)).*1e4,cols{i});
        hp(i).FaceAlpha = 0.5;
        hold on
        legstr{i} = sprintf('%s (%1.0f cm)',locname{i},z(i));
    end
    xlim([0 maxE2])
    ylim([0 maxE1])
    xlabel('Boulder erosion rate, E2 (m/Myr)')
    ylabel('Soil erosion rate, E1 (m/Myr)')
    box on
    axis square
    legend(hp,legstr,'location','northeast');

end


% export_fig 'boulders_2_new.png' '-m2'


% figure
% for k=1:6
%
%     imagesc(B)
%     hold on
%     plot(all_E1_mMyr(c),all_E2_mMyr(r), 'r.')
%
% end

%The idea behind Matrix M:
%Matrix M is a grid of 1 million pixels, which show the differences between
%expected and observed 10Be concentrations in the corestones. The predicted concentrations
%are based on our imposed set of 1000 possible soil erosion rates (E1, ranging from 0 to 40 m/Myr) and 1000
%possible corestone erosion rates (E2, ranging from 0 to 30 m/Myr), which are input into the classic equation that
%calculates 10Be concentrations (here represented by the internal function
%N10calc).

%Instead of single values for E1 and E2, we offer the function N10calc
%vectors of 1000 possibilities for each. a For-loop iterates over each of them to
%create this_N, which contains all of the expected concentrations given the
%full range of imposed possibilities for E1 and E2, but only one at a time - the
%value of this_N is redifined over and over as the for-loop goes through it
%(??? I think!?)

%N10calc also requires other inputs, which are listed with this_N in the
%parentheses and defined elsewehre in the script (N0, this_e, etc).
%each of these inputs have (k) so that they are used for all 6 samples
%(k=6, 6 samples and 6 plots).
%But we also have our measured concentrations and
%we want to see where there is the smallest / acceptable difference between expected
%and measured so we take the absolute values of this_N - N. Matrix M is
%composed of cells that contain the absolute values of the differencs
%between expected and observed. The red block shows the cells in
%Matrix M that are < the measured N10 error (Nsig), aka where expected and observed
%concentrations are closest to each other.

%Matrix A is just for fun and it is just N10calc, so just the
%concentrations, to see what its actually calculating