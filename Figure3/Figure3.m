%TITLE: MATLAB code for "Motion-robust gradient-echo T2 mapping in the liver using readout alternation"
%Matlab & Python scripts to generate table S1 in the paper
%Author: Daiki Tamada
%Affiliation: Department of Radiology, University of Wisconsin-Madison
%Date: 9/21/2026
%Email: dtamada@wisc.edu

%By downloading, installing, or otherwise accessing or using the Software , you ("Recipient") agree to receive and use the above-
%identified SOFTWARE subject to the following terms, obligations and restrictions. If you do not agree to all of the following terms, 
%obligations and restrictions you are not permitted to download, install,
%execute, access, or use the SOFTWARE:

%1.	Originators of the SOFTWARE.  Provider is willing to license its rights in the SOFTWARE ("Provider's Rights") to academic researchers to use free of charge solely for academic, non-commercial research purposes subject to the terms and conditions outlined herein. The SOFTWARE was created at the University of Wisconsin ("UW") by Dakai Tamada. Please note Provider's Rights may include, but are not limited to, certain patents or patent applications owned by the Wisconsin Alumni Research Foundation ("WARF"). 
%2.	Limited License.  Provider hereby grants to Recipient a non-commercial, non-transferable, royalty-free, non-exclusive license, without the right to sublicense, under Provider's Rights to  download, install, access, execute and use the SOFTWARE solely for academic, non-commercial research purposes. SOFTWARE may not be used, directly or indirectly, to perform services for a fee or for the production or manufacture of products for sale to third parties. The foregoing license does not include any license to third party intellectual property that may be contained in the SOFTWARE; obtaining a license to such rights is Recipient's responsibility. 
%3.	Restrictions on SOFTWARE use and distribution.  Recipient shall not take, authorize, or permit any of the following actions with the SOFTWARE: (1) modify, translate or otherwise create any derivative works; or (2) publicly display (e.g., Internet) or publicly perform (e.g., present at a press conference); or (3) sell, lease, rent or lend; or (4) use it for any commercial purposes whatsoever. Recipient must fully reproduce and not obscure, alter or remove any of the Provider's proprietary notices that appear on the SOFTWARE, including copyright notices or additional license terms included with any the third party software contained in the SOFTWARE. Recipient may not provide any third party with access to the SOFTWARE or use the SOFTWARE on a timeshare or service bureau basis. Recipient represents that it is compliance with all applicable export control provisions and is not prohibited from receiving the SOFTWARE. 
%4.	Reservation of rights.  Provider retains all rights and title in the SOFTWARE, including without limitation all intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future exist in the SOFTWARE, regardless of form or medium. Provider retains ownership and all of Its rights in the SOFTWARE, including all of its intellectual property rights (e.g., patent, copyright and trade secret rights) that may now or in the future cover the SOFTWARE or any uses of the SOFTWARE, regardless of form or medium; title remains with Provider and the SOFTWARE is merely being loaned to Recipient for the specific purposes and under the specific restrictions stated herein. Nothing in this Agreement grants Recipient any additional rights to the SOFTWARE, any right to obtain any updates or new releases of the SOFTWARE, any commercial license for the SOFTWARE, or any other intellectual property owned or licensed by Provider. Provider has no obligation to provide any support, updates, or bug fixes.
%5.	Disclaimer of Warranty. PROVIDER IS PROVIDING THE SOFTWARE TO RECIPIENT ON AN "AS IS" BASIS. PROVIDER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING THE SOFTWARE OR ANY OUTCOME THAT MAY BE OBTAINED BY USING THE SOFTWARE, AND EXPRESSLY DISCLAIMS ALL SUCH WARRANTIES, INCLUDING WITHOUT LIMITATION ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS. PROVIDER MAKES NO REMEDY THAT THE SOFTWARE WILL OPERATE ERROR FREE OR UNINTERRUPTED.
%6.	Limitation of Liability; Indemnity.  TO THE FULLEST EXTENT PERMITTED BY LAW, IN NO EVENT SHALL PROVIDER BE LIABLE TO RECIPIENT FOR ANY LOST PROFITS OR ANY DIRECT, INDIRECT, EXEMPLARY, PUNITIVE, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING FROM THE SOFTWARE OR ITS USE. FURTHERMORE, IN NO EVENT WILL PROVIDER'S LIABILITY TO RECIPIENT EXCEED $100. PROVIDER HAS NO LIABILITY FOR ANY DECISION, ACT OR OMISSION MADE BY RECIPIENT AS A RESULT OF USE OF THE SOFTWARE. To the extent permitted by applicable law, Recipient agrees to indemnify, defend and hold harmless Provider, UW, and the SOFTWARE authors against all claims and expenses, including legal expenses and reasonable attorneys fees, arising from Recipient's use of the SOFTWARE.
%7.	No use of names/trademarks.  Recipient shall not use Provider's name, or the name of any author of the SOFTWARE or that of UW, in any manner without the prior written approval of the entity or person whose name is being used.
%8.	Termination.  Without prejudice to any other rights, Provider may terminate this Agreement if Recipient fails to comply with the terms of this Agreement for any reason. Upon termination for any reason, Recipient must immediately destroy all copies of the SOFTWARE in Recipient's possession, custody, or control.	
% clc
% clear variables
%De Bazelaire, Cedric MJ, et al. "MR imaging relaxation times of abdominal and pelvic tissues measured in vivo at 3.0 T: preliminary results." Radiology 230.3 (2004): 652-659.
T1      = 809;  % ms
addpath('../tools/')
mov_int = 0.0:0.025:1.0;

%Load simulation results from Julia code
load('simulation_results.mat');

%%
%Build dictionary
TR_DESS      = 6.8;   % ms
TR_PB      = 6.3;   % ms
flip_DESS = 35;
flip_PB = 20;
dphi = 2;

LUTs = build_LUTs_DESS_Bloch(0, flip_DESS, TR_DESS, T1, 1);
LUTs_ROA = build_LUTs_DESS_Bloch(0, flip_DESS, TR_DESS, T1, -1);

LUTs_PB = build_LUTs_DESS_Bloch(dphi, flip_PB, TR_PB, T1, 1);
LUTs_PB_ROA = build_LUTs_DESS_Bloch(dphi, flip_PB, TR_PB, T1, -1);


%%
%Motion simulation with variable significancy

T2map_DESS = T2map_phase(SR, LUTs, 'magnitude');
T2map_DESS_ROA = T2map_phase(SR_ROA, LUTs_ROA, 'magnitude');
T2_mean_DESS = mean(T2map_DESS,2);
T2_mean_DESS_ROA = mean(T2map_DESS_ROA,2);
T2_std_DESS = std(T2map_DESS,[],2);
T2_std_DESS_ROA = std(T2map_DESS_ROA,[],2);

T2map_PB = T2map_phase(PH, LUTs_PB, 'phase');
T2map_PB_ROA = T2map_phase(PH_ROA, LUTs_PB_ROA, 'phase');
T2_mean_PB = mean(T2map_PB,2);
T2_mean_PB_ROA = mean(T2map_PB_ROA,2);
T2_std_PB = std(T2map_PB,[],2);
T2_std_PB_ROA = std(T2map_PB_ROA,[],2);


%%

% Assuming mov_int, R2_mean_OPT, R2_std_OPT, etc. are already in your workspace

figure('Color', 'w');
hold on; grid on;

% Define colors (Using professional color palettes)
blueCol = [0, 0.447, 0.741];
redCol  = [0.850, 0.325, 0.098];

R2_mean_OPT = T2_mean_DESS(:,1)'; R2_std_OPT = T2_std_DESS(:,1)';
R2_mean_ROA = T2_mean_DESS_ROA(:,1)'; R2_std_ROA = T2_std_DESS_ROA(:,1)';

% 1. Plot ROA Shading (Red)
fill([mov_int, fliplr(mov_int)], ...
     [(R2_mean_ROA + R2_std_ROA), fliplr(R2_mean_ROA - R2_std_ROA)], ...
     redCol, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% 2. Plot OPT Shading (Blue)
fill([mov_int, fliplr(mov_int)], ...
     [(R2_mean_OPT + R2_std_OPT), fliplr(R2_mean_OPT - R2_std_OPT)], ...
     blueCol, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% 3. Plot Mean Lines (Plotted last to stay on top)
p1 = plot(mov_int, R2_mean_OPT, '-o', 'Color', blueCol, 'LineWidth', 2, ...
    'MarkerSize', 5, 'MarkerFaceColor', blueCol);
p2 = plot(mov_int, R2_mean_ROA, '-s', 'Color', redCol, 'LineWidth', 2, ...
    'MarkerSize', 5, 'MarkerFaceColor', redCol);

% --- Final Polish ---
xlabel('Peak Velocity (cm/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('T2 (ms)', 'FontSize', 12, 'FontWeight', 'bold');
%title('Liver R_2 estimation vs. Motion Intensity', 'FontSize', 14);

legend([p1, p2], {'DESS without ROA', 'DESS with ROA'}, 'Location', 'NorthWest', 'Box', 'off');
set(gca, 'XLim', [min(mov_int) max(mov_int)], 'FontSize', 11);
set(gca, 'fontname', 'Arial', 'FontSize',16,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);axis square;



%%

% Assuming mov_int, R2_mean_OPT, R2_std_OPT, etc. are already in your workspace

figure('Color', 'w');
hold on; grid on;

% Define colors (Using professional color palettes)
blueCol = [0, 0.447, 0.741];
redCol  = [0.850, 0.325, 0.098];

R2_mean_OPT = T2_mean_PB(:,1)'; R2_std_OPT = T2_std_PB(:,1)';
R2_mean_ROA = T2_mean_PB_ROA(:,1)'; R2_std_ROA = T2_std_PB_ROA(:,1)';

% 1. Plot ROA Shading (Red)
fill([mov_int, fliplr(mov_int)], ...
     [(R2_mean_ROA + R2_std_ROA), fliplr(R2_mean_ROA - R2_std_ROA)], ...
     redCol, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% 2. Plot OPT Shading (Blue)
fill([mov_int, fliplr(mov_int)], ...
     [(R2_mean_OPT + R2_std_OPT), fliplr(R2_mean_OPT - R2_std_OPT)], ...
     blueCol, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% 3. Plot Mean Lines (Plotted last to stay on top)
p1 = plot(mov_int, R2_mean_OPT, '-o', 'Color', blueCol, 'LineWidth', 2, ...
    'MarkerSize', 5, 'MarkerFaceColor', blueCol);
p2 = plot(mov_int, R2_mean_ROA, '-s', 'Color', redCol, 'LineWidth', 2, ...
    'MarkerSize', 5, 'MarkerFaceColor', redCol);

% --- Final Polish ---
xlabel('Peak Velocity (cm/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('T2 (ms)', 'FontSize', 12, 'FontWeight', 'bold');
%title('Liver R_2 estimation vs. Motion Intensity', 'FontSize', 14);

legend([p1, p2], {'PBT2 without ROA', 'PBT2 with ROA'}, 'Location', 'NorthWest', 'Box', 'off');
set(gca, 'XLim', [min(mov_int) max(mov_int)], 'FontSize', 11);
set(gca, 'fontname', 'Arial', 'FontSize',16,'FontWeight','normal','LineWidth',2);
set(gcf,'units','inches','position',[0,0,5,5]);axis square;