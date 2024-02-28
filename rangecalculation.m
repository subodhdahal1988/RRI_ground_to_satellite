data=struct('lat',[],'lon',[],'elev',[],'azm',[],'alt',[],'slant_range',[],'GSP',[]);
wgs84 = wgs84Ellipsoid;
for i=1:length(satdata.geo_lat)
   [az,elev,slantRange]=geodetic2aer(satdata.geo_lat(i),satdata.geo_lon(i),satdata.alt(i),52.16,-106.53,494,wgs84);
   [arclen,az1] = distance(52.16,-106.53,satdata.geo_lat(i),satdata.geo_lon(i));
   GSD=deg2km(arclen);
   data(i).slant_range=slantRange;
   data(i).elev=elev;
   data(i).GSP=GSD;
   data(i).azm=az;
   data(i).lat=satdata.geo_lat(i);
   data(i).lon=satdata.geo_lon(i);
   data(i).alt=satdata.alt(i);
 
end
% Convert structure array to table
dataTable = struct2table(data);

% Specify the file path for saving the CSV file
csvFilePath = '/volumes/subodh/codingphd/slantrange_matlab/20150401.csv';

% Write the table to a CSV file with headers
writetable(dataTable, csvFilePath);
%% 


% Convert structure array to table
dataTable = struct2table(data);

% Specify the file path for saving the CSV file
csvFilePath = '/volumes/subodh/codingphd/slantrange_matlab/20150401.csv';

% Write the table to a CSV file with headers
writetable(dataTable, csvFilePath);

%%
H_values = [110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330];
time_1 = datetime(time1.time, 'InputFormat', 'yyyy MM dd HH mm ss.SSS');
figure('Position', [100, 100, 1000, 1500]);

hold on;

% Store legend labels in a cell array
legend_labels = cell(1, length(H_values));
for i = 1:length(H_values)
    col_name = sprintf('VarName%d', i + 1);
    plot(time_1, frszn.(col_name), 'LineWidth', 7);
    
    % Store legend label for each color
    legend_labels{i} = ['H = ' num2str(H_values(i)) ' km'];
end

hold off;

grid on;
xtickangle(45);
xlabel('Time (UT)', 'FontSize', 30); % Adjusted font size
ylabel('First Fresnel Radius (m)', 'FontSize', 30); % Adjusted font size

% Set the x-axis limits
xlim([time_1(1), time_1(end)]);

% Set tick label font size
set(gca, 'FontSize', 30); % Adjusted font size


% Add legend at the top of the plot
legend(legend_labels, 'Location', 'eastoutside', 'Orientation', 'vertical', 'FontSize', 40);




%%
H_values = [110, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310, 330];
time_1 = datetime(time1.time, 'InputFormat', 'yyyy MM dd HH mm ss.SSS');
figure('Position', [100, 100, 1000, 1500]);

hold on;

% Use Viridis-like colormap
viridis_colors = viridis(length(H_values));

% Store legend labels in a cell array
legend_labels = cell(1, length(H_values));
for i = 1:length(H_values)
    col_name = sprintf('VarName%d', i + 1);
    plot(time_1, frszn.(col_name), 'LineWidth', 3, 'Color', viridis_colors(i, :));
    
    % Store legend label for each color
    legend_labels{i} = ['H = ' num2str(H_values(i)) ' km'];
end

hold off;

grid on;
gridLineWidth = 1.5;  % Adjust as needed
set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.8)

xtickangle(45);
xlabel('Time (UT)', 'FontSize', 45); % Adjusted font size
ylabel('First Fresnel Radius (m)', 'FontSize', 45); % Adjusted font size

% Set the x-axis limits
xlim([time_1(1), time_1(end)]);
ylim([200,1400]);
% Set tick label font size
set(gca, 'FontSize', 30); % Adjusted font size
%set title 
title("First Fresnel Radius at different altitude of irregularity layer",'FontSize', 40)

% Add legend at the top of the plot
legend(legend_labels, 'Location', 'eastoutside', 'Orientation', 'vertical', 'FontSize', 40);

%%
 col_name = sprintf('VarName%d', 8 + 1);
 plot(time_1, frszn.(col_name))
 xlim([time_1(1), time_1(end)]);
 %%
