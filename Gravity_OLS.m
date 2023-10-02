%% Distance_Data    
            
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 14);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["iso_o", "iso_d", "contig", "comlang_off", "comlang_ethno", "colony", "comcol", "curcol", "col45", "smctry", "dist", "distcap", "distw", "distwces"];
opts.VariableTypes = ["categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "iso_d", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["iso_o", "iso_d"], "EmptyFieldRule", "auto");

% Import the data
dist_cepii2 = readtable("C:\Users\alimu\Desktop\Economics PHD\All_semesters\5th Semester\International Economics\/dist_cepii.csv", opts);

% Convert to output type
iso_o2 = dist_cepii2.iso_o;
iso_d2 = dist_cepii2.iso_d;
dist2 = dist_cepii2.dist;


% Clear temporary variables
clear dist_cepii2

% Clear temporary variables
clear opts

% Display results
iso_o2, iso_d2, dist2





% Select Countries

countries = {'USA','NLD','JPN','DEU','GBR','IND','NOR','ITA','CAN','BEL','HUN','BRA','AUS','ESP','POL'};
[~,idxO] = ismember(iso_o2, countries);
[~,idxD] = ismember(iso_d2, countries);

iso_sele = sort(iso_o2(idxO & idxD),1)
iso_d_sel = sort(iso_d2(idxO & idxD),1)
dist_sel = dist2(idxO & idxD)

% reshape(iso_sele,[],15)
% 
% reshape(iso_d_sel,15,[])


% 
matrix_dis = reshape(dist_sel,[],15)
% matrix = matrix'
% % 
% % 
% [~,idx] = sort(countries)
% matrix = matrix(idx,idx)
% % 
% disp(matrix)
% 
% 
% 
% % Create table 
% T = cell2table(num2cell(matrix), 'VariableNames', countries');
% T.Properties.RowNames = countries;
% 
% % Display table 
% disp(T)

%% Export_Data


% Import the Datasets 

Trade_data =readtable('C:\Users\alimu\Desktop\Economics PHD\All_semesters\5th Semester\International Economics\Trade_data.csv','TextType','string');

% Cleaning Trade_data
Trade_data = removevars(Trade_data,{'FLOW','FREQUENCY','TIME','UnitCode',...
    'PowerCodeCode','ReferencePeriodCode','ReferencePeriod','Flags','FlagCodes','Unit','PowerCode'});

    Trade_data =Trade_data(Trade_data.Flow =='Exports',:);
    Trade_data =Trade_data(Trade_data.Frequency =='Annual',:);
    Trade_data =Trade_data(Trade_data.Time == 2020,:);


    Trade_data = removevars(Trade_data,{'Flow','Frequency'})
    Trade_data = removevars(Trade_data,{'PartnerCountry','ReporterCountry','Time'})
    Trade_data = sortrows(Trade_data,"LOCATION")
    Trade_data = renamevars(Trade_data,'LOCATION',"Ori")
    Trade_data = renamevars(Trade_data,'PARTNER',"Dest")
    
   
    %countries = {'USA','NLD','JPN','DEU','GBR','IND','NOR','ITA','CAN','BEL','HUN','BRA','AUS','ESP','POL'};


    Trade_data_r = readtable("Pivot_Expo.csv")
   
    myMatrix_exp = table2array(Trade_data_r)
% % 
% % Example dataMatrix (replace with your actual data)
% dataMatrix = myMatrix_exp;
% % Replace with your actual data
% 
% % Get unique countries 
% countries = unique(dataMatrix(:,1));
% 
% % Extract columns for country1, country2, value
% c1 = dataMatrix(:,1);
% c2 = dataMatrix(:,2); 
% val =dataMatrix(:,3);
% 
% % Create index vectors for reshaping
% idx1 = ismember(c1, countries);
% idx2 = ismember(c2, countries);
% 
% % Reshape and transpose 
% matrix_exp = reshape(val(idx1 & idx2), [], length(countries))'; 
% 
% % Make symmetric
% for i = 1:size(matrix_exp,1)
%     for j = i+1:size(matrix_exp,2)
%         matrix_exp(j,i) = matrix_exp(i,j); 
%     end
% end
% 
% % Display
% disp(matrix_exp)
% 
% 
% 
% size(matrix_exp)

%matrix_exp = cell2mat(matrix_exp)

% 
% % Get unique origin countries
% oriCountries = unique(dataMatrix(:,1)); 
% 
% % Create map from country name to index
% [~,idx] = ismember(oriCountries,oriCountries);
% 
% % Get origin index for each row 
% oriIdx = arrayfun(@(c) find(ismember(oriCountries,c)), dataMatrix(:,1));
% 
% % Get destination index for each row
% desIdx = arrayfun(@(c) find(ismember(oriCountries,c)), dataMatrix(:,2)); 
% 
% % Convert to linear indices
% oriLinIdx = sub2ind([15 15], oriIdx, ones(size(oriIdx)));
% desLinIdx = sub2ind([15 15], ones(size(desIdx)), desIdx);
% 
% % Reshape into 15x15 matrix
% % Reshape into 15x15 matrix 
% vector = reshape(str2double(dataMatrix(:,3)), 225, 1);
% 
% % Finally reshape 225x1 into 15x15
% exportMatrix = reshape(vector, 15, 15);
% 
% % Display reshaped matrix
% disp(exportMatrix)




GDP_Data
GDP_data = readtable('C:\Users\alimu\Desktop\Economics PHD\All_semesters\5th Semester\International Economics\GDP_data.csv','TextType','string') ;



% Cleaning Trade_data


%% GDP_Data
% Clean GDP data
GDP_data = removevars(GDP_data,{'TRANSACT','MEASURE','TIME',...
    'UnitCode','PowerCodeCode','ReferencePeriodCode','FlagCodes','Flags','Measure','ReferencePeriod','Transaction','PowerCode','Unit'});



GDP_data =GDP_data(GDP_data.Year == 2020,:);


countries = {'USA','NLD','JPN','DEU','GBR','IND','NOR','ITA','CAN','BEL','HUN','BRA','AUS','ESP','POL'};

% Use the logical index to select the rows or specific countries  from the table
logical_index_gdp = ismember(GDP_data.LOCATION,countries);
GDP_data = GDP_data(logical_index_gdp,:);

GDP_data = removevars(GDP_data,"Year")

GDP_data = removevars(GDP_data,"Country")



% Converting the Table into Martrix
myMatrix_GDP = table2array(GDP_data)

% Check if it's a matrix
if ismatrix(myMatrix_GDP)
    disp('It is a matrix.');
else
    disp('It is not a matrix.');
end

% Extract specific columns as vectors
vectorColumn1gdp = myMatrix_GDP(:, 1);
vectorColumn2gdp = myMatrix_GDP(:, 2);



vectorColumn2gdp =str2double(vectorColumn2gdp)
vectorColumn2gdp =1000000*vectorColumn2gdp


% Example data (replace with your actual data)
countryNames = vectorColumn1gdp
gdpValues = vectorColumn2gdp; % GDP values in the same order as countryNames



reshape(vectorColumn2gdp,[],15)




% Example cell array (replace with your actual data)
% cellArray = matrix_exp; % Create an empty 15x15 cell array
% 
% double_array = cell2mat(cellfun(@str2double, cellArray, 'UniformOutput', false));






concatenated_matrix = horzcat(matrix_dis,myMatrix_exp,gdpValues)



%% OLS Regression Model 

% Label the dependent and independent varaibles 


% % Extract Varaibles
X = concatenated_matrix(:,16:end-1)
d = concatenated_matrix(:,1:15)
GDP = concatenated_matrix(:,end)

% 
% % Take logs
lnX = log(X);
lnX = lnX(:);
lnGDPi = log(GDP);
lnGDPj=lnGDPi'
lnd = log(d);
lnd  = lnd(:);
% 
% % Construct regression matrix
% N = 15; 
% 
% lnGDP = log(GDP); 
% lnGDP_rep = repmat(lnGDP,N,1); 
% lnGDP_transpose = lnGDP';
% 
% lnd = log(d); 
% lnd = lnd(:); 
% 
% A = [ones(N^2,1) lnGDP_rep lnGDP_transpose lnd];

% N = 15;
% 
% % lnX = log(X);  
% % lnd = log(d); 
% % lnGDP = log(GDP);
% 
% ones_col = ones(N^2, 1); 
% 
% lnGDP_rep = repmat(lnGDPi, N, 1); 
% lnGDP_transpose = lnGDPi' 
% 
% lnGDP_rep_1 = repmat(lnGDP_transpose, 15, 15)
% 
% lnGDP_col = lnGDP_rep(:,1)
% 
% lnd_vec = lnd(:); 
% 
% 
% 
% % Regression Matrix
% A = [ones_col lnGDP_rep lnGDP_col lnd_vec];
% 
% 
% [b, bint, ~, ~, stats] = regress(lnX, A)
% 
% % Calculate residuals
% epsilon_hat = lnX - A * b
% 
% % Calculate fitted values
% X_hat = A * b
% 
% % Calculate mean squared error
% mse = sum(epsilon_hat.^2) / (225 - 4)
% 
% % Calculate standard errors of coefficients
% se_alpha_b = sqrt(mse * diag(inv(X' * X)))
% 
% % Calculate 95% confidence intervals
% % alpha_lower = b - tinv(0.975, 225 - 4) * se_alpha_b;
% % alpha_upper = b + tinv(0.975, 225 - 4) * se_alpha_b;



%%

GDP_i = repmat(lnGDPi, 1, 15);
GDP_j = repmat(lnGDPj, 15, 1);
GDP_j = reshape(GDP_j, [], 1);
GDP_i = reshape(GDP_i, [], 1);


b = fitlm([GDP_i,GDP_j, lnd], lnX, "linear");


