
%% Import the Datasets 

Trade_data =readtable('C:\Users\alimu\Desktop\Economics PHD\All_semesters\5th Semester\International Economics\Trade_data.csv','TextType','string');
Dist_data = readtable('C:\Users\alimu\Desktop\Economics PHD\All_semesters\5th Semester\International Economics\dist_cepii.csv','TextType','string');
GDP_data = readtable('C:\Users\alimu\Desktop\Economics PHD\All_semesters\5th Semester\International Economics\GDP_data.csv','TextType','string') ;



%% Cleaning Trade_data
Trade_data = removevars(Trade_data,{'FLOW','FREQUENCY','TIME','UnitCode',...
    'PowerCodeCode','ReferencePeriodCode','ReferencePeriod','Flags','FlagCodes','Unit','PowerCode',});

Trade_data =Trade_data(Trade_data.Flow =='Exports',:);
Trade_data =Trade_data(Trade_data.Frequency =='Annual',:);

countries_to_select = {'AUS', 'AUT', 'BEL', 'BRA', 'CAN', 'CHE', 'CHL', 'CHN', 'DEU', 'FIN', 'FRA', 'GBR'};

logical_index = ismember(Trade_data.LOCATION,countries_to_select);

% Use the logical index to select the rows or specific countries  from the table
Trade_data = Trade_data(logical_index, :);
%Trade_data = removevars(Trade_data,"Flow")



% Clean GDP data
GDP_data = removevars(GDP_data,{'TRANSACT','MEASURE','TIME',...
    'UnitCode','PowerCodeCode','ReferencePeriodCode','FlagCodes','Flags','Measure','ReferencePeriod','Transaction','PowerCode','Unit'});

GDP_data;
% Use the logical index to select the rows or specific countries  from the table
logical_index_gdp = ismember(GDP_data.LOCATION,countries_to_select);
GDP_data = GDP_data(logical_index_gdp,:);




% Clean Distance data 
Dist_data = removevars(Dist_data,{'contig','comlang_off','distw',...
    'colony','comcol','curcol','col45','smctry','distwces','distcap','comlang_ethno'});

Dist_data;
% Use the logical index to select the rows or specific countries  from the table
logical_index_dist = ismember(Dist_data.iso_o,countries_to_select);
Dist_data = Dist_data(logical_index_dist,:);

% missing_dist = ismissing(Dist_data.iso_o)
% missing_dist_2 = ismissing(Dist_data.iso_d)
% 
% Dist_data(missing_dist, :) = [];
% Dist_data(missing_dist_2, :) = [];

%%
% Panel Data Formulation



% % Join tables
GDP_trade_Data = outerjoin(GDP_data,Trade_data,"LeftKeys",["LOCATION","Year"],...
    "RightKeys",["LOCATION","Time"]);



% missing_values = ismissing(GDP_trade_Data.LOCATION_GDP_data)
% missing_values_1 = ismissing(GDP_trade_Data.PARTNER)
% 
% GDP_trade_Data(missing_values, :) = [];
% GDP_trade_Data(missing_values_1, :) = [];
% 
% ismissing(GDP_trade_Data.LOCATION_GDP_data) 



% Join tables
Merged_Data = innerjoin(GDP_trade_Data,Dist_data,"LeftKeys",...
    ["LOCATION_Trade_data","PARTNER"],"RightKeys",["iso_o","iso_d"]);