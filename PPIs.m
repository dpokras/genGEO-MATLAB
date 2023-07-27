classdef PPIs
    methods(Static)
        function PPI = load_PPI()

            % A persistent variable won't be loaded into memory every time, only the
            % first time. It's scope is specific to the function.
            persistent data;
            persistent headers;
            if (isempty(data))
                [data, headers] = xlsread(fullfile('data', 'PPI_Table.xlsx'),'Sheet1');
  
            end
            PPI.data = data;
            PPI.headers = headers;
        end
        function PPI_value = return_PPI(indexName, costYear, params);

            idxPPI = find(strcmp(params.PPI.headers(1,:), indexName));
            if (isempty(idxPPI))
                throw(MException('PPI:unknownPPIIndex','Cant find cost index PPI'));
            end
            
            idxYear = find(params.PPI.data(:,1) == costYear);
            if (isempty(idxYear))
                throw(MException('PPI:unknownYearIndex','Cant find cost index year'));
            end
            
            PPI_value = params.PPI.data(idxYear, idxPPI);

        end
    end
end
