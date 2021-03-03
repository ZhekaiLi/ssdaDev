function trainIndex = prePredict(allu, ssda_Data)
    % 这只是有一个想法，可以先试试看
    % 或者再之后用聚类来做
    for i=1:length(allu)
        ssda_Data.U - allu(i)
    end
end