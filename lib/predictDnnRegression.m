%%
% predict objective variables by DNN Regression
% input:
%  mdl          trained DNN Regression mld
%  X            regression explanatory variable (data point x expNum)

function Y = predictDnnRegression(mdl,X)
    Y = predict(mdl.network, X.');
    Y = Y.';
end
