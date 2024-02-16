function y_model = weightedExpCosSumModel(params, lag )
    numTerms = floor(numel(params) / 3);
    weights = params(1:numTerms);
    decayConstants = params(numTerms + 1:2*numTerms);
    angularFrequencies = params(2*numTerms + 1:end);
    expCosTerms = zeros(size(xData));
    for i = 1:numTerms
        expCosTerms = expCosTerms + weights(i) * exp(-decayConstants(i) * xData) .* cos(angularFrequencies(i) * xData);
    end
    y_model = expCosTerms;
end