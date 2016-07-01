function s=isempty(model)
s=isempty(model.mean) && isempty(model.variance);