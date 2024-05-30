function [Xh,idx] = patch_pruning(Xh,threshold)

pvars = var(Xh, 0, 1);

idx = pvars > threshold;

Xh = Xh(:, idx);