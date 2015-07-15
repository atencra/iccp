function iccpairs_ccc_pos_neg_summary(posdata, negdata, posnegdata, label)


fprintf('\n');
pos = simple_stats(posdata);
neg = simple_stats(negdata);
posneg = simple_stats(posnegdata);
fprintf('\n');
fprintf('%s: Pos: mean = %.3f\n', label, pos.mn);
fprintf('%s: Pos: sd = %.3f\n', label, pos.sd);
fprintf('%s: Pos: se = %.3f\n', label, pos.se);
fprintf('%s: Pos: median = %.3f\n', label, pos.md);
fprintf('%s: Pos: mad = %.3f\n', label, pos.mad);
fprintf('\n');
fprintf('%s: Neg: mean = %.3f\n', label, neg.mn);
fprintf('%s: Neg: sd = %.3f\n', label, neg.sd);
fprintf('%s: Neg: se = %.3f\n', label, neg.se);
fprintf('%s: Neg: median = %.3f\n', label, neg.md);
fprintf('%s: Neg: mad = %.3f\n', label, neg.mad);
fprintf('\n');
fprintf('%s: Pos+Neg: mean = %.3f\n', label, posneg.mn);
fprintf('%s: Pos+Neg: sd = %.3f\n', label, posneg.sd);
fprintf('%s: Pos+Neg: se = %.3f\n', label, posneg.se);
fprintf('%s: Pos+Neg: median = %.3f\n', label, posneg.md);
fprintf('%s: Pos+Neg: mad = %.3f\n', label, posneg.mad);

[p_pos_vs_neg, h] = ranksum(posdata, negdata);
[p_pos_vs_pos_neg, h] = ranksum(posdata, posnegdata);
[p_neg_vs_pos_neg, h] = ranksum(negdata, posnegdata);
fprintf('\n');
fprintf('%s: Pos vs. Neg: p = %.5f\n', label, p_pos_vs_neg);
fprintf('%s: Pos vs. Pos+Neg: p = %.5f\n', label, p_pos_vs_pos_neg);
fprintf('%s: Neg vs. Pos+Neg: p = %.5f\n', label, p_neg_vs_pos_neg);
fprintf('\n');