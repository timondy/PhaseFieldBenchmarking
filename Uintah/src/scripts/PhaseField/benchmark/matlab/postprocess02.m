function postprocess02(outpath)

out = [];

var = split(ls(outpath));
for a=1:length(var)
    if isempty(var{a})
        continue
    end
    entry.var = var{a};

    eps = split(ls(fullfile(outpath,var{a})));
    for b=1:length(eps)
        if isempty(eps{b})
            continue
        end
        val = textscan(eps{b}, 'eps%d');
        entry.eps = double(val{1})/100;

        n = split(ls(fullfile(outpath,var{a},eps{b})));
        for c=1:length(n)
            if isempty(n{c})
                continue
            end
            val = textscan(n{c}, 'n%d');
            entry.n = double(val{1});

            k = split(ls(fullfile(outpath,var{a},eps{b},n{c})));
            for d=1:length(k)
                if isempty(k{d})
                    continue
                end
                val = textscan(k{d}, 'k%f');
                entry.k = val{1};

                info1 = fullfile(outpath,var{a},eps{b},n{c},k{d},'u1.dat');
                info2 = fullfile(outpath,var{a},eps{b},n{c},k{d},'u2.dat');
                if ~exist(info1, 'file') || ~exist(info2, 'file')
                    continue
                end

                entry.T1 = computeT(info1);
                entry.T2 = computeT(info2);

                fprintf('%s: %f %f\n', info1, entry.T1, entry.T2);
                out = [out entry];
            end
        end
    end
end

writetable(struct2table(out), 'benchmark02.csv')

end
