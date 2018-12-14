function postprocess03(outpath)

out = [];

var = split(ls(outpath));
for a=1:length(var)
    if isempty(var{a})
        continue
    end
    entry.var = var{a};

    n = split(ls(fullfile(outpath,var{a})));
    for c=1:length(n)
        if isempty(n{c})
            continue
        end
        val = textscan(n{c}, 'n%d');
        entry.n = double(val{1});

        k = split(ls(fullfile(outpath,var{a},n{c})));
        for d=1:length(k)
            if isempty(k{d})
                continue
            end
            val = textscan(k{d}, 'k%f');
            entry.k = val{1};

            info = fullfile(outpath,var{a},n{c},k{d},'u0.dat');
            if ~exist(info, 'file')
                continue
            end

            entry.T = computeT(info);

            disp(sprintf('%s: %f', info, entry.T));
            out = [out entry];
        end
    end
end

writetable(struct2table(out), 'benchmark03.csv')
