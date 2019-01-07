function postprocess01(outpath)

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
                
                info = fullfile(outpath,var{a},eps{b},n{c},k{d},'u0.dat');
                if ~exist(info, 'file') 
                    continue
                end
                
                entry.T = computeT(info);
                
                disp(sprintf('%s: %f', info, entry.T));
                out = [out entry];
            end
        end
    end
end

writetable(struct2table(out), 'benchmarkI_sge.csv')

end
