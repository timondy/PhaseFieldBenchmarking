function postprocess04(outpath)

out = [];

figure;
hold on;

cmp = [];
lnt1 = linspace(-5,7,1000);
lnt2 = linspace(-5,2,1000);

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
        
        info = fullfile(outpath,var{a},n{c},'energy.dat');
        if ~exist(info, 'file')
            continue
        end
        
        dat = load(info);
        lnt = log(dat(:,1));
        lnE = log(dat(:,2));
        
        tmp.lnE1 = interp1(lnt,lnE,lnt1,'linear');
        tmp.lnE2 = interp1(lnt,lnE,lnt2,'linear');
        cmp = [cmp tmp];
        if c>1
            dE1 = abs(cmp(c).lnE1-cmp(c-1).lnE1);
            dE2 = abs(cmp(c).lnE2-cmp(c-1).lnE2);
            entry.D1 = sum((dE1(2:end)+dE1(1:end-1)).*(lnt1(2:end)-lnt1(1:end-1)))/2;
            entry.D2 = sum((dE2(2:end)+dE2(1:end-1)).*(lnt1(2:end)-lnt1(1:end-1)))/2;
            plot (lnt,lnE,'DisplayName',sprintf('N=%d',entry.n))
        else
            entry.D1 = NaN;
            entry.D2 = NaN;
            plot (lnt,lnE);
            legend (sprintf('N=%d',entry.n))
        end
        
        disp(sprintf('%s: %f %f', info, entry.D1, entry.D2));
        out = [out entry];
    end
end

writetable(struct2table(out), 'benchmarkIV_sge.csv')

end
