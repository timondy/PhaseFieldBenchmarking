function postprocess04(outpath,benchmark)

if nargin<2
    benchmark='run9.mat';
end

out = [];

figure;
hold on;

load(benchmark)
lt2 = log(t(2:end));
lE2 = log(E(2:end));

plot (lt2,lE2, 'k-');
legend ('benchmark')

theta1 = linspace(-5,7,1000); % Use for benchmark D1
theta2 = linspace(-5,2,1000); % Use for benchmark D2

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
        lt1 = log(dat(2:end,1));
        lE1 = log(dat(2:end,2));

        % Linearly interpolate values to a uniform grid
        lE1interp = interp1(lt1,lE1,theta1);
        lE2interp = interp1(lt2,lE2,theta1);
        entry.D1 = trapz(theta1, abs(lE1interp-lE2interp));

        lE1interp = interp1(lt1,lE1,theta2);
        lE2interp = interp1(lt2,lE2,theta2);
        entry.D2 = trapz(theta2, abs(lE1interp-lE2interp));

        plot (lt1,lE1,'DisplayName',sprintf('N=%d',entry.n))

        fprintf('%s: %f %f\n', info, entry.D1, entry.D2);
        out = [out entry];
    end
end

writetable(struct2table(out), 'benchmark04.csv')

xlabel('ln t')
ylabel('ln E(t)')
xlim([-8,8])
saveas(gcf,'benchmark04.eps','epsc')

end
