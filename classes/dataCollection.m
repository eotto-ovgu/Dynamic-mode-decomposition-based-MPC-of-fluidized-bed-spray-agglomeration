classdef dataCollection
    properties
        datasets cell
    end
    methods
        function obj = dataCollection()
        end
        function obj = addData(obj,varargin)
            for i = 1:numel(varargin)
                obj.datasets{end+1} = varargin{i};
            end
        end
        function plotData(obj,index)
            fignum = randi(10000);
            figure(fignum)
            t = tiledlayout(6,2);
            
            for i = index
                dat = obj.datasets{i};
                figure(i+9687)
                surf(dat.p.xgrid,dat.tspan/3600,dat.x','EdgeColor','none')
                xlabel('d / mm'); ylabel('time / h'); zlabel('n_3 / mm^{-1}')
                xlim([obj.datasets{i}.p.xgrid(1),1.5]); title(dat.label)

                figure(fignum)
                nexttile(1,[3,1]); hold on
                plot(dat.tspan/3600,sauter(dat.p.xgrid,dat.p.dx,dat.x))
                xlabel('time / h'); ylabel('d_{32} / h')
                nexttile(7,[3,1]); hold on
                plot(dat.tspan(1:end-1)/3600,dat.u)
                xlabel('time / h'); ylabel('u / h')
                nexttile(8,[3,1]); hold on
                plot(dat.tspan/3600,dat.p.rho_p*pi/6*Moments(dat.x,dat.p.xgrid,dat.p.dx,3))
                xlabel('time / h'); ylabel('m / h')
                
                nexttile(2,[1 1]); hold on
                plot(dat.tspan/3600,Moments(dat.x,dat.p.xgrid,dat.p.dx,0))
                nexttile(4,[1 1]); hold on
                plot(dat.tspan/3600,Moments(dat.x,dat.p.xgrid,dat.p.dx,1))
                nexttile(6,[1 1]); hold on
                plot(dat.tspan/3600,Moments(dat.x,dat.p.xgrid,dat.p.dx,2))
                xlabel('time / h');
            end
        end
        function exportData(obj,index)
            for i = index
                dat = obj.datasets{i};
                d32 = sauter(dat.p.xgrid,dat.p.dx,dat.x);
                m = bedmass(dat.p.rho_p,dat.p.xgrid,dat.p.dx,dat.x);
                writematrix([dat.tspan'/3600 d32'],['FilesForAgglo\results\data\d32',dat.label,'.csv'])
                writematrix([dat.tspan(1:end-1)'/3600 dat.u'],['FilesForAgglo\results\data\u',dat.label,'.csv'])
                writematrix([dat.tspan'/3600 m'],['FilesForAgglo\results\data\m',dat.label,'.csv'])

                % PSD
                n = dat.p.xgrid'.^(-3).*dat.x.*dat.p.dx'./dat.p.dv';
                % reduction for pgflpots
                x = dat.p.xgrid(2:2:end);
                t = dat.tspan(1:4:end);
                n3 = dat.x(2:2:end,1:4:end);
                
                v = dat.p.v(2:2:end);
                n = n(2:2:end,1:4:end);

                figure()
                surf(v,t,n')
                h = gca
                set(h,'xscale','log')

                [T,V] = meshgrid(t,v);
                [T,X] = meshgrid(t,x);

                writematrix([X(:), T(:)/3600, n3(:)],['FilesForAgglo\results\data\n3',dat.label,'.csv'])
                writematrix([V(:), T(:)/3600, n(:)],['FilesForAgglo\results\data\n',dat.label,'.csv'])
            end
        end
    end
end