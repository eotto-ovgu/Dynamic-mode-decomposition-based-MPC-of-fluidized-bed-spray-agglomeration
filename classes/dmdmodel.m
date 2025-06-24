classdef dmdmodel
    properties
        A 
        B
        U
        stateScale double
        inputScale double
        N_x int32 
        N_u int32
        N_del int32 = 0
        N_x_transformed int32
        ext logical = false
        xop                     % operating point
        uop                     % operating point
        stateSpace              % matlab state space object for application in control
        p                       % parameters 
        ts                      % sample interval length
    end
    methods
        function obj = dmdmodel()
            % empty constructor
        end
        function obj = train(obj,dat,N_del,alpha1,alpha2,ext)
            trainState = dat.x;
            trainInput = dat.u;

            % 0) get model/data parameters
            obj.N_x = size(trainState,1);
            obj.N_u = size(trainInput,1);
            obj.N_del = N_del;
            obj.ext = ext;
            obj.xop = dat.xss;
            obj.uop = dat.uss;
            obj.ts = dat.ts;

            obj.p.xgrid = dat.p.xgrid;
            obj.p.dx = dat.p.dx;
            obj.p.v = dat.p.v;
            obj.p.dv = dat.p.dv;
            obj.p.rho_p = dat.p.rho_p;
            
            % 1) coordinate transformations (deviation variables, scaling, 
            % time delay embedding, nonlinear measurements)
            [X,Xp,Gam,obj.N_x_transformed,obj.stateScale,obj.inputScale] = CoordinateTransformation(obj,trainState,trainInput,N_del,true);
            
            % 4) train DMD model
            [At, Bt, Uh, ~] = DMDc(X,Xp,Gam,obj.N_x_transformed,obj.N_u,alpha1,alpha2);

            % 5) save DMD model
            obj.A = At;
            obj.B = Bt;
            obj.U = Uh;

            % 6) save state space object for usage in simulink
            obj.stateSpace = createStateSpaceObject(obj);
        end
        function [obj, predDat] = test(obj,testdat)
            % extract data
            testState = testdat.x;
            testInput = testdat.u;
            tspan = testdat.tspan;

            % define prediction horizon (default = test data length)
            N_t = size(testState,2) - 1;

            % 1) coordinate transformations (time delay embedding/
            % nonlinear measurements)
            [X,~,Gam,~,~,~] = CoordinateTransformation(obj,testState,testInput,obj.N_del,false);

            % 2) compute prediction with DMD model
            x = X(:,1);
            xt = obj.U'*x;

            for i = 1:N_t-obj.N_del
                xt(:,i+1) = obj.A*xt(:,i) + obj.B*Gam(:,i);
                x(:,i+1) = obj.U*xt(:,i+1);
            end
            
            % 3) compute inverse coordinate transfomration back into original
            % space
            predState = inverseCoordinateTransformation(obj.N_x,obj.N_x_transformed,...
                            obj.N_del,obj.xop,obj.stateScale,obj.inputScale,obj.ext,x);

            % 7) gather and plot prediction data
            predDat = data();
            predDat = predDat.importData(obj.xop,obj.uop,predState,testInput,obj.p,tspan,obj.ts,[testdat.label,'_pred']);
            col = dataCollection();
            col = col.addData(testdat,predDat);
            col.plotData([1 2])
        end
        function stateSpace = createStateSpaceObject(obj)  
            if obj.N_del == 0
                temp = [];
                temp1 = [];
            elseif obj.N_del == 1
                temp = diag(ones(1,obj.N_del-1),1);
                temp1 = 0;
            else
                temp = diag(ones(1,obj.N_del-1),1);
                temp1 = [zeros(obj.N_del-1,1); 1];
            end


            Ac = [obj.A, obj.B(:,1:obj.N_del); ...
                 zeros(obj.N_del,size(obj.A,2)), temp];
            Bc = [obj.B(:,obj.N_del+1); temp1];

            % measurement function
            g = @(x,u) sauter(obj.p.xgrid,obj.p.dx,computeDistribution(obj,x(1:size(obj.A,1)))) - sauter(obj.p.xgrid,obj.p.dx,obj.xop);
    
            % compute jacobians
            Cc = jacobianest(@(x) g(x(1:size(obj.A,1)),obj.uop), zeros(size(Ac,1),1));
            Dc = jacobianest(@(u) g(zeros(size(obj.A,1)+obj.N_del,1),u), 0);
            
%             Cc = eye(size(Ac,1));
%             Dc = zeros(size(Ac,1),1);

            % compute state space
            stateSpace = ss(Ac,Bc,Cc,Dc,obj.ts);
        end
    end
end

function [X,Xp,Gam,N_x_transformed,stateScale,inputScale] = CoordinateTransformation(obj,state,input,N_del,trainflag)
    N_x = size(state,1);
    N_u = size(input,1);
    N_t = size(state,2) - 1;

    % 1) compute deviation variable
    state = state - obj.xop;
    input = input - obj.uop;

    % 2) scaling 
    [state, input, stateScale, inputScale] = scaleData(obj,state,input,trainflag);

    % 3) nonlinear transformations
    if obj.ext == true
        N_x_transformed = (obj.N_x + obj.N_x^2)*(N_del+1);
        stateNew = zeros(N_x+N_x^2,N_t+1); 
        for i = 1:N_t+1 
            temp = state(:,i)*state(:,i).';
            stateNew(:,i) = [state(:,i); reshape(temp,numel(temp),1)];
        end
    else
        N_x_transformed = obj.N_x*(N_del+1);
        stateNew = state;
    end

    % 4) time delay embedding (N_del = 0 => no time delays)
    X = zeros(size(stateNew,1)*(N_del+1),N_t-N_del);
    Xp = zeros(size(X));
    Gam = zeros(N_u*(N_del+1),N_t-N_del);

    for i = 1:N_t-N_del
        t1 = stateNew(:,i:i+N_del);
        t2 = stateNew(:,i+1:i+N_del+1);
        t3 = input(:,i:i+N_del);

        X(:,i) = t1(:);
        Xp(:,i) = t2(:);
        Gam(:,i) = t3(:);
    end
end

function [state, input, stateScale, inputScale] = scaleData(obj,state,input,trainflag)
    if trainflag == true
        stateScale = max(max(state));
        inputScale = max(max(input));

        state = state/stateScale;
        input = input/inputScale;
    else
        state = state/obj.stateScale;
        input = input/obj.inputScale;

        stateScale = [];
        inputScale = [];
    end
end
function [At, Bt, Uh, Phi] = DMDc(X,Xp,Gam,N_x_transformed,N_u,alpha1,alpha2)
    % SVD of input space
    Om      =   [X; Gam]; 
    [U,S,V]    =   svd(Om,'econ');
    
    alpha   =   alpha1;
    cdS     =   cumsum(diag(S)./sum(diag(S)));
    ptr     =   min(find(cdS>alpha));
    
    Ut      =   U(:,1:ptr);
    St      =   S(1:ptr,1:ptr);
    Vt      =   V(:,1:ptr);
    
    Ut1     =   Ut(1:N_x_transformed,:);
    Ut2     =   Ut(N_x_transformed+N_u:end,:);
    
    % SVD of output space
    [U,S,V]    =   svd(Xp,'econ');
    
    alpha   =   alpha2;
    cdS     =   cumsum(diag(S)./sum(diag(S)));
    rtr     =   min(find(cdS>alpha));
    
    Uh      =   U(:,1:rtr);
    Sh      =   S(1:rtr,1:rtr);
    Vh      =   V(:,1:rtr);
    
    At = Uh'*Xp*Vt*inv(St)*Ut1'*Uh;
    Bt = Uh'*Xp*Vt*inv(St)*Ut2';

    % Compute eigenmodes
%     [W, Lambda] = eig(At);

%     Phi = Xp*Vt*inv(St)*Ut1'*Uh*W;
    Phi = 0; 
end
function n3 = computeDistribution(obj,xt)
    x = obj.U*xt;
    n3 = inverseCoordinateTransformation(obj.N_x,obj.N_x_transformed,obj.N_del,...
                obj.xop,obj.stateScale,obj.inputScale,obj.ext,x);
    n3 = n3(:,end);
end