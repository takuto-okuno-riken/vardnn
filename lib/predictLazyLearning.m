%%
% predict of Lazy Learning (multi-input-multi-output)
% input:
%  X            input data matrix (input node x input data)
%  Xt           training data matrix (input node x training data)
%  Yt           training target matrix (output node x training data)
%  LL           LL structure
%  kn           number of k of nearest neighbor (default:1)

function Y = predictLazyLearning(X, Xt, Yt, kn)
    if nargin < 4, kn = 1; end
    inLen = size(X,2);
    outNum = size(Yt,1);

    Y = nan(outNum,inLen);
    for t=1:inLen
        % calc distance between input and data set and find nearest neighbor
        diffs = Xt - repmat(X(:,t),[1 size(Xt,2)]);
        distances = sum(diffs .* diffs,1);
        [B,I] = sort(distances);

        % predict next time step
        if kn==1
            Y(:,t) = Yt(:,I(1));
        else
            % calculate E_LOO error (S.B.Taieb 2012)
            Yq = zeros(outNum,kn);
            E = zeros(outNum,kn);
            ELOO = zeros(1,kn);
            Yq(:,1) = Yt(:,I(1)); E(:,1) = 0; ELOO(1) = 0;
            for k=2:kn
                Yq(:,k) = sum(Yt(:,I(1:k))) / k;
                for j=2:k, E(:,j) = k * (Yt(:,I(j)) - Yq(:,k)) / (k-1); end
                T = sum(E(:,1:k) .* E(:,1:k),2) / k;
                ELOO(k) = mean(T);
            end
            [Be,Ie] = sort(ELOO);
            Kastr = Ie(2);
            Y(:,t) =  Yq(:,Kastr);
        end
    end
end
