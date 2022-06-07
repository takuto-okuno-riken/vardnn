function testXcorr
    maxlag = 0;
    N = 3;
    X = rand(100,N);
    Y = rand(100,N);
    Xm = X - mean(X);
    Ym = Y - mean(Y);
    XY = [Xm, Ym];
    C0 = xcorr(XY,maxlag,'normalized');
    C1 = zeros(N*2,N*2);
    for i=1:N*2
        for j=1:N*2
            C1(i,j) = xcorr(XY(:,i),XY(:,j),maxlag,'normalized');
        end
    end
    C2 = C1(:)';
    s = sum(abs(C0-C2),'all')
    if s >= 1e-10, disp('bad for loop with C0'); else, disp('for loop with C0 : ok'); end

    % check N vs. 1 pattern
    a1 = xcorrMat(Xm,Xm(:,1),maxlag,'normalized');
    a3 = xcorrMat(Xm,Xm,maxlag,'normalized');
    s = sum(abs(a1-a3(1:N)));
    if s >= 1e-10, disp('bad a1 vs. a3'); else, disp('for a1 vs. a3 : ok'); end

    % check xcorrMat
    a = xcorrMat(Xm,Xm,maxlag,'normalized');
    b11 = reshape(a,N,N);
    a = xcorrMat(Xm,Ym,maxlag,'normalized');
    b21 = reshape(a,N,N);
    a = xcorrMat(Ym,Xm,maxlag,'normalized');
    b12 = reshape(a,N,N);
    a = xcorrMat(Ym,Ym,maxlag,'normalized');
    b22 = reshape(a,N,N);
    B = [b11, b12; b21, b22];
    S = C1 - B;
    s = sum(S,'all');
    if s >= 1e-10, disp('bad xcorrMat with C1'); else, disp('xcorrMat with C1 : ok'); end

    C3 = xcorrMatSep(XY,XY,maxlag,'normalized',1);
    s = sum(abs(C0-C3),'all');
    if s >= 1e-10, disp('bad xcorrMatSep1 with C0'); else, disp('xcorrMatSep1 with C0 : ok'); end
    
    C3 = xcorrMatSep(XY,XY,maxlag,'normalized',2);
    s = sum(abs(C0-C3),'all');
    if s >= 1e-10, disp('bad xcorrMatSep with C0'); else, disp('xcorrMatSep with C0 : ok'); end

    maxlag = 1;
    C0 = xcorr(XY,maxlag,'normalized');
    C3 = xcorrMatSep(XY,XY,maxlag,'normalized',2);
    s = sum(abs(C0-C3),'all')
    if s >= 1e-10, disp('bad xcorrMatSep with C0 (m1)'); else, disp('xcorrMatSep with C0 (m1) : ok'); end


    % check odd pattern
    maxlag = 0;
    X = rand(100,7);
    Xm = X - mean(X);
    C1 = xcorr(Xm,maxlag,'normalized');
    C2 = xcorrMat(Xm,Xm,maxlag,'normalized');
    s = sum(abs(C1-C2),'all');
    if s >= 1e-10, disp('bad xcorrMat with C1 (odd)'); else, disp('xcorrMat with C1 (odd) : ok'); end

    C3 = xcorrMatSep(Xm,Xm,maxlag,'normalized',2);
    s = sum(abs(C1-C3),'all');
    if s >= 1e-10, disp('bad xcorrMatSep with C0 (odd)'); else, disp('xcorrMatSep with C1 (odd) : ok'); end

    maxlag = 1;
    C1 = xcorr(Xm,maxlag,'normalized');
    C2 = xcorrMat(Xm,Xm,maxlag,'normalized');
    s = sum(abs(C1-C2),'all');
    if s >= 1e-10, disp('bad xcorrMat with C1 (odd2)'); else, disp('xcorrMat with C1 (odd2) : ok'); end

    C3 = xcorrMatSep(Xm,Xm,maxlag,'normalized',2);
    s = sum(abs(C1-C3),'all')
    if s >= 1e-10, disp('bad xcorrMatSep with C0 (odd2)'); else, disp('xcorrMatSep with C1 (odd2) : ok'); end


    % check sep 3
    N = 8;
    X = rand(100,N);
    maxlag = 1;
    C0 = xcov(X,maxlag,'normalized');
    Xm = X - mean(X);
    C3 = xcorrMatSep(Xm,Xm,maxlag,'normalized',3);
    s = sum(abs(C0-C3),'all')
    if s >= 1e-10, disp('bad xcorrMatSep with C0 (r2)'); else, disp('xcorrMatSep with C0 (r2) : ok'); end


    % check auto-correlation (a bit big)
    X = rand(200,101);
    maxlag = 2;
    C0 = xcov(X,maxlag,'normalized');
    Xm = X - mean(X);
    C1 = xcorr(Xm,maxlag,'normalized');
    C2 = xcorrMat(Xm,Xm,maxlag,'normalized');
    s = sum(abs(C0-C2),'all');
    if s >= 1e-10, disp('bad xcorrMat with C0'); else, disp('xcorrMat with C0 : ok'); end
    s = sum(abs(C1-C2),'all');
    if s >= 1e-10, disp('bad xcorrMat with C1'); else, disp('xcorrMat with C1 : ok'); end
    C3 = xcorrMatSep(Xm,Xm,maxlag,'normalized',2);
    s = sum(abs(C0-C3),'all');
    if s >= 1e-10, disp('bad xcorrMatSep with C0'); else, disp('xcorrMatSep with C0 : ok'); end

    X = rand(1000,1551);
    maxlag = 3;
    tic;
    C0 = xcov(X,maxlag,'normalized');
    toc
    Xm = X - mean(X);
    tic;
    C3 = xcorrMatSep(Xm,Xm,maxlag,'normalized',5);
    toc
    s = sum(abs(C0-C3),'all');
    if s >= 1e-10, disp('bad xcorrMatSep with C0 (big)'); else, disp('xcorrMatSep with C0 (big) : ok'); end
end
