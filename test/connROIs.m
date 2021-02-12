
function connROIs
    scFilename = 'results/ad-dlcm_ex-cn-roi132.mat';
    roiFilename = '../conn/rois/atlas.nii';
    scf = load(scFilename); % read tractography structure data
    info = niftiinfo(roiFilename);
    v = niftiread(roiFilename);
    roiNum = 132;
    x = []; y = []; z = [];

    % get ROI centroid
    for i=1:roiNum
        BW = logical(zeros(size(v,1),size(v,2),size(v,3)));
        idx = find(v==i);
        BW(idx) = 1;
        s = regionprops3(BW, 'Centroid');
        x = [x; s.Centroid(1,1)];
        y = [y; s.Centroid(1,2)];
        z = [z; s.Centroid(1,3)];
    end
    % re-order ROIs
    ReLOCATE = [1	3	5	7	9	11	13	15	17	19	21	23	25	27	29	31	33	35	37	39	41	43	45	47	50	53	58	60	62	64	66	68	70	72	74	76	78	80	82	84	86	88	90	92	94	96	98	100	102	104	106	108	110	112	114	116	118	120	122	124	49	52	55	56	57	2	4	6	8	10	12	14	16	18	20	22	24	26	28	30	32	34	36	38	40	42	44	46	48	51	54	59	61	63	65	67	69	71	73	75	77	79	81	83	85	87	89	91	93	95	97	99	101	103	105	107	109	111	113	115	117	119	121	123	125	126	127	128	129	130	131	132];
    x = x(ReLOCATE);
    y = y(ReLOCATE);
    z = z(ReLOCATE);
    
    % show ROI centroid
    figure;
    scatter3(x,y,z);

    % show whole connection
    m = max(max(scf.meanWeights));
    sc = scf.meanWeights / m;

    figure;
    hold on;
    scatter3(x,y,z);
    for i=1:roiNum-1
        for j=i+1:roiNum
            x1 = [x(i); x(j)];
            y1 = [y(i); y(j)];
            z1 = [z(i); z(j)];
            s = plot3(x1,y1,z1);
            %{
            if sc(i,j) >= 0
            else
                s = plot3(x1,y1,z1,'b');
            end
            alpha(s,sc(i,j));
            %}
        end
    end
    hold off;
end
