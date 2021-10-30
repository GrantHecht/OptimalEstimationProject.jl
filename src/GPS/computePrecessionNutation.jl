
function computePrecessionNutation(ttt, Ax, Ay, As)

    # Define Constants
    ARCSEC2RAD = π / 64800.0
    MUARCSEC2RAD = π / 648000000000

    # Precompute powers of ttt
    ttt2 = ttt*ttt
    ttt3 = ttt2*ttt
    ttt4 = ttt3*ttt

    # Compute Delaunay Arguments
    DelArgs = @SVector[
        rem((485868.249036 + 1717915923.2178*ttt + 31.8792*ttt2 + 0.051635*ttt3 - 0.00024470*ttt4)*ARCSEC2RAD,2*pi),
        rem((1287104.79305 + 129596581.0481*ttt - 0.5532*ttt2 + 0.000136*ttt3 - 0.00001149*ttt4)*ARCSEC2RAD,2*pi),
        rem((335779.526232 + 1739527262.8478*ttt - 12.7512*ttt2 - 0.001037*ttt3 + 0.00000417*ttt4)*ARCSEC2RAD,2*pi),
        rem((1072260.70369 + 16029616001.2090*ttt - 6.3706*ttt2 + 0.006593*ttt3 - 0.00003169*ttt4)*ARCSEC2RAD,2*pi),
        rem((450160.398036 - 6962890.5431*ttt + 7.4722*ttt2 + 0.007702*ttt3 - 0.00005939*ttt4)*ARCSEC2RAD,2*pi),
        rem(4.402608842 + 2608.7903141574*ttt,2*pi),
        rem(3.176146698 + 1021.3285546211*ttt,2*pi),
        rem(1.753470314 + 628.3075849991*ttt,2*pi),
        rem(6.203480913 + 334.0612426700*ttt,2*pi),
        rem(0.599546497 + 52.9690962641*ttt,2*pi),
        rem(0.874016757 + 21.3299104960*ttt,2*pi),
        rem(5.481293872 + 7.4781598567*ttt,2*pi),
        rem(5.311886287 + 3.8133035638*ttt,2*pi),
        rem(0.02438175*ttt + 0.00000538691*ttt2,2*pi)]

    # Compute terms
    XT0 = 0.0; YT0 = 0.0; ST0 = 0.0
    XT1 = 0.0; YT1 = 0.0; ST1 = 0.0 
    XT2 = 0.0; YT2 = 0.0; ST2 = 0.0
    XT3 = 0.0; YT3 = 0.0; ST3 = 0.0
    XT4 = 0.0; YT4 = 0.0; ST4 = 0.0
    for i in 1:15
        if i == 1
            for j in 1:1306
                ap   = transpose(view(Ax,j,4:17))*DelArgs
                XT0 += Ax[j,2]*sin(ap) + Ax[j,3]*cos(ap)
            end
        elseif i == 2
            for j in 1:253
                ap   = transpose(view(Ax,1306 + j,4:17))*DelArgs;
                XT1 += Ax[1306 + j,2]*sin(ap) + Ax[1306 + j,3]*cos(ap);
            end
        elseif i == 3
            for j in 1:36
                ap   = transpose(view(Ax,1559 + j,4:17))*DelArgs;
                XT2 += Ax[1559 + j,2]*sin(ap) + Ax[1559 + j,3]*cos(ap);
            end
        elseif i == 4
            for j = 1:4
                ap   = transpose(view(Ax,1595 + j,4:17))*DelArgs;
                XT3 += Ax[1595 + j,2]*sin(ap) + Ax[1595 + j,3]*cos(ap);
            end
        elseif i == 5
            ap   = transpose(view(Ax,1600,4:17))*DelArgs;
            XT4 += Ax[1600,2]*sin(ap) + Ax[1600,3]*cos(ap);
        elseif i == 6
            for j in 1:962
                ap   = transpose(view(Ay,j,4:17))*DelArgs;
                YT0 += Ay[j,2]*sin(ap) + Ay[j,3]*cos(ap);
            end
        elseif i == 7
            for j in 1:277
                ap   = transpose(view(Ay,962 + j,4:17))*DelArgs;
                YT1 += Ay[962 + j,2]*sin(ap) + Ay[962 + j,3]*cos(ap);
            end
        elseif i == 8
            for j in 1:30
                ap   = transpose(view(Ay,1239 + j,4:17))*DelArgs;
                YT2 += Ay[1239 + j,2]*sin(ap) + Ay[1239 + j,3]*cos(ap);
            end
        elseif i == 9
            for j in 1:5
                ap   = transpose(view(Ay,1269 + j,4:17))*DelArgs;
                YT3 += Ay[1269 + j,2]*sin(ap) + Ay[1269 + j,3]*cos(ap);
            end
        elseif i == 10
            ap   = transpose(view(Ay,1275,4:17))*DelArgs;
            YT4 += Ay[1275,2]*sin(ap) + Ay[1275,3]*cos(ap);
        elseif i == 11
            for j in 1:33
                ap   = transpose(view(As,j,4:17))*DelArgs;
                ST0 += As[j,2]*sin(ap) + As[j,3]*cos(ap);
            end
        elseif i == 12
            for j in 1:3
                ap   = transpose(view(As,33 + j,4:17))*DelArgs;
                ST1 += As[33 + j,2]*sin(ap) + As[33 + j,3]*cos(ap);
            end 
        elseif i == 13
            for j in 1:25
                ap   = transpose(view(As,36 + j,4:17))*DelArgs;
                ST2 += As[36 + j,2]*sin(ap) + As[36 + j,3]*cos(ap);
            end 
        elseif i == 14
            for j in 1:4
                ap  = transpose(view(As,61 + j,4:17))*DelArgs;
                ST3 += As[61 + j,2]*sin(ap) + As[61 + j,3]*cos(ap);
            end 
        else
            ap   = transpose(view(As,66,4:17))*DelArgs;
            ST4 += As[66,2]*sin(ap) + As[66,3]*cos(ap);
        end
    end

    # Get pole corrdinate corrections 
    JD = ttt*36525 + 2451545.0
    ΔX = getdx(JD)
    ΔY = getdy(JD)

    # Compute X, Y, and S
    X = (-0.016617 + 2004.191898*ttt - 0.4297829*ttt2 - 0.19861834*ttt3 + 
         0.000007578*ttt4 + 0.0000059285*ttt4*ttt)*ARCSEC2RAD + 
         (XT0 + XT1*ttt + XT2*ttt2 + XT3*ttt3 + XT4*ttt4)*MUARCSEC2RAD + ΔX;

    Y = (-0.006951 - 0.025896*ttt - 22.4072747*ttt2 + 0.00190059*ttt3 + 
         0.001112526*ttt4 + 0.0000001358*ttt4*ttt)*ARCSEC2RAD + 
         (YT0 + YT1*ttt + YT2*ttt2 + YT3*ttt3 + YT4*ttt4)*MUARCSEC2RAD + ΔY;

    S = -(X*Y)/2 + (0.000094 + 0.00380865*ttt - 0.00012268*ttt2 - 0.072574*ttt3 + 
         0.00002798*ttt4 + 0.00001562*ttt4*ttt)*ARCSEC2RAD + 
         (ST0 + ST1*ttt + ST2*ttt2 + ST3*ttt3 + ST4*ttt4)*MUARCSEC2RAD;

    # Compute rotation matrix
    d = atan(sqrt((X^2 + Y^2)/(1 - X^2 - Y^2)))
    a = 1/(1 + cos(d))

    m1 = @SMatrix [1-a*X^2  -a*X*Y  X;
                   -a*X*Y  1-a*Y^2  Y;
                   -X  -Y  1-a*(X^2+Y^2)]
    m2 = @SMatrix [ cos(S) sin(S) 0;
                   -sin(S) cos(S) 0;
                    0 0 1]
    PN = SMatrix{3,3}(m1*m2)

    return PN
end