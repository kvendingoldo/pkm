Window = figure("figure_name", "GraphApp", "BackgroundColor", [1, 1, 1])
set(Window, "Position", [300, 50, 900, 1000])

eqappLabel = uicontrol("Style", "text", "Position", [20, 580, 100, 35], "String", "The equation: ", "BackgroundColor", [1, 1, 1])
eqLabel = uicontrol("Style", "text", "Position", [130, 570, 270, 55], "String", "$\frac{d^2\varphi}{dx^2} + \alpha\frac{d\varphi}{dx} + \beta\varphi = 0$", "BackgroundColor", [1, 1, 1])

alphaLabel = uicontrol("Style", "text", "Position", [20, 540, 70, 25], "String", "$\alpha = $", "BackgroundColor", [1, 1, 1])
alphaEdit = uicontrol("Style", "edit", "Position", [60, 540, 70, 25], "String", "-6")

betaLabel = uicontrol("Style", "text", "Position", [145, 540, 70, 25], "String", "$\beta = $", "BackgroundColor", [1, 1, 1])
betaEdit = uicontrol("Style", "edit", "Position", [185, 540, 70, 25], "String", "9")

MLabel = uicontrol("Style", "text", "Position", [272, 540, 70, 25], "String", "$M = $", "BackgroundColor", [1, 1, 1])
MEdit = uicontrol("Style", "edit", "Position", [320, 540, 70, 25], "String", "3")

DirichletLabel = uicontrol("Style", "text", "Position", [20, 510, 154, 25], "String", "Dirichlet problem: ", "BackgroundColor", [1, 1, 1])
DirichletAppLabel = uicontrol("Style", "text", "Position", [20, 490, 210, 25], "String", "$\varphi(0) = 0,\ \varphi(1) = 1$", "BackgroundColor", [1, 1, 1], "fontsize", 15)
DirichletRadioButton = uicontrol("Style", "radiobutton", "Position", [250, 500, 20, 25], "BackgroundColor", [1, 1, 1], "Value", 1, "Callback", "DirichletRadioButtonCallback")

NeumannLabel1 = uicontrol("Style", "text", "Position", [20, 460, 275, 25], "String", "Neumann problem: ", "BackgroundColor", [1, 1, 1])
NeumannAppLabel2 = uicontrol("Style", "text", "Position", [20, 440, 275, 25], "String", "$\varphi(0) = 0,\ (d\varphi/dx)_{x = 1} = 1$", "BackgroundColor", [1, 1, 1], "fontsize", 15)
NeumannRadioButton = uicontrol("Style", "radiobutton", "Position", [250, 450, 20, 25], "BackgroundColor", [1, 1, 1], "Value", 0, "Callback", "NeumannRadioButtonCallback")

DirichletDiffLabel = uicontrol("Style", "text", "Position", [300, 510, 280, 25], "String", "Error for the solution of", "BackgroundColor", [1, 1, 1]);
DirichletDiffLabel = uicontrol("Style", "text", "Position", [300, 490, 280, 25], "String", "Dirichlet problem", "BackgroundColor", [1, 1, 1]);
DirichletDiffRadioButton = uicontrol("Style", "radiobutton", "Position", [500, 500, 20, 25], "BackgroundColor", [1, 1, 1], "Value", 0, "Callback", "DDRadioButtonCallback")

NeumannDiffLabel1 = uicontrol("Style", "text", "Position", [300, 460, 275, 25], "String", "Error for the solution of", "BackgroundColor", [1, 1, 1]);
NeumannDiffLabel2 = uicontrol("Style", "text", "Position", [300, 440, 275, 25], "String", "Neumann problem ", "BackgroundColor", [1, 1, 1]);
NeumannDiffRadioButton = uicontrol("Style", "radiobutton", "Position", [500, 450, 20, 25], "BackgroundColor", [1, 1, 1], "Value", 0,  "Callback", "NDRadioButtonCallback")

realSolve = uicontrol("Style", "text", "Position", [20, 400, 175, 25], "String", "The exact solution", "BackgroundColor", [1, 1, 1])
realSolveCheckBox = uicontrol("Style", "checkbox", "Position", [155, 403, 20, 19], "BackgroundColor", [0, 0, 0], "Value", 1)

method1Label = uicontrol("Style", "text", "Position", [20, 370, 195, 25], "String", "1) The pointwise collocation method", "BackgroundColor", [1, 1, 1])
method1CheckBox = uicontrol("Style", "checkbox", "Position", [275, 373, 20, 19], "BackgroundColor", [1, 0, 0], "Value", 1)

method2Label = uicontrol("Style", "text", "Position", [20, 340, 255, 25], "String", "2) The collocation method by region", "BackgroundColor", [1, 1, 1])
method2CheckBox = uicontrol("Style", "checkbox", "Position", [275, 343, 20, 19], "BackgroundColor", [0, 1, 0], "Value", 1)

method3Label = uicontrol("Style", "text", "Position", [20, 310, 175, 25], "String", "3) The Galerkin Method", "BackgroundColor", [1, 1, 1])
method3CheckBox = uicontrol("Style", "checkbox", "Position", [275, 313, 20, 19], "BackgroundColor", [0, 0, 1], "Value", 1)

method4Label = uicontrol("Style", "text", "Position", [310, 370, 220, 25], "String", "4) The method 1 without", "BackgroundColor", [1, 1, 1])
method4psiLabel = uicontrol("Style", "text", "Position", [485, 373, 40, 25], "String", "$\psi$", "BackgroundColor", [1, 1, 1], "fontsize", 17)
method4CheckBox = uicontrol("Style", "checkbox", "Position", [510, 373, 20, 19], "BackgroundColor", [1, 0, 1], "Value", 0)

method5Label = uicontrol("Style", "text", "Position", [310, 340, 220, 25], "String", "5) The method 2 without ", "BackgroundColor", [1, 1, 1])
method5psiLabel = uicontrol("Style", "text", "Position", [485, 343, 40, 25], "String", "$\psi$", "BackgroundColor", [1, 1, 1], "fontsize", 17)
method5CheckBox = uicontrol("Style", "checkbox", "Position", [510, 343, 20, 19], "BackgroundColor", [1, 0.7, 0], "Value", 0)

method6Label = uicontrol("Style", "text", "Position", [310, 310, 220, 25], "String", "6) The method 3 without ", "BackgroundColor", [1, 1, 1])
method6psiLabel = uicontrol("Style", "text", "Position", [485, 313, 40, 25], "String", "$\psi$", "BackgroundColor", [1, 1, 1], "fontsize", 17)
method6CheckBox = uicontrol("Style", "checkbox", "Position", [510, 313, 20, 19], "BackgroundColor", [0, 0.7, 1], "Value", 0)

method7Label = uicontrol("Style", "text", "Position", [550, 370, 220, 25], "String", "7) The weak formulation", "BackgroundColor", [1, 1, 1])
method7CheckBox = uicontrol("Style", "checkbox", "Position", [770, 373, 20, 19], "BackgroundColor", [1, 0.5, 0], "Value", 0)
set(method7CheckBox, "Enable", "off");

method8Label = uicontrol("Style", "text", "Position", [550, 340, 220, 25], "String", "8) The finite element method", "BackgroundColor", [1, 1, 1])
method8CheckBox = uicontrol("Style", "checkbox", "Position", [770, 343, 20, 19], "BackgroundColor", [0.0, 0.7, 0.7], "Value", 0)

Button = uicontrol("Style", "pushbutton", "Position", [360, 280, 70, 25], "String", "Run", "Callback", "run");

q = 101;
points = 0:(1/(q - 1)):1;
diffMode = 0;

function DirichletRadioButtonCallback()
	set(DirichletRadioButton, "Value", 1);
    set(NeumannRadioButton, "Value", 0);
    set(DirichletDiffRadioButton, "Value", 0);
    set(NeumannDiffRadioButton, "Value", 0);
    set(realSolveCheckBox, "Value", 1);
    set(realSolveCheckBox, "Enable", "on");
    set(method7CheckBox, "Value", 0);
    set(method7CheckBox, "Enable", "off");
endfunction;

function NeumannRadioButtonCallback()
    set(NeumannRadioButton, "Value", 1);
	set(DirichletRadioButton, "Value", 0);
    set(DirichletDiffRadioButton, "Value", 0);
    set(NeumannDiffRadioButton, "Value", 0);
    set(realSolveCheckBox, "Value", 1);
    set(realSolveCheckBox, "Enable", "on");
    set(method7CheckBox, "Enable", "on");
endfunction;

function DDRadioButtonCallback()
    set(DirichletDiffRadioButton, "Value", 1);
    set(NeumannDiffRadioButton, "Value", 0);
    set(DirichletRadioButton, "Value", 0);
	set(NeumannRadioButton, "Value", 0);
    set(realSolveCheckBox, "Value", 0);
    set(realSolveCheckBox, "Enable", "off");
    set(method7CheckBox, "Value", 0);
    set(method7CheckBox, "Enable", "off");
endfunction;

function NDRadioButtonCallback()
    set(DirichletDiffRadioButton, "Value", 0);
    set(NeumannDiffRadioButton, "Value", 1);
    set(DirichletRadioButton, "Value", 0);
	set(NeumannRadioButton, "Value", 0);
    set(realSolveCheckBox, "Value", 0);
    set(realSolveCheckBox, "Enable", "off");
    set(method7CheckBox, "Enable", "on");
endfunction;

function phi = phiValueD(x, aArray, MValue)
    phi = x;
	for k = 1:1:MValue
		phi = phi + aArray(k)*(1 - x)*x^k;
	end;
endfunction;

function phi = phiValueN(x, aArray, MValue)
    phi = x;
	for k = 1:1:MValue
        phi = phi + aArray(k)*(x^(k + 1)/(k + 1) - x^k/k);
	end;
endfunction

function phi = phiDiffValueN(x, aArray, MValue)
    phi = 1;
    for k = 1:1:MValue
        phi = phi + aArray(k)*(x^k - x^(k - 1));
    end
endfunction

function phi = phiValueWithoutPsiD(x, aArray, MValue)
    phi = 0;
    for k = 1:1:MValue
		phi = phi + aArray(k)*x^k;
	end;
endfunction;

function phi = phiValueWithoutPsiN(x, aArray, MValue)
    phi = 0;
    for k = 1:1:MValue
		phi = phi + aArray(k)*x^k;
	end;
endfunction;

function phi = phiDiffValueWithoutPsiN(x, aArray, MValue)
    phi = 0;
    for k = 1:1:MValue
        phi = phi + aArray(k)*k*(x^(k - 1));
    end
endfunction;

function Nk = getNkOfFiniteElements(kValue, xValue, MValue)
     MPoints = 0:(1/MValue):1;Nk = 0;
     if kValue = 1 then
         kValue = 2;
     end
     if (xValue <= MPoints(kValue - 1)) then
         Nk = 0;
     elseif ((xValue > MPoints(kValue - 1)) & (xValue < MPoints(kValue)))
         Nk = (xValue - MPoints(kValue - 1)) / (MPoints(kValue) - MPoints(kValue - 1));
     elseif (xValue == MPoints(k))
         Nk = 1;
     elseif ((xValue > MPoints(kValue)) & (xValue < MPoints(kValue + 1)))
        Nk = (MPoints(kValue + 1) - xValue) / (MPoints(kValue + 1) - MPoints(kValue));
    elseif (xValue > MPoints(kValue + 1))
        Nk = 0;
    end
endfunction

function phi = phiValueOfFiniteElements(x, aArray, MValue)
    phi = 0;
    for k = 1:1:(MValue + 1)
        phi = phi + aArray(k) * getNkOfFiniteElements(k,  x, MValue);
    end
endfunction

function phiArray = analiticSolutionD()
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = points(i)*exp(3*points(i))/exp(3);
	end;
endfunction;

function phiArray = analiticSolutionN()
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = exp(3*points(i) - 3)*points(i)/4;
	end;
endfunction

function realSolutionArray = DirichletODE2(alphaValue, betaValue, points)
    len = length(points);
    h = 1/(len - 1);
    coeffsArray = zeros(len, len);
    freeArray = zeros(len, 1);
    for s = 1:1:len
        for k = 1:1:len
            if (s == 1) then
                if k == 2 then
                    coeffsArray(s, k) = 2*betaValue*h^2 - 4;
                    coeffsArray(s, k + 1) = 2 + alphaValue*h;
                end
            elseif s == len
                if k == len - 1 then
                    coeffsArray(s, k) = 2*betaValue*h^2 - 4;
                    coeffsArray(s, k - 1) = 2 - alphaValue*h;
                    freeArray(s) = (+1)*(2 + alphaValue*h);
                end
            elseif s == (k + 1)
                if k < len - 1 then
                    if k > 1 then
                        coeffsArray(s, k - 1) = 2 - alphaValue*h;
                        coeffsArray(s, k) = 2*betaValue*h^2 - 4;
                        coeffsArray(s, k + 1) = 2 + alphaValue*h;
                    end
                end
            end
        end
    end
    solution = linsolve(coeffsArray, freeArray);
    realSolutionArray = ones(len, 1);
    realSolutionArray(1) = 0;
    realSolutionArray(len) = 1;
    for i = 2:1:len - 1
        realSolutionArray(i) = solution(i);
    end
endfunction;

function phiArray = methodD1(alphaValue ,betaValue, MValue)
    h = 1/(MValue + 1);
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
		fArray(s) = (+1)*(alphaValue + betaValue*xs);
		for k = 1:1:MValue
			KArray(s, k) = k*(k - 1)*xs^(k - 2) - k*(k + 1)*xs^(k - 1) + alphaValue*k*xs^(k - 1) - alphaValue*(k + 1)*xs^k + betaValue*(1 - xs)*xs^k;
		end;
	end;
    disp(1);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = phiValueD(points(i), aArray, MValue);
	end;
endfunction;

function phiArray = methodN1(alphaValue ,betaValue, MValue)
    h = 1/(MValue + 1);
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
		fArray(s) = (+1)*(alphaValue + betaValue*xs);
		for k = 1:1:MValue
			KArray(s, k) = k*xs^(k - 1) - (k - 1)*xs^(k - 2) + alphaValue*xs^k - alphaValue*xs^(k - 1) + (betaValue*xs^(k + 1))/(k + 1) - (betaValue*xs^k)/k;
		end;
	end;
    disp(1);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
        if diffMode == 1 then
            phiArray(i) = phiDiffValueN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueN(points(i), aArray, MValue);
        end
	end;
endfunction;

function phiArray = methodD2(alphaValue ,betaValue, MValue)
    h = 1/MValue;
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
		xsl = (s - 1)*h;
		fArray(s) = alphaValue*xs - alphaValue*xsl + (betaValue*xs^2)/2 - (betaValue*xsl^2)/2;
		for k = 1:1:MValue
			KArray(s, k) = k*xs^(k - 1) - k*xsl^(k - 1) - (k + 1)*xs^k + (k + 1)*xsl^k + alphaValue*xs^k - alphaValue*xsl^k - alphaValue*xs^(k + 1) + alphaValue*xsl^(k + 1) + (betaValue*xs^(k + 1))/(k + 1) - (betaValue*xsl^(k + 1))/(k + 1) - (betaValue*xs^(k + 2))/(k + 2) + (betaValue*xsl^(k + 2))/(k + 2);
		end;
	end;
    disp(2);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = phiValueD(points(i), aArray, MValue);
	end;
endfunction;

function phiArray = methodN2(alphaValue ,betaValue, MValue)
    h = 1/MValue;
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
		xsl = (s - 1)*h;
        fArray(s) = alphaValue*xs - alphaValue*xsl + (betaValue*xs^2)/2 - (betaValue*xsl^2)/2;
		for k = 1:1:MValue
            KArray(s, k) = xs^k - xsl^k - xs^(k - 1) + xsl^(k - 1) + (alphaValue*xs^(k + 1))/(k + 1) - (alphaValue*xsl^(k + 1))/(k + 1) - (alphaValue*xs^k)/k + (alphaValue*xsl^k)/k + (betaValue*xs^(k + 2))/((k + 1)*(k + 2)) - (betaValue*xsl^(k + 2))/((k + 1)*(k + 2)) - (betaValue*xs^(k + 1))/(k*(k + 1)) + (betaValue*xsl^(k + 1))/(k*(k + 1));
		end;
	end;
    disp(2);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		if diffMode == 1 then
            phiArray(i) = phiDiffValueN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueN(points(i), aArray, MValue);
        end
	end;
endfunction;

function phiArray = methodD3(alphaValue, betaValue, MValue)
    fArray = zeros(MValue, 1);
    KArray = zeros(MValue, MValue);
    for s = 1:1:MValue
        fArray(s) = (+1)*((alphaValue/(s + 1) + (betaValue - alphaValue)/(s + 2) - betaValue/(s + 3)));
        for k = 1:1:MValue
            KArray(s, k) = alphaValue*k/(k + s) - alphaValue/(k + s + 1) - 2*alphaValue*k/(k + s + 1) + alphaValue/(k + s + 2) + alphaValue*k/(k + s + 2) + betaValue/(k + s + 1) - 2*betaValue/(k + s + 2) + betaValue/(k + s + 3) + (k^2)/(k + s - 1) - (2*k^2)/(k + s) + (k^2)/(k + s + 1) - k/(k + s - 1) + k/(k + s + 1);
        end
    end
    disp(3);disp(det(KArray));
    aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = phiValueD(points(i), aArray, MValue);
	end;
endfunction;

function phiArray = methodN3(alphaValue, betaValue, MValue)
    fArray = zeros(MValue, 1);
    KArray = zeros(MValue, MValue);
    for s = 1:1:MValue
        fArray(s) = -alphaValue/(s*(s + 1)) + alphaValue/((s + 1)*(s + 2)) - betaValue/((s + 1)*(s + 2)) - betaValue/(s*(s + 1)*(s + 2)) + betaValue/((s + 1)*(s + 3));
        for k = 1:1:MValue
            KArray(s, k) = alphaValue/((s + 1)*(k + s)) + alphaValue/(s*(s + 1)*(k + s)) - 2*alphaValue/((s + 1)*(k + s + 1)) - alphaValue/(s*(s + 1)*(k + s + 1)) + alphaValue/((k + 1)*(s + 1)*(k + s + 2)) + alphaValue*k/((k + 1)*(s + 1)*(k + s + 2)) + betaValue/(k*(s + 1)*(k + s + 1)) + betaValue/(k*s*(s + 1)*(k + s + 1)) - 2*betaValue/((k + 1)*(s + 1)*(k + s + 2)) - betaValue/(k*(k + 1)*(s + 1)*(k + s + 2)) - betaValue/(s*(k + 1)*(s + 1)*(k + s + 2)) + betaValue/((k + 1)*(s + 1)*(k + s + 3)) + k/(s*(k + s - 1)) - 1/(s*(k + s - 1)) - 2*k/((s + 1)*(k + s)) - k/(s*(s + 1)*(k + s)) + 1/((s + 1)*(k + s)) + k/((s + 1)*(k + s + 1));
        end
    end
    disp(3);disp(det(KArray));
    aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		if diffMode == 1 then
            phiArray(i) = phiDiffValueN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueN(points(i), aArray, MValue);
        end
	end;
endfunction;

function phiArray = methodD4(alphaValue, betaValue, MValue)
    h = 1/(MValue);
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
        fArray(s) = 0;
        if s == MValue then
            fArray(s) = (+1);
        end
		for k = 1:1:MValue
            KArray(s, k) = k*(k - 1)*xs^(k - 2) + alphaValue*k*xs^(k - 1) + betaValue*xs^(k);
            if s == MValue then
                KArray(s, k) = KArray(s, k) - 1;
            end
		end;
	end;
    disp(4);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = phiValueWithoutPsiD(points(i), aArray, MValue);
	end;
endfunction

function phiArray = methodN4(alphaValue, betaValue, MValue)
    h = 1/(MValue);
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
        fArray(s) = 0;
        if s == MValue then
            fArray(s) = (+1);
        end
		for k = 1:1:MValue
            KArray(s, k) = k*(k - 1)*xs^(k - 2) + alphaValue*k*xs^(k - 1) + betaValue*xs^(k);
            if s == MValue then
                KArray(s, k) = KArray(s, k) - k;
            end
		end;
	end;
    disp(4);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
        if diffMode == 1 then
            phiArray(i) = phiDiffValueWithoutPsiN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueWithoutPsiN(points(i), aArray, MValue);
        end
	end;
endfunction

function phiArray = methodD5(alphaValue, betaValue, MValue)
    h = 1/MValue;
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
		xsl = (s - 1)*h;
        if s == MValue then
            fArray(s) = (+1);
        end
		for k = 1:1:MValue
            KArray(s, k) = k*xs^(k - 1) - k*xsl^(k - 1) + alphaValue*xs^k - alphaValue*xsl^k + betaValue*(xs^(k + 1))/(k + 1) - betaValue*(xsl^(k + 1))/(k + 1);
            if s == MValue then
                KArray(s, k) = KArray(s, k) - 1;
            end
		end;
	end;
    disp(5);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = phiValueWithoutPsiD(points(i), aArray, MValue);
	end;
endfunction

function phiArray = methodN5(alphaValue, betaValue, MValue)
    h = 1/MValue;
	fArray = zeros(MValue, 1);
	KArray = zeros(MValue, MValue);
	for s = 1:1:MValue
		xs = s*h;
		xsl = (s - 1)*h;
        if s == MValue then
            fArray(s) = (+1);
        end
		for k = 1:1:MValue
            KArray(s, k) = k*xs^(k - 1) - k*xsl^(k - 1) + alphaValue*xs^k - alphaValue*xsl^k + betaValue*(xs^(k + 1))/(k + 1) - betaValue*(xsl^(k + 1))/(k + 1);
            if s == MValue then
                KArray(s, k) = KArray(s, k) - k;
            end
		end;
	end;
    disp(5);disp(det(KArray));
	aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
	phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		if diffMode == 1 then
            phiArray(i) = phiDiffValueWithoutPsiN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueWithoutPsiN(points(i), aArray, MValue);
        end
	end;
endfunction

function phiArray = methodD6(alphaValue, betaValue, MValue)
    fArray = zeros(MValue, 1);
    KArray = zeros(MValue, MValue);
    for s = 1:1:MValue
        fArray(s) = (+1);
        for k = 1:1:MValue
            KArray(s, k) = k*(k - 1)/(k + s - 1) + alphaValue*k/(k + s) + betaValue/(k + s + 1) - 1;
        end
    end
    disp(6);disp(det(KArray));
    aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		phiArray(i) = phiValueWithoutPsiD(points(i), aArray, MValue);
	end;
endfunction

function phiArray = methodN6(alphaValue, betaValue, MValue)
    fArray = zeros(MValue, 1);
    KArray = zeros(MValue, MValue);
    for s = 1:1:MValue
        fArray(s) = (+1);
        for k = 1:1:MValue
            KArray(s, k) = k*(k - 1)/(k + s - 1) + alphaValue*k/(k + s) + betaValue/(k + s + 1) - k;
        end
    end
    disp(6);disp(det(KArray));
    aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		if diffMode == 1 then
            phiArray(i) = phiDiffValueWithoutPsiN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueWithoutPsiN(points(i), aArray, MValue);
        end
	end;
endfunction

function phiArray = weakFormulation(alphaValue, betaValue, MValue)
    fArray = zeros(MValue, 1);
    KArray = zeros(MValue, MValue);
    for s = 1:1:MValue
        fArray(s) = (+1);
        for k = 1:1:MValue
            KArray(s, k) = -(k*s)/(k + s - 1) + alphaValue*k/(k + s) + betaValue/(k + s + 1);
        end
    end
    disp(7);disp(det(KArray));
    aArray = linsolve(KArray, fArray);
    aArray = -inv(KArray)* fArray;
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
		if diffMode == 1 then
            phiArray(i) = phiDiffValueWithifoutPsiN(points(i), aArray, MValue);
        else
            phiArray(i) = phiValueWithoutPsiN(points(i), aArray, MValue);
        end
	end;
endfunction

function phiArray = finiteElementsD(alphaValue, betaValue, MValue)
    fArray = zeros(MValue + 1, 1);
    fArray(MValue + 1, 1) = 1;

    h = 1/MValue;
    a = -1/h - alphaValue/2 + betaValue*h/3;
    b = 1/h + alphaValue/2 + betaValue*h/6;
    c = 1/h - alphaValue/2 + betaValue*h/6;
    d = -1/h + alphaValue/2 + betaValue*h/3;
    ad = a + d;

    KArray = zeros(MValue + 1, MValue + 1);
    KArray(1,1) = 1;
    KArray(MValue + 1, MValue + 1) = 1;
    for i=2:1:MValue
        KArray(i,i-1) = c;
        KArray(i,i)= ad;
        KArray(i,i+1) = b;
    end
    disp(8);disp(KArray);disp(fArray);
    aArray = linsolve(KArray, -fArray);
    disp(aArray);
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
        phiArray(i) = phiValueOfFiniteElements(points(i), aArray, MValue);
	end;
endfunction

function phiArray = finiteElementsN(alphaValue, betaValue, MValue)
    fArray = zeros(MValue + 1, 1);
    fArray(MValue + 1, 1) = -1;

    h = 1/MValue;
    a = -1/h - alphaValue/2 + betaValue*h/3;
    b = 1/h + alphaValue/2 + betaValue*h/6;
    c = 1/h - alphaValue/2 + betaValue*h/6;
    d = -1/h + alphaValue/2 + betaValue*h/3;
    ad = a + d;

    KArray = zeros(MValue + 1, MValue + 1);
    KArray(1,1) = 1;
    KArray(MValue + 1, MValue) = c;
    KArray(MValue + 1, MValue + 1) = d;
    for i=2:1:MValue
        KArray(i,i-1) = c;
        KArray(i,i)= ad;
        KArray(i,i+1) = b;
    end
    disp(8);disp(KArray);disp(fArray);
    aArray = linsolve(KArray, -fArray);
    disp(aArray);
    phiArray = zeros(q, 1);
	for i = 1:1:length(points)
        phiArray(i) = phiValueOfFiniteElements(points(i), aArray, MValue);
	end;
endfunction

function run()
	alphaValue = evstr(get(alphaEdit, "String"));
	betaValue = evstr(get(betaEdit, "String"));
	MValue = evstr(get(MEdit, "String"));
	realSolveValue = evstr(get(realSolveCheckBox, "Value"));
	DirichletValue = evstr(get(DirichletRadioButton, "Value"));
    NeumannValue = evstr(get(NeumannRadioButton, "Value"));
    DirichletDiffValue = evstr(get(DirichletDiffRadioButton, "Value"));
    NeumannDiffValue = evstr(get(NeumannDiffRadioButton, "Value"));
	m1 = get(method1CheckBox, "Value");
	m2 = get(method2CheckBox, "Value");
	m3 = get(method3CheckBox, "Value");
	m4 = get(method4CheckBox, "Value");
	m5 = get(method5CheckBox, "Value");
	m6 = get(method6CheckBox, "Value");
    m7 = get(method7CheckBox, "Value");
    m8 = get(method8CheckBox, "Value");

	delete(gca());
	subplot(2, 1, 2);
	xgrid();

	if DirichletValue == 1 then

		if realSolveValue == 1 then
            //phA = DirichletODE2(alphaValue, betaValue, points);
            phA = analiticSolutionD();
            plot2d(points, phA, style = color(0, 0, 0));
		end;

		if m1 == 1 then
            phA = methodD1(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(255, 0, 0));
		end;

		if m2 == 1 then
            phA = methodD2(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 255, 0));
		end;

        if m3 == 1 then
            phA = methodD3(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 0, 255));
        end

        if m4 == 1 then
            phA = methodD4(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(255, 0, 255));
        end

        if m5 == 1 then
            phA = methodD5(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(255, 176, 0));
        end

        if m6 == 1 then
            phA = methodD6(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 176, 255));
        end

        if m8 == 1 then
            phA = finiteElementsD(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 179, 179));
        end

    elseif DirichletDiffValue == 1

        //realSolution = DirichletODE2(alphaValue, betaValue, points);
        realSolution = analiticSolutionD();

        if m1 == 1 then
            phA = methodD1(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 0, 0));
		end;

		if m2 == 1 then
            phA = methodD2(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 255, 0));
		end;

        if m3 == 1 then
            phA = methodD3(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 0, 255));
        end

        if m4 == 1 then
            phA = methodD4(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 0, 255));
        end

        if m5 == 1 then
            phA = methodD5(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 176, 0));
        end

        if m6 == 1 then
            phA = methodD6(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 176, 255));
        end

        if m8 == 1 then
            phA = finiteElementsD(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 179, 179));
        end

    elseif NeumannValue == 1

		if realSolveValue == 1 then
            phA = analiticSolutionN();
            plot2d(points, phA, style = color(0, 0, 0));
		end;

		if m1 == 1 then
            phA = methodN1(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(255, 0, 0));
		end;

        if m2 == 1 then
            phA = methodN2(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 255, 0));
		end;

        if m3 == 1 then
            phA = methodN3(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 0, 255));
		end;

        if m4 == 1 then
            phA = methodN4(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(255, 0, 255));
		end;

        if m5 == 1 then
            phA = methodN5(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(255, 176, 0));
		end;

        if m6 == 1 then
            phA = methodN6(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 176, 255));
		end;

        if m7 == 1 then
             phA = weakFormulation(alphaValue, betaValue, MValue);
             plot2d(points, phA, style = color(255, 128, 0));
		end;

        if m8 == 1 then
            phA = finiteElementsN(alphaValue, betaValue, MValue);
            plot2d(points, phA, style = color(0, 179, 179));
        end

     elseif NeumannDiffValue == 1

        //realSolution = NeumannODE2(alphaValue, betaValue, points);
        realSolution = analiticSolutionN();

        if m1 == 1 then
            phA = methodN1(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 0, 0));
		end;

		if m2 == 1 then
            phA = methodN2(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 255, 0));
		end;

        if m3 == 1 then
            phA = methodN3(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 0, 255));
        end

        if m4 == 1 then
            phA = methodN4(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 0, 255));
        end

        if m5 == 1 then
            phA = methodN5(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 176, 0));
        end

        if m6 == 1 then
            phA = methodN6(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 176, 255));
        end

        if m7 == 1 then
            phA = weakFormulation(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(255, 128, 0));
        end

        if m8 == 1 then
            phA = finiteElementsN(alphaValue, betaValue, MValue);
            phiArray = zeros(q, 1);
            for i = 1:1:length(points)
                phiArray(i) = abs(phA(i) - realSolution(i));
            end;
            plot2d(points, phiArray, style = color(0, 179, 179));
        end

	end;
endfunction;
